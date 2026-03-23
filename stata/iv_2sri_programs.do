* --------------------------------------------------------------------
* NOTE: weakivtest implements Montiel Olea-Pflueger effective F statistic.
* Designed for linear IV (2SLS) — used here as a diagnostic approximation
* for the probit first stage. 
* --------------------------------------------------------------------
* iv_2sri_programs.do
* Reusable program definitions for 2SRI IV analysis with recycled predictions,
* subgroup bootstrapping, IV strength, and forest plots.
*
* Programs defined:
*   recycled_predictions  — counterfactual predictions under D=0 and D=1
*   subgroupboot          — 2SRI bootstrap, configurable Stage 1 version
*   compute_iv_strength   — Montiel Olea-Pflueger effective F (clustered)
*   forest_plot           — forest plot of treatment effects by version/outcome

* --------------------------------------------------------------------
* Recycled predictions
* --------------------------------------------------------------------
capture program drop recycled_predictions
program define recycled_predictions, rclass
syntax, model(str) type(str) ///
    [a(numlist min=1 max=1 missingokay)] ///
    [b(numlist min=1 max=1 missingokay)]

    if inlist("`type'", "ystar", "pr", "e") & "`a'" == "" {
        noi disp "parameters a and b must be specified for type `type'"
        exit
    }
    if inlist("`type'", "ystar", "pr", "e") & "`b'" == "" {
        noi disp "parameters a and b must be specified for type `type'"
        exit
    }

    if "`type'" == "xb"                          local opt "xb"
    if "`type'" == "ystar"                       local opt "ystar(`a',`b')"
    if "`type'" == "pr" & "`model'" == "tobit"   local opt "pr(`a',`b')"
    if "`type'" == "pr" & "`model'" == "probit"  local opt "pr"
    if "`type'" == "e"                           local opt "e(`a',`b')"

    gen orig_treated = $Treated
    replace $Treated = 0
    predict y0`type'$OutcomeNumber, `opt'
    replace $Treated = 1
    predict y1`type'$OutcomeNumber, `opt'
    gen ite_`type'$OutcomeNumber = y1`type'$OutcomeNumber - y0`type'$OutcomeNumber
    replace $Treated = orig_treated
    drop orig_treated
end

* --------------------------------------------------------------------
* Derive outcome number suffix from outcome name programmatically
* --------------------------------------------------------------------
capture program drop SetOutcomeNumber
program define SetOutcomeNumber
    * Extracts a short numeric suffix from outcome name for variable naming.
    * Format: OutN where N is position in $Outcomelist
    local i = 1
    foreach o of global Outcomelist {
        if "$OutcometoUse" == "`o'" {
            global OutcomeNumber = "Out`i'"
        }
        local i = `i' + 1
    }
    if "$OutcomeNumber" == "" {
        noi disp "WARNING: outcome $OutcometoUse not found in Outcomelist"
    }
end

* --------------------------------------------------------------------
* 2SRI bootstrap program — configurable Stage 1 version
*
* version options (aligned with R globals CSV):
*   z_only                — Stage 1: Xlist + Z only
*   z_x_stage2_treatment  — Stage 1: Xlist + Z + Z*stage2_treatment_interactions
*   z_x_stage1_full       — Stage 1: Xlist_stage1 + Z + Z*stage2_main_effects_and_xx
*   z_x_stage1_instrument — Stage 1: Xlist_stage1 + Z + Z*stage1_instrument_interactions
* --------------------------------------------------------------------
capture program drop subgroupboot
program define subgroupboot, rclass
syntax, ytouse(str) version(str)

    global OutcometoUse `ytouse'
    SetOutcomeNumber

    * ── Stage 1: exposure model ──────────────────────────────────────
    if "`version'" == "z_only" {
        probit $Treated $Xlist_stage2 $Z, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage2_treatment" {
        probit $Treated $Xlist_stage2 $Z c.$Z#($Mlist_bin) c.$Z#c.ageatsurgery, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage1_full" {
        probit $Treated $Xlist_stage2 $Z c.$Z#($Xlist_stage2_bin) c.$Z#c.ageatsurgery c.$Z#c.age_sq, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage1_instrument" {
        probit $Treated $Xlist_stage1 $Z c.$Z#($Zlist_bin) c.$Z#c.ageatsurgery, cluster($clustervar)
    }
    else {
        noi disp "ERROR: unknown version `version'"
        exit 198
    }

    capture drop gen_resid
    predict gen_resid, score
    if $withIV == 1 {
        global GenResidual gen_resid
    }
    else {
        global GenResidual
    }

    * ── Stage 2: outcome model ───────────────────────────────────────
    * Tobit: continuous outcomes (daoh, total_los) with horizon-specific limits
    if "$OutcomeFamily" == "gaussian" {
        local ul = $Horizon
        tobit $OutcometoUse $Xlist_stage2 $Treated ///
            i.$Treated#($Mlist_bin) c.$Treated#c.ageatsurgery $GenResidual, ///
            ll(0) ul(`ul') vce(cluster $clustervar)

        recycled_predictions, type(xb)    model(tobit)
        recycled_predictions, type(ystar) model(tobit) a(0) b(`ul')
        recycled_predictions, type(pr)    model(tobit) a(0) b(`ul')
        recycled_predictions, type(e)     model(tobit) a(0) b(`ul')
    }

    * Probit: binary outcomes (mortality, readmission)
    if "$OutcomeFamily" == "binomial" {
        probit $OutcometoUse $Xlist_stage2 $Treated ///
            i.$Treated#($Mlist_bin) c.$Treated#c.ageatsurgery $GenResidual, ///
            cluster($clustervar)

        recycled_predictions, type(pr) model(probit)
    }

    * ── Collect subgroup means ───────────────────────────────────────
    foreach pred in xb ystar e pr {
        foreach sub of global Subgrouplist {
            quietly capture summarize y0`pred'$OutcomeNumber if `sub' == 1
            if _rc == 0 return scalar mean0`pred'_`sub'$OutcomeNumber = r(mean)

            quietly capture summarize y1`pred'$OutcomeNumber if `sub' == 1
            if _rc == 0 return scalar mean1`pred'_`sub'$OutcomeNumber = r(mean)

            quietly capture summarize ite_`pred'$OutcomeNumber if `sub' == 1
            if _rc == 0 return scalar mean`pred'_`sub'$OutcomeNumber = r(mean)
        }
        capture drop ite_`pred'$OutcomeNumber
        capture drop y0`pred'$OutcomeNumber
        capture drop y1`pred'$OutcomeNumber
    }
end

* --------------------------------------------------------------------
* IV strength — Montiel Olea-Pflueger effective F (clustered)
* Run once per outcome x version (Stage 1 differs across both)
* Requires -weakivtest- package: net install weakivtest
* --------------------------------------------------------------------
capture program drop compute_iv_strength
program define compute_iv_strength
syntax, ytouse(str) version(str)

    global OutcometoUse `ytouse'
    SetOutcomeNumber

    * Run Stage 1 matching the version
    if "`version'" == "z_only" {
        probit $Treated $Xlist_stage2 $Z, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage2_treatment" {
        probit $Treated $Xlist_stage2 $Z c.$Z#($Mlist_bin) c.$Z#c.ageatsurgery, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage1_full" {
        probit $Treated $Xlist_stage2 $Z c.$Z#($Xlist_stage2) c.$Z#c.ageatsurgery, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage1_instrument" {
        probit $Treated $Xlist_stage1 $Z c.$Z#($Zlist_bin) c.$Z#c.ageatsurgery, cluster($clustervar)
    }

    * Montiel Olea-Pflueger effective F
    weakivtest

    * Export to CSV
    local eff_f  = r(eff_F)
    local crit5  = r(crit5)
    post iv_strength ///
        ("$OutcometoUse") ("`version'") (`eff_f') (`crit5')
end

* --------------------------------------------------------------------
* Forest plot — one per timepoint
* Versions on y-axis, outcomes within version, CI from bootstrap
* Effect measure: ystar for tobit, pr for probit
* --------------------------------------------------------------------
capture program drop forest_plot
program define forest_plot
syntax, horizon(int)

    * Collect point estimates and CIs from bootstrap results
    * Expects matrices: forest_est, forest_lo, forest_hi, forest_labels
    * populated in the experiment script after bootstrap

    local versions z_only z_x_stage2_treatment z_x_stage1_full z_x_stage1_instrument
    local outcomes $Outcomelist
    local nrow = 0

    * Count rows
    foreach v of local versions {
        foreach o of local outcomes {
            local nrow = `nrow' + 1
        }
    }

    * Build plot data from stored bootstrap scalars
    tempfile plotdata
    tempname fh
    file open `fh' using `plotdata', write replace
    file write `fh' "outcome,version,est,lo,hi" _n

    foreach v of local versions {
        foreach o of local outcomes {
            global OutcometoUse `o'
            SetOutcomeNumber

            * Select pred type based on family
            local pred = cond("$OutcomeFamily" == "gaussian", "ystar", "pr")

            capture {
                local est = _b[m`pred'_All$OutcomeNumber]
                local lo  = _b[m`pred'_All$OutcomeNumber] - ///
                    1.96 * _se[m`pred'_All$OutcomeNumber]
                local hi  = _b[m`pred'_All$OutcomeNumber] + ///
                    1.96 * _se[m`pred'_All$OutcomeNumber]
                file write `fh' "`o',`v',`est',`lo',`hi'" _n
            }
        }
    }
    file close `fh'

    * Import and plot
    preserve
    import delimited using `plotdata', clear

    gen row = _n
    gen zero_line = 0

    twoway ///
        (rcap hi lo row, horizontal lcolor(gs8)) ///
        (scatter row est, mcolor(navy) msymbol(circle)) ///
        (line zero_line row, lcolor(red) lpattern(dash)) ///
        , ylabel(1/`nrow', valuelabel angle(0)) ///
          xlabel(, format(%6.2f)) ///
          ytitle("") xtitle("Treatment effect estimate") ///
          title("Clinical effectiveness — `horizon'-day outcomes") ///
          legend(off) ///
          graphregion(color(white))

    graph export "forest_plot_`horizon'd.png", replace width(2400)
    restore
end

* --------------------------------------------------------------------
* Extract bootstrap results from saved DTA
* Point estimate: from _b_ prefixed variables (original sample run)
* CI: 2.5th / 97.5th percentiles of replication columns
* --------------------------------------------------------------------
capture program drop extract_bootstrap_results
program define extract_bootstrap_results
syntax, outcome(str) version(str) horizon(int) results_dir(str)

    preserve
    use "`results_dir'bootstrap_`outcome'_`version'_`horizon'd.dta", clear

    tempname fh
    file open `fh' using ///
        "`results_dir'bootstrap_results_`outcome'_`version'_`horizon'd.csv", ///
        write replace
    file write `fh' "outcome,version,horizon,stat,est,ci_lo,ci_hi" _n

    foreach var of varlist _b_* {
        * Stat name: strip _b_ prefix to get replication column name
        local stat = subinstr("`var'", "_b_", "", 1)

        * Point estimate from original sample run
        local est = `var'[1]

        * CI from percentiles of replication columns
        quietly centile `stat', centile(2.5 97.5)
        local lo = r(c_1)
        local hi = r(c_2)

        file write `fh' ///
            "`outcome',`version',`horizon',`stat',`est',`lo',`hi'" _n
    }

    file close `fh'
    restore
end