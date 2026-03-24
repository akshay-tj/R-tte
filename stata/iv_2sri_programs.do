* iv_2sri_programs.do
* Reusable program definitions for 2SRI IV analysis with recycled predictions,
* subgroup bootstrapping, IV strength, and forest plots.
*
* Programs defined:
*   recycled_predictions  — counterfactual predictions under D=0 and D=1
*   SetOutcomeNumber      — derives outcome suffix from position in Outcomelist
*   subgroupboot          — 2SRI bootstrap, configurable Stage 1 version
*   compute_iv_strength   — Montiel Olea-Pflueger effective F (clustered)
*   forest_plot           — forest plot of treatment effects by version/outcome
*   extract_bootstrap_results — extracts point estimates and CIs from bootstrap DTA
*
* NOTE: weakivtest implements Montiel Olea-Pflueger effective F statistic.
* Designed for linear IV (2SLS) — used here as a diagnostic approximation
* for the probit first stage. Interpret with caution.
* Install: ssc install weakivtest
* Fallback: net install weakivtest, from("https://raw.githubusercontent.com/Alalalalaki/weakivtest/master/") replace
* Also required: ssc install labutil (for forest plot y-axis labels)

* --------------------------------------------------------------------
* Recycled predictions
* --------------------------------------------------------------------
capture program drop recycled_predictions
program define recycled_predictions, rclass
syntax, model(str) type(str) ///
    [a(numlist min=1 max=1 missingokay)] ///
    [b(numlist min=1 max=1 missingokay)]

    if inlist("`type'", "ystar", "e") & "`a'" == "" {
        noi disp "parameters a and b must be specified for type `type'"
        exit
    }
    if inlist("`type'", "ystar", "e") & "`b'" == "" {
        noi disp "parameters a and b must be specified for type `type'"
        exit
    }
    if "`type'" == "pr" & "`model'" == "tobit" & "`a'" == "" {
        noi disp "parameters a and b must be specified for tobit pr"
        exit
    }
    if "`type'" == "pr" & "`model'" == "tobit" & "`b'" == "" {
        noi disp "parameters a and b must be specified for tobit pr"
        exit
    }

    if "`type'" == "xb"                          local opt "xb"
    if "`type'" == "ystar"                       local opt "ystar(`a',`b')"
    if "`type'" == "pr" & "`model'" == "tobit"   local opt "pr(`a',`b')"
    if "`type'" == "pr" & "`model'" == "probit"  local opt "pr"
    if "`type'" == "e"                           local opt "e(`a',`b')"

    * Store original treatment and D*X interaction values
    gen orig_treated = $Treated
    foreach var of varlist d_x_* {
        gen orig_`var' = `var'
    }

    * Counterfactual: treatment = 0
    * D*X cols must be updated to reflect flipped treatment
    replace $Treated = 0
    foreach var of varlist d_x_* {
        quietly replace `var' = 0
    }
    predict y0`type'$OutcomeNumber, `opt'

    * Counterfactual: treatment = 1
    replace $Treated = 1
    foreach var of varlist d_x_* {
        local base = subinstr("`var'", "d_x_", "", 1)
        quietly replace `var' = `base'
    }
    predict y1`type'$OutcomeNumber, `opt'

    gen ite_`type'$OutcomeNumber = y1`type'$OutcomeNumber - y0`type'$OutcomeNumber

    * Restore original treatment and D*X values
    replace $Treated = orig_treated
    foreach var of varlist d_x_* {
        quietly replace `var' = orig_`var'
        drop orig_`var'
    }
    drop orig_treated
end

* --------------------------------------------------------------------
* Derive outcome number suffix from outcome name programmatically
* --------------------------------------------------------------------
capture program drop SetOutcomeNumber
program define SetOutcomeNumber
    * Format: OutN where N is position in $Outcomelist
    global OutcomeNumber ""
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
* D*X and Z*X interactions are pre-generated in R and included directly
* in the appropriate Xlist global — no factor variable syntax needed.
*
* version options (aligned with R globals CSV):
*   z_only                — Stage 1: Xlist + Z only
*   z_x_stage2_treatment  — Stage 1: Xlist + Z + Z*Mlist cols (pre-generated)
*   z_x_stage1_full       — Stage 1: Xlist + Z + Z*all-main-effects cols
*   z_x_stage1_instrument — Stage 1: Xlist_stage1 + Z + Z*Zlist cols (LASSO-optimised)
*
* Each version reads its Xlist from the corresponding global set in the
* experiment script from the globals CSV.
* --------------------------------------------------------------------
capture program drop subgroupboot
program define subgroupboot, rclass
syntax, ytouse(str) version(str)

    global OutcometoUse `ytouse'
    SetOutcomeNumber

    * ── Stage 1: exposure model ──────────────────────────────────────
    * Each version uses a different Stage 1 Xlist (main effects + X*X +
    * version-specific Z*X cols). D*X cols are NOT in Stage 1.
    * Z is added separately via $Z.
    if "`version'" == "z_only" {
        probit $Treated $Xlist_s1_v1 $Z, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage2_treatment" {
        probit $Treated $Xlist_s1_v2 $Z, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage1_full" {
        probit $Treated $Xlist_s1_v3 $Z, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage1_instrument" {
        probit $Treated $Xlist_s1_v4 $Z, cluster($clustervar)
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
    * Same across all versions — uses Xlist_s2 which includes
    * main effects + X*X outcome terms + pre-generated D*X cols.
    * Tobit: continuous outcomes (daoh, total_los) with horizon-specific limits
    if "$OutcomeFamily" == "gaussian" {
        local ul = $Horizon
        tobit $OutcometoUse $Xlist_s2 $Treated $GenResidual, ///
            ll(0) ul(`ul') vce(cluster $clustervar)

        recycled_predictions, type(xb)    model(tobit)
        recycled_predictions, type(ystar) model(tobit) a(0) b(`ul')
        recycled_predictions, type(pr)    model(tobit) a(0) b(`ul')
        recycled_predictions, type(e)     model(tobit) a(0) b(`ul')
    }

    * Probit: binary outcomes (mortality, readmission)
    if "$OutcomeFamily" == "binomial" {
        probit $OutcometoUse $Xlist_s2 $Treated $GenResidual, ///
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
* Uses linear first stage via ivreg2 (weakivtest post-estimation command
* requires ivreg2, not probit). Kleibergen-Paap rk Wald F also reported.
*
* Versions 1-3: run once per version (Xlist same across outcomes)
* Version 4:    run once per outcome x version (Xlist differs per outcome)
*
* Requires: ssc install ivreg2
*           ssc install weakivtest (or net install from GitHub)
* --------------------------------------------------------------------
capture program drop compute_iv_strength

program define compute_iv_strength

    syntax, version(str) outcome(str)

    * Use short alias to avoid temp variable name length issues in ivreg2
    quietly gen _iv_depvar = `outcome'

    if "`version'" == "z_only" {

        ivreg2 _iv_depvar $Xlist_s1_v1 ($Treated = $Z), ///
            cluster($clustervar) first

    }
    else if "`version'" == "z_x_stage2_treatment" {

        ivreg2 _iv_depvar $Xlist_s1_v2 ($Treated = $Z), ///
            cluster($clustervar) first

    }
    else if "`version'" == "z_x_stage1_full" {

        ivreg2 _iv_depvar $Xlist_s1_v3 ($Treated = $Z), ///
            cluster($clustervar) first

    }
    else if "`version'" == "z_x_stage1_instrument" {

        ivreg2 _iv_depvar $Xlist_s1_v4 ($Treated = $Z), ///
            cluster($clustervar) first

    }

    * Montiel Olea-Pflueger effective F (post ivreg2)
    weakivtest
    local eff_f  = r(eff_F)
    local crit5  = r(crit5)

    * Kleibergen-Paap rk Wald F from ivreg2
    local kp_f   = e(rkf)

    post iv_strength ///
        ("`outcome'") ("`version'") (`eff_f') (`crit5') (`kp_f')

    * Clean up
    drop _iv_depvar

end
* --------------------------------------------------------------------
* Forest plot — one per timepoint
* Reads bootstrap result CSVs produced by extract_bootstrap_results.
* Effect measure: ystar for tobit (gaussian), pr for probit (binomial)
* Family is read from the globals CSV per outcome — not from $OutcomeFamily
* global which reflects only the last outcome set in the loop.
* Versions on y-axis, outcomes within version, CI from bootstrap.
* --------------------------------------------------------------------
capture program drop forest_plot
program define forest_plot
syntax, horizon(int) results_dir(str) data_dir(str)

    local versions z_only z_x_stage2_treatment z_x_stage1_full z_x_stage1_instrument

    * Read family per outcome from globals CSV
    preserve
    local globals_path = "`data_dir'" + "non_elective_`horizon'd_globals.csv"
    import delimited using "`globals_path'", clear varnames(1)
    local n_outcomes = _N
    forvalues i = 1/`n_outcomes' {
        local fam_`i'     = family[`i']
        local outcome_`i' = outcome[`i']
    }
    restore

    * Collect results from bootstrap CSVs into one tempfile
    tempfile plotdata
    tempname fh
    file open `fh' using `plotdata', write replace
    file write `fh' "outcome,version,est,ci_lo,ci_hi,label" _n

    foreach v of local versions {
        forvalues i = 1/`n_outcomes' {
            local o   = "`outcome_`i''"
            local fam = "`fam_`i''"

            * Select pred type from family — not from global
            local pred_prefix = cond("`fam'" == "gaussian", "mystar", "mpr")

            * Shorten outcome name for label (strip horizon suffix e.g. _90d)
            local o_short = regexr("`o'", "_[0-9]+d$", "")

            local csvfile = "`results_dir'bootstrap_results_`o'_`v'_`horizon'd.csv"
            capture {
                preserve
                import delimited using "`csvfile'", clear varnames(1)

                * Find the All subgroup ITE row for correct pred type
                quietly keep if strpos(stat, "`pred_prefix'") & strpos(stat, "All") & !strpos(stat, "m0") & !strpos(stat, "m1")
                if _N > 0 {
                    local est = est[1]
                    local lo  = ci_lo[1]
                    local hi  = ci_hi[1]
                    file write `fh' "`o',`v',`est',`lo',`hi',`o_short' | `v'" _n
                }
                restore
            }
        }
    }
    file close `fh'

    * Import and plot
    preserve
    import delimited using `plotdata', clear

    gen row = _n

    * Build value label from label column
    tostring row, gen(rowstr)
    labmask row, values(label)

    twoway ///
        (rcap ci_hi ci_lo row, horizontal lcolor(gs8)) ///
        (scatter row est, mcolor(navy) msymbol(circle)) ///
        , ylabel(1/`=_N', valuelabel angle(0)) ///
          xline(0, lcolor(red) lpattern(dash)) ///
          xlabel(, format(%6.2f)) ///
          ytitle("") xtitle("Treatment effect estimate") ///
          title("Clinical effectiveness — `horizon'-day outcomes") ///
          legend(off) ///
          graphregion(color(white))

    graph export "`results_dir'forest_plot_`horizon'd.png", replace width(2400)
    restore
end

* --------------------------------------------------------------------
* Extract bootstrap results from saved DTA
* Point estimate: mean of bootstrap replications
* SE: SD of bootstrap replications
* CI: normal-based (mean ± 1.96 * SE)
* NOTE: point estimate is mean of replications, not full-sample estimate.
*       Full-sample estimate to be added in a future update.
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
    file write `fh' "outcome,version,horizon,stat,est,se,ci_lo,ci_hi" _n

    foreach var of varlist _all {
        quietly summarize `var'
        local est = r(mean)
        local se  = r(sd)
        local lo  = `est' - 1.96 * `se'
        local hi  = `est' + 1.96 * `se'

        file write `fh' ///
            "`outcome',`version',`horizon',`var',`est',`se',`lo',`hi'" _n
    }

    file close `fh'
    restore
end
