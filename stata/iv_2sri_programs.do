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
*   model1: homogeneous treatment effects — Z only, no Z*X
*   model2: heterogeneous treatment effects — Z + Z*X
*   no_iv:  no IV — same Xlist as model 1, but Z and Z*X terms removed
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
    if "`version'" == "model1" {
        * Model 1: homogeneous — Z only, no Z*X
        probit $Treated $Xlist_s1_m1 $Z, cluster($clustervar)
    }
    else if "`version'" == "model2" {
        * Model 2: heterogeneous — Z + Z*X
        probit $Treated $Xlist_s1_m2 $Z $Zlist_s1_m2, cluster($clustervar)
    }
    else if "`version'" == "no_iv" {
        probit $Treated $Xlist_s1_v1 $Z, cluster($clustervar)
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
    * Resolve Stage 2 Xlist based on version
    if "`version'" == "model1" {
        local s2_xlist $Xlist_s2_m1
    }
    else if "`version'" == "model2" {
        local s2_xlist $Xlist_s2_m2 $Dxlist_s2_m2
    }
    else if "`version'" == "no_iv" {
        local s2_xlist $Xlist_s2
    }
        else {
            noi disp "ERROR: unknown version `version'"
            exit 198
        }

    * Tobit: continuous outcomes (daoh, total_los) with horizon-specific limits
    if "$OutcomeFamily" == "gaussian" {
        local ul = $Horizon
        tobit $OutcometoUse `s2_xlist' $Treated $GenResidual, ///
            ll(0) ul(`ul') vce(cluster $clustervar)

        recycled_predictions, type(xb)    model(tobit)
        recycled_predictions, type(ystar) model(tobit) a(0) b(`ul')
        recycled_predictions, type(pr)    model(tobit) a(0) b(`ul')
        recycled_predictions, type(e)     model(tobit) a(0) b(`ul')
    }

    * Probit: binary outcomes (mortality, readmission)
    if "$OutcomeFamily" == "binomial" {
        probit $OutcometoUse `s2_xlist' $Treated $GenResidual, ///
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

    quietly gen _iv_depvar = `outcome'
	
	if "`version'" == "model1" {
    ivreg2 _iv_depvar $Xlist_s1_m1 ($Treated = $Z), ///
        cluster($clustervar) first ffirst
	}
	
	else if "`version'" == "model2" {
    ivreg2 _iv_depvar $Xlist_s1_m2 ($Treated = $Z $Zlist_s1_m2), ///
        cluster($clustervar) first ffirst
}

    weakivtest
    local eff_f  = r(F_eff)
    local crit5  = r(c_TSLS_5)
    local kp_f   = e(rkf)

    post iv_strength ///
        ("`outcome'") ("`version'") (`eff_f') (`crit5') (`kp_f')

    drop _iv_depvar

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
syntax, outcome(str) version(str) horizon(int) results_dir(str) dta_path(str)

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

    use "`dta_path'", clear
    encode nvrhospitalname, gen(HospitalCluster)
    gen All = 1

end
* --------------------------------------------------------------------
* Propensity score overlap plots
* For the specified version and outcome, runs the Stage 1 probit on the
* full sample to get predicted probabilities of treatment (propensity scores)
* Plots histograms of propensity scores by treatment group to visualize
* overlap. Only implemented for version 4 (LASSO-optimised) and 90-day outcomes
* --------------------------------------------------------------------
capture program drop plot_ps_overlap
program define plot_ps_overlap
    syntax, version(str) outcome(str) horizon(int) results_dir(str)

    * Run Stage 1 probit on full sample to get propensity scores
    if "`version'" == "z_only" {
        probit $Treated $Xlist_s1_v1 $Z, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage2_treatment" {
        probit $Treated $Xlist_s1_v2 $Z $Zlist_s1_v2, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage1_full" {
        probit $Treated $Xlist_s1_v3 $Z $Zlist_s1_v3, cluster($clustervar)
    }
    else if "`version'" == "z_x_stage1_instrument" {
        probit $Treated $Xlist_s1_v4 $Z $Zlist_s1_v4, cluster($clustervar)
    }

    tempvar ps
    predict `ps', pr

    local out_label = regexr("`outcome'", "_[0-9]+d$", "")

    twoway ///
        (hist `ps' if $Treated == 1, color(red%30) frequency) ///
        (hist `ps' if $Treated == 0, color(blue%60) frequency), ///
        legend(order(1 "Treated (early)" 2 "Control (late)")) ///
        xtitle("Predicted probability of early surgery") ///
        ytitle("Frequency") ///
        title("`out_label' | `version'") ///
        graphregion(color(white))

    graph export "`results_dir'ps_overlap_`version'.png", replace width(1600)

end
