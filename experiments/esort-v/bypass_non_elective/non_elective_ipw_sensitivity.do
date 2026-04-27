* ipw_sensitivity.do
*
* IPW-weighted sensitivity analysis (no IV) for non-elective bypass cohort.
* Uses LASSO-selected covariate lists from globals CSV.
* Doubly-robust: IPW weight + direct covariate adjustment in outcome model.
* Analytic CIs via margins (clustered sandwich SE).
* One log file per horizon.
*
* Inputs:
*   - non_elective_{H}d.dta
*   - non_elective_{H}d_globals.csv
*
* Output:
*   - ipw_sensitivity_{H}d.smcl (log, one per horizon)

* ============================================================================
* CONFIGURATION
* ============================================================================

local  horizons    90 180 365

global data_dir    "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/lasso_outputs/"
global results_dir "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_non_elective_240426/clinical_effectiveness_results/"

global treated     early_surgery
global clustervar  HospitalCluster

* Year range for surgyr reference category derivation
local  surgyr_min  2015
local  surgyr_max  2023

* IPW weight diagnostic percentiles to report
local  diag_pctiles 95 99 99.5 99.9

* ============================================================================
* SETUP
* ============================================================================

capture mkdir "$results_dir"

* ============================================================================
* MAIN LOOP — horizons
* ============================================================================

foreach horizon of local horizons {

    local dta_path     "$data_dir/non_elective_`horizon'd.dta"
    local globals_path "$data_dir/non_elective_`horizon'd_globals.csv"

    * ── Log ──────────────────────────────────────────────────────────
    capture log close
    log using "$results_dir/ipw_sensitivity_`horizon'd.smcl", replace

    * ── Load data ────────────────────────────────────────────────────
    use "`dta_path'", clear

    sort nvrhospitalname
    encode nvrhospitalname, gen(HospitalCluster)

    * Convenience indicator
    gen byte All = 1

    * Binary complements
    gen byte gender_M         = 1 - gender_F
    gen byte fontaine_not4    = 1 - fontaine_4
    gen byte comorbidity_not1 = 1 - comorbidity_1
    gen byte krt_no           = 1 - krt_yn

    * Reference year: all other surgyr dummies == 0
    local surgyr_conditions ""
    forvalues y = `=`surgyr_min'+1'/`surgyr_max' {
        local surgyr_conditions "`surgyr_conditions' & surgyr_`y' == 0"
    }
    gen byte surgyr_`surgyr_min' = (1 == 1 `surgyr_conditions')

    * ── Read globals CSV ─────────────────────────────────────────────
    * PS covariates are outcome-invariant (model1), so row 1 suffices.
    preserve
        import delimited using "`globals_path'", clear varnames(1)
        local n_outcomes = _N
        forvalues i = 1/`n_outcomes' {
            local outcome_`i'       = outcome[`i']
            local family_`i'        = family[`i']
            local xlist_s2_noiv_`i' = xlist_s2_noiv[`i']
        }
        global Xlist_ps = xlist_s1_m1[1]
    restore

    * ── Propensity score model ───────────────────────────────────────
    display _newline "===== Propensity score model | `horizon'd ====="
    probit $treated $Xlist_ps, vce(cluster $clustervar)
    predict double ps, pr

    * ── IPW weights ──────────────────────────────────────────────────
    gen double ipw = cond($treated == 1, 1/ps, 1/(1-ps))

    * ── Weight diagnostics ───────────────────────────────────────────
    display _newline "===== IPW weight diagnostics | `horizon'd ====="

    summ ipw, detail
    by $treated, sort: summ ipw

    * Tail percentiles
    local plist : subinstr local diag_pctiles " " ",", all
    _pctile ipw, p(`plist')
    local j = 1
    foreach p of local diag_pctiles {
        display "  `p'th percentile: " r(r`j')
        local ++j
    }

    * Effective sample size (Kish approximation)
    quietly {
        summ ipw
        local n     = r(N)
        local wmean = r(mean)
        gen double ipw2 = ipw^2
        summ ipw2
        local ess = (`n' * `wmean')^2 / (r(sum))
        drop ipw2
    }
    display "  ESS (Kish): " %9.1f `ess' "  (N = `n')"

    * ── Outcome loop ─────────────────────────────────────────────────
    display _newline "===== IPW-weighted ATE | `horizon'd ====="

    forvalues i = 1/`n_outcomes' {

        local outcome  "`outcome_`i''"
        local family   "`family_`i''"
        local s2_xlist "`xlist_s2_noiv_`i''"

        display _newline "--- Outcome: `outcome' ---"

        if "`family'" == "gaussian" {

            tobit `outcome' `s2_xlist' $treated [pweight = ipw], ///
                ll(0) ul(`horizon') vce(cluster $clustervar)

            margins, dydx($treated)

            matrix b  = r(b)
            matrix V  = r(V)
            scalar ate = b[1,1]
            scalar se  = sqrt(V[1,1])
            display "ATE (early - later):       " %6.3f ate ///
                    "  95% CI: (" %6.3f ate - 1.96*se ", " %6.3f ate + 1.96*se ")"
        }

        if "`family'" == "binomial" {

            probit `outcome' `s2_xlist' $treated [pweight = ipw], ///
                vce(cluster $clustervar)

            margins, dydx($treated) predict(pr)

            matrix b  = r(b)
            matrix V  = r(V)
            scalar ate = b[1,1]
            scalar se  = sqrt(V[1,1])
            display "Risk diff (early - later): " %6.3f ate ///
                    "  95% CI: (" %6.3f ate - 1.96*se ", " %6.3f ate + 1.96*se ")"
        }
    }

    display _newline "===== End | `horizon'd ====="
    log close
}