********************************************************************************
* elective_sens_trimmed_window_clinical_effectiveness.do
*
* Sensitivity analysis for elective bypass cohort.
*
* Sensitivity trimmed window treatment definition:
*   early_surg_sens = 1 if daystosurgery is 1–14 days inclusive
*   early_surg_sens = 0 if daystosurgery is 15–28 days inclusive
*
* This script:
*   1. Reads sensitivity .dta files produced by elective_sens_trimmed_window_setup.R
*   2. Reads outcome-specific covariate lists from the matching globals CSV
*   3. Calls programs defined in stata/iv_2sri_programs.do
*   4. Runs bootstrapped 2SRI models for all outcomes x time horizons
********************************************************************************

* ============================================================================
* CONFIGURATION — edit here only
* ============================================================================

* Timepoints to analyse; must match files produced by R
local horizons 90 180 365

* Bootstrap settings
global nreps = 300
global seed = 37563845

* Versions to run
* model1: homogeneous effects
* model2: heterogeneous effects
* no_iv: non-IV regression
global versions model2
* global versions no_iv model1 model2

* Paths
global data_dir "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/sensitivity_trimmed_window/lasso_outputs/"
global results_dir "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/sensitivity_trimmed_window/clinical_effectiveness_results/"
global programs_dir "C:/Users/LSHAJ82/Documents/GitHub/R-tte/stata/"

* Subgroups for recycled predictions
global Subgrouplist All gender_F gender_M fontaine_4 fontaine_not4 comorbidity_1 comorbidity_not1 krt_yn krt_no surgyr_2016 surgyr_2017 surgyr_2018 surgyr_2019 surgyr_2020 surgyr_2021 surgyr_2022 surgyr_2023

* Cluster variable
global clustervar HospitalCluster

* Treatment and instrument
global Treated early_surg_sens
global Z instrumental_variable

* ============================================================================
* SETUP
* ============================================================================

capture mkdir "$results_dir"
cd "$results_dir"

do "${programs_dir}iv_2sri_programs.do"

* Initialise IV strength postfile
tempfile iv_strength_file
capture postclose iv_strength

postfile iv_strength str50 outcome str50 version ///
    double eff_f double crit5 double kp_f ///
    using `iv_strength_file', replace

* ============================================================================
* MAIN LOOP — time horizons
* ============================================================================

foreach horizon of local horizons {

    global Horizon = `horizon'

    local dta_path = "$data_dir" + "elective_sens_`horizon'd.dta"
    local globals_path = "$data_dir" + "elective_sens_`horizon'd_globals.csv"

    * ------------------------------------------------------------------------
    * Load data and create derived subgroup / cluster variables
    * ------------------------------------------------------------------------

    use "`dta_path'", clear

    encode nvrhospitalname, gen(HospitalCluster)

    gen All = 1
    gen gender_M = 1 - gender_F
    gen fontaine_not4 = 1 - fontaine_4
    gen comorbidity_not1 = 1 - comorbidity_1
    gen krt_no = 1 - krt_yn

    gen surgyr_2016 = ///
        (surgyr_2017 == 0 & ///
         surgyr_2018 == 0 & ///
         surgyr_2019 == 0 & ///
         surgyr_2020 == 0 & ///
         surgyr_2021 == 0 & ///
         surgyr_2022 == 0 & ///
         surgyr_2023 == 0)

    * Optional quick check: treatment group counts in analysis data
    di "============================================================"
    di "Sensitivity cohort counts, horizon `horizon'd"
    tab early_surg_sens
    di "============================================================"

    * ------------------------------------------------------------------------
    * Read globals CSV to get outcome names and families
    * ------------------------------------------------------------------------

    preserve
        import delimited using "`globals_path'", clear varnames(1)

        local n_outcomes = _N

        forvalues i = 1/`n_outcomes' {
            local outcome_`i' = outcome[`i']
            local family_`i' = family[`i']
        }
    restore

    * Build global outcome list
    global Outcomelist

    forvalues i = 1/`n_outcomes' {
        global Outcomelist $Outcomelist `outcome_`i''
    }

    * ------------------------------------------------------------------------
    * Loop over outcomes
    * ------------------------------------------------------------------------

    forvalues i = 1/`n_outcomes' {

        local outcome = "`outcome_`i''"
        local family = "`family_`i''"

        * Re-read globals CSV for this outcome.
        * Direct global assignment avoids Stata local macro length limits.
        preserve
            import delimited using "`globals_path'", clear varnames(1)

            global Xlist_s2_m1 = xlist_s2_m1[`i']
            global Xlist_s1_m1 = xlist_s1_m1[`i']

            global Xlist_s2_m2 = xlist_s2_m2[`i']
            global Xlist_s1_m2 = xlist_s1_m2[`i']
            global Zlist_s1_m2 = zlist_s1_m2[`i']
            global Dxlist_s2_m2 = dxlist_s2_m2[`i']

            global Xlist_s2_noiv = xlist_s2_noiv[`i']
        restore

        global OutcomeFamily `family'

        * --------------------------------------------------------------------
        * IV strength diagnostics
        * --------------------------------------------------------------------

        compute_iv_strength, version("model1") outcome("`outcome'")
        compute_iv_strength, version("model2") outcome("`outcome'")

        * --------------------------------------------------------------------
        * Propensity score overlap
        * --------------------------------------------------------------------

        plot_ps_overlap, version("model2") outcome("`outcome'") ///
            horizon(`horizon') results_dir("$results_dir")

        * --------------------------------------------------------------------
        * Bootstrap loop over selected model versions
        * --------------------------------------------------------------------

        foreach version of global versions {

            global OutcometoUse `outcome'
            SetOutcomeNumber

            if "$OutcomeFamily" == "gaussian" local predlist xb ystar pr e
            if "$OutcomeFamily" == "binomial" local predlist pr

            local bootstats

            foreach pred of local predlist {
                foreach sub of global Subgrouplist {

                    local bootstats `bootstats' ///
                        m0`pred'_`sub'$OutcomeNumber = ///
                        r(mean0`pred'_`sub'$OutcomeNumber)

                    local bootstats `bootstats' ///
                        m1`pred'_`sub'$OutcomeNumber = ///
                        r(mean1`pred'_`sub'$OutcomeNumber)

                    local bootstats `bootstats' ///
                        m`pred'_`sub'$OutcomeNumber = ///
                        r(mean`pred'_`sub'$OutcomeNumber)
                }
            }

            capture log close

            log using "bootstrap_`outcome'_`version'_`horizon'd.smcl", replace

            bootstrap `bootstats', ///
                reps($nreps) seed($seed) ///
                cluster($clustervar) idcluster(id_$clustervar) ///
                saving("bootstrap_`outcome'_`version'_`horizon'd.dta", replace): ///
                subgroupboot, ytouse("`outcome'") version("`version'")

            log close

            * Extract point estimates and CIs from bootstrap DTA
            extract_bootstrap_results, ///
                outcome("`outcome'") ///
                version("`version'") ///
                horizon(`horizon') ///
                results_dir("$results_dir") ///
                dta_path("`dta_path'")

            * Reload clean analysis data before next model
            use "`dta_path'", clear

            encode nvrhospitalname, gen(HospitalCluster)

            gen All = 1
            gen gender_M = 1 - gender_F
            gen fontaine_not4 = 1 - fontaine_4
            gen comorbidity_not1 = 1 - comorbidity_1
            gen krt_no = 1 - krt_yn

            gen surgyr_2016 = ///
                (surgyr_2017 == 0 & ///
                 surgyr_2018 == 0 & ///
                 surgyr_2019 == 0 & ///
                 surgyr_2020 == 0 & ///
                 surgyr_2021 == 0 & ///
                 surgyr_2022 == 0 & ///
                 surgyr_2023 == 0)
        }
    }
}

* ============================================================================
* SAVE IV STRENGTH RESULTS
* ============================================================================

postclose iv_strength

use `iv_strength_file', clear

tostring outcome version, replace force

export delimited using "${results_dir}iv_strength_all.csv", replace

di "Sensitivity 2SRI clinical effectiveness analysis complete."