* elective_run_clinical_effectiveness.do
*
* Runs 2SRI IV clinical effectiveness analysis for elective bypass cohort.
* Reads LASSO-selected variable lists from globals CSV (produced by R).
* Runs bootstrapped recycled predictions for all outcomes x timepoints.
* Produces IV strength diagnostics.
*
* D*X and Z*X interaction columns are pre-generated in R and included
* directly in the appropriate Xlist — no factor variable syntax in Stata.

* ============================================================================
* CONFIGURATION — edit here only
* ============================================================================

* Timepoints to analyse (must match DTA/globals files produced by R)
local horizons 90 180 365 //   

* Bootstrap settings
global nreps  = 300       // number of bootstrap replications
global seed   = 37563845  // bootstrap seed 

* Versions to run (model1: homogeneous, model2: heterogeneous, no_iv: no instrument)
global versions model2 // no_iv model1

* Paths
global data_dir     "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/lasso_outputs/"
global results_dir  "Z:/PHP/HSR/ESORT-V/ESORT-V/bypass_elective_270426/clinical_effectiveness_results/"
global programs_dir "C:/Users/LSHAJ82/Documents/GitHub/R-tte/stata/"

* Subgroups for recycled predictions
global Subgrouplist All gender_F gender_M fontaine_4 fontaine_not4 comorbidity_1 comorbidity_not1 krt_yn krt_no surgyr_2016 surgyr_2017 surgyr_2018 surgyr_2019 surgyr_2020 surgyr_2021 surgyr_2022 surgyr_2023

* Cluster variable
global clustervar HospitalCluster

* Treatment and instrument
global Treated early_surgery
global Z       instrumental_variable

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
* MAIN LOOP — timepoints
* ============================================================================

foreach horizon of local horizons {

    global Horizon = `horizon'

    * ── Load data ────────────────────────────────────────────────────
    local dta_path     = "$data_dir" + "elective_`horizon'd.dta"
    local globals_path = "$data_dir" + "elective_`horizon'd_globals.csv"

    use "`dta_path'", clear
    encode nvrhospitalname, gen(HospitalCluster)
    gen All = 1
    gen gender_M         = 1 - gender_F
    gen fontaine_not4    = 1 - fontaine_4
    gen comorbidity_not1 = 1 - comorbidity_1
    gen krt_no           = 1 - krt_yn
    gen surgyr_2016 = (surgyr_2017 == 0 & surgyr_2018 == 0 & surgyr_2019 == 0 & surgyr_2020 == 0 & surgyr_2021 == 0 & surgyr_2022 == 0 & surgyr_2023 == 0)

    * ── Read globals CSV ─────────────────────────────────────────────
    import delimited using "`globals_path'", clear varnames(1)

    local n_outcomes = _N
    forvalues i = 1/`n_outcomes' {
        local outcome_`i' = outcome[`i']
        local family_`i'  = family[`i']
    }

    * Reload analysis data
    use "`dta_path'", clear
    encode nvrhospitalname, gen(HospitalCluster)
    gen All = 1
    gen gender_M         = 1 - gender_F
    gen fontaine_not4    = 1 - fontaine_4
    gen comorbidity_not1 = 1 - comorbidity_1
    gen krt_no           = 1 - krt_yn
    gen surgyr_2016 = (surgyr_2017 == 0 & surgyr_2018 == 0 & surgyr_2019 == 0 & surgyr_2020 == 0 & surgyr_2021 == 0 & surgyr_2022 == 0 & surgyr_2023 == 0)

    * Build Outcomelist global
    global Outcomelist
    forvalues i = 1/`n_outcomes' {
        global Outcomelist $Outcomelist `outcome_`i''
    }

    * ── Loop over outcomes ───────────────────────────────────────────
    forvalues i = 1/`n_outcomes' {

        local outcome = "`outcome_`i''"
        local family  = "`family_`i''"

        * Set Xlist globals by re-reading globals CSV for this outcome
        * (direct string assignment avoids Stata local macro length limits)
        preserve
        import delimited using "`globals_path'", clear varnames(1)
        global Xlist_s2_m1   = xlist_s2_m1[`i']
        global Xlist_s1_m1   = xlist_s1_m1[`i']
        global Xlist_s2_m2   = xlist_s2_m2[`i']
        global Xlist_s1_m2   = xlist_s1_m2[`i']
        global Zlist_s1_m2   = zlist_s1_m2[`i']
        global Dxlist_s2_m2  = dxlist_s2_m2[`i']
        global Xlist_s2_noiv = xlist_s2_noiv[`i']
        restore
        global OutcomeFamily `family'
		
        * ── IV strength ──────────────────────────────────────────────
        compute_iv_strength, version("model1") outcome("`outcome'")
		compute_iv_strength, version("model2") outcome("`outcome'")
        * no_iv has no first stage — IV strength not applicable

		* Propensity score overlap — all horizons
        * Both model1 and model2 have outcome-specific Stage 1 specs
        *plot_ps_overlap, version("model1") outcome("`outcome'") /// * TO DO: these plots are currently overwriting each other — need to add version and horizon to filename
            horizon(`horizon') results_dir("$results_dir")
        plot_ps_overlap, version("model2") outcome("`outcome'") ///
            horizon(`horizon') results_dir("$results_dir")
            
        * ── Bootstrap loop over versions ─────────────────────────────
        foreach version of global versions {

            global OutcometoUse `outcome'
            SetOutcomeNumber

            * Select pred types by family
            if "$OutcomeFamily" == "gaussian" local predlist xb ystar pr e
			if "$OutcomeFamily" == "binomial" local predlist pr

            * Build bootstats local
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

            use "`dta_path'", clear
            encode nvrhospitalname, gen(HospitalCluster)
            gen All = 1
            gen gender_M         = 1 - gender_F
            gen fontaine_not4    = 1 - fontaine_4
            gen comorbidity_not1 = 1 - comorbidity_1
            gen krt_no           = 1 - krt_yn
            gen surgyr_2016 = (surgyr_2017 == 0 & surgyr_2018 == 0 & surgyr_2019 == 0 & surgyr_2020 == 0 & surgyr_2021 == 0 & surgyr_2022 == 0 & surgyr_2023 == 0)
        }
    }
}

* ── Save IV strength results ──────────────────────────────────────────
postclose iv_strength
use `iv_strength_file', clear
* Cast to string to ensure proper export
tostring outcome version, replace force
export delimited using "${results_dir}iv_strength_all.csv", replace