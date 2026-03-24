* non_elective_run_clinical_effectiveness.do
*
* Runs 2SRI IV clinical effectiveness analysis for non-elective bypass cohort.
* Reads LASSO-selected variable lists from globals CSV (produced by R).
* Runs four Stage 1 versions, bootstrapped, for all outcomes x timepoints.
* Produces forest plots and IV strength diagnostics.
*
* D*X and Z*X interaction columns are pre-generated in R and included
* directly in the appropriate Xlist — no factor variable syntax in Stata.

* ============================================================================
* CONFIGURATION — edit here only
* ============================================================================

* Timepoints to analyse (must match DTA/globals files produced by R)
local horizons 90 180 365 

* Bootstrap settings
global nreps  = 300       // number of bootstrap replications
global seed   = 37563845  // bootstrap seed
global withIV = 1         // 1 = include generalised residual; 0 = naive

* Stage 1 versions to run
global versions z_only z_x_stage2_treatment z_x_stage1_full z_x_stage1_instrument // no_iv

* Paths
global data_dir     "Z:/PHP/HSR/ESORT-V/ESORT-V/Akshay_Scripts_Bypass_TTE_180226/analysable_subsets/march23_lasso_outputs/"
global results_dir  "Z:/PHP/HSR/ESORT-V/ESORT-V/Akshay_Scripts_Bypass_TTE_180226/clinical_effectiveness_results/"
global programs_dir "C:/Users/LSHAJ82/Documents/GitHub/R-tte/stata/"

* Subgroups for recycled predictions
global Subgrouplist All

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
    local dta_path     = "$data_dir" + "non_elective_`horizon'd.dta"
    local globals_path = "$data_dir" + "non_elective_`horizon'd_globals.csv"

    use "`dta_path'", clear
    encode nvrhospitalname, gen(HospitalCluster)
    gen All = 1

    * ── Read globals CSV ─────────────────────────────────────────────
    * Columns: outcome, family, xlist_s2, xlist_s1_v1, xlist_s1_v2,
    *          xlist_s1_v3, xlist_s1_v4
    * See mapping table in CLAUDE.md for full contents of each xlist.
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
        global Xlist_s2    = xlist_s2[`i']
        global Xlist_s1_v1 = xlist_s1_v1[`i']
        global Xlist_s1_v2 = xlist_s1_v2[`i']
        global Xlist_s1_v3 = xlist_s1_v3[`i']
        global Xlist_s1_v4 = xlist_s1_v4[`i']
        restore
        global OutcomeFamily `family'

        * ── IV strength ──────────────────────────────────────────────
        * Versions 1-3: first stage Xlist same across outcomes — run once
        * per version using first outcome only (F-stat unaffected by outcome)
        * Version 4: Xlist differs per outcome — run for each outcome
        if `i' == 1 {
            foreach version in z_only z_x_stage2_treatment z_x_stage1_full {
                compute_iv_strength, version("`version'") outcome("`outcome'")
            }
        }
        compute_iv_strength, version("z_x_stage1_instrument") outcome("`outcome'")

        * ── Bootstrap loop over versions ─────────────────────────────
        foreach version of global versions {

            global OutcometoUse `outcome'
            SetOutcomeNumber

            * Select pred types by family
            if "`family'" == "gaussian" local predlist xb ystar pr e
            if "`family'" == "binomial" local predlist pr

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
                results_dir("$results_dir")
        }
    }

    * ── Forest plot for this timepoint ───────────────────────────────
    forest_plot, horizon(`horizon') results_dir("$results_dir") data_dir("$data_dir")
}

* ── Save IV strength results ──────────────────────────────────────────
postclose iv_strength
use `iv_strength_file', clear
export delimited using "${results_dir}iv_strength_all.csv", replace
