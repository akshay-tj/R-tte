* non_elective_run_clinical_effectiveness.do
*
* Runs 2SRI IV clinical effectiveness analysis for non-elective bypass cohort.
* Reads LASSO-selected variable lists from globals CSV (produced by R).
* Runs four Stage 1 versions, bootstrapped, for all outcomes x timepoints.
* Produces forest plots and IV strength diagnostics.

* ============================================================================
* CONFIGURATION — edit here only
* ============================================================================

* Timepoints to analyse (must match DTA/globals files produced by R)
local horizons 90 // changed from 90 180 365

* Bootstrap settings
global nreps      = 20        // number of bootstrap replications
global seed       = 37563845   // bootstrap seed
global withIV     = 1          // 1 = include generalised residual; 0 = naive

* Stage 1 versions to run
global versions z_only z_x_stage2_treatment z_x_stage1_full z_x_stage1_instrument

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

cd "$results_dir"
do "${programs_dir}iv_2sri_programs.do"

* Initialise IV strength postfile
tempfile iv_strength_file
capture postclose iv_strength
postfile iv_strength str50 outcome str50 version ///
    double eff_f double crit5 ///
    using `iv_strength_file', replace

* ============================================================================
* MAIN LOOP — timepoints
* ============================================================================

foreach horizon of local horizons {

    global Horizon = `horizon'

    * ── Load data ────────────────────────────────────────────────────
    use "${data_dir}non_elective_`horizon'd.dta", clear
    encode nvrhospitalname, gen(HospitalCluster)
    gen All = 1

    * ── Read globals CSV — sets lists per outcome ────────────────────
    * globals CSV columns:
    *   outcome, family,
    *   stage2_main_effects_and_xx, stage1_main_effects_and_xx,
    *   stage2_treatment_interactions, stage1_instrument_interactions

    import delimited using ///
        "${data_dir}non_elective_`horizon'd_globals.csv", ///
        clear varnames(1)
    * Store in matrices/locals for lookup below
    local n_outcomes = _N
    forvalues i = 1/`n_outcomes' {
        local outcome_`i'  = outcome[`i']
        local family_`i'   = family[`i']
        local xlist2_`i'   = stage2_main_effects_and_xx[`i']
        local xlist1_`i'   = stage1_main_effects_and_xx[`i']
        local mlist_`i'    = stage2_treatment_interactions[`i']
        local zlist_`i'    = stage1_instrument_interactions[`i']
    }

    * Reload analysis data
    use "${data_dir}non_elective_`horizon'd.dta", clear
    encode nvrhospitalname, gen(HospitalCluster)
    gen All = 1

    * Build Outcomelist global from globals CSV outcomes
    global Outcomelist
    forvalues i = 1/`n_outcomes' {
        global Outcomelist $Outcomelist `outcome_`i''
    }

    * ── Loop over outcomes ───────────────────────────────────────────
    forvalues i = 1/`n_outcomes' {

        local outcome  = "`outcome_`i''"
        local family   = "`family_`i''"

        * Set globals for this outcome
        global Xlist_stage2  `xlist2_`i''
        global Xlist_stage1  `xlist1_`i''
        global Mlist         `mlist_`i''
        global Zlist         `zlist_`i''
        global OutcomeFamily `family'
		
		* Derive binary-only Mlist and Zlist — ageatsurgery handled
        * separately via c.$Z#c.ageatsurgery and c.$Treated#c.ageatsurgery
		local cont_vars ageatsurgery age_sq ageatsurgery_x_age_sq
        global Mlist_bin       : list global(Mlist)       - local(cont_vars)
        global Zlist_bin       : list global(Zlist)       - local(cont_vars)
        global Xlist_stage2_bin : list global(Xlist_stage2) - local(cont_vars)

        * ── IV strength ──────────────────────────────────────────────
        * NOTE: weakivtest requires linear IV — skipping for now, revisit
		*foreach version of global versions {
        *    compute_iv_strength, ytouse("`outcome'") version("`version'")
        *}

        * ── Bootstrap loop over versions ────────────────────────────
        foreach version of global versions {

            * Build bootstats local for this outcome x version
            global OutcometoUse `outcome'
            SetOutcomeNumber

            * Select pred type
            if "`family'" == "gaussian" local predlist xb ystar pr e
            if "`family'" == "binomial" local predlist pr

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

            * Save bootstrap results to CSV
            preserve
            use "bootstrap_`outcome'_`version'_`horizon'd.dta", clear
            export delimited using ///
                "bootstrap_`outcome'_`version'_`horizon'd.csv", replace
            restore
        }
    }

    * ── Forest plot for this timepoint ──────────────────────────────
    forest_plot, horizon(`horizon')
}

* ── Save IV strength results ─────────────────────────────────────────
postclose iv_strength
use `iv_strength_file', clear
export delimited using "${results_dir}iv_strength_all.csv", replace