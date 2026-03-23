# Charlson_SCARF_scoring.R
#
# ICD-10 code dictionaries and deficit lists for:
#   - Charlson comorbidity index (RCS version)
#   - SCARF frailty index
#
# Source this file before calling add_flags_from_concat(), add_rcs_score(),
# add_scarf_score(), or add_scarf_admin_flags().

# =============================================================================
# 1. CHARLSON
# =============================================================================

icd10_groups_charlson <- list(
  mi     = c("I252"),
  premi  = c("I21", "I22", "I23"),
  ccf    = c("I11", "I13", "I42", "I43", "I50", "I255", "I517"),
  pvd    = c("I70", "I71", "I72", "I73", "I770", "I771", "K551", "K558", "K559",
             "R02", "Z958", "Z959"),
  cvd    = c("G45", "G46", "I6"),
  dem    = c("A810", "F00", "F01", "F02", "F03", "F051", "G30", "G31"),
  cpd    = c("I26", "I27", "J40", "J41", "J42", "J43", "J44", "J45", "J47",
             "J60", "J61", "J62", "J63", "J64", "J65", "J66", "J67", "J684",
             "J701", "J703"),
  precpd = c("J46"),
  rhd    = c("M05", "M06", "M09", "M120", "M315", "M32", "M33", "M34", "M35", "M36"),
  liv    = c("B18", "I85", "I864", "I982", "K70", "K71", "K721", "K729", "K76",
             "R162", "Z944"),
  dbm    = c("E10", "E11", "E12", "E13", "E14"),
  ple    = c("G114", "G81", "G82", "G83"),
  ren    = c("I12", "I13", "N01", "N03", "N05", "N07", "N08", "N18", "N25",
             "Z49", "Z940", "Z992"),
  preren = c("N171", "N172", "N19"),
  mal    = c("C0", "C1", "C2", "C3", "C40", "C41", "C43", "C45", "C46", "C47",
             "C48", "C49", "C5", "C6", "C70", "C71", "C72", "C73", "C74", "C75",
             "C76", "C80", "C81", "C82", "C83", "C84", "C85", "C88", "C9"),
  mst    = c("C77", "C78", "C79"),
  aid    = c("B20", "B21", "B22", "B23", "B24")   # defined but excluded from scoring
)

# Comorbidities scored at index admission
charlson_all_list  <- c("mi", "ccf", "pvd", "cvd", "dem", "cpd", "rhd", "liv",
                        "dbm", "ple", "ren", "mal", "mst")

# Comorbidities that should only be counted in pre-index admissions
# (defined for completeness; not used in index-admission-only analyses)
charlson_pre_list  <- c("premi", "precpd", "preren")

# =============================================================================
# 2. SCARF
# =============================================================================

icd10_groups_frailty <- list(
  dis_scarf    = c("R26", "S78", "S88", "S98", "T136", "Y83", "Z993", "G11",
                   "G81", "G82", "G83", "M62"),
  pd_scarf     = c("G122", "G20", "G21", "G22", "G23", "G25", "G26", "G32",
                   "G35", "R25"),
  care_scarf   = c("R40", "Z50", "Z74", "Z755", "Z998", "Z999"),
  soc_scarf    = c("F1", "R460", "R468", "Y06", "Z59", "Z60", "Z63", "Z73"),
  hear_scarf   = c("H833", "H90", "H91", "Z453", "Z461", "Z974"),
  vis_scarf    = c("H25", "H28", "H35", "H40", "H43", "H53", "H54"),
  fall_scarf   = c("S00", "S01", "R296", "W00", "W01", "W04", "W05", "W06",
                   "W07", "W08", "W10", "W18", "W19"),
  ulc_scarf    = c("I83", "I89", "L03", "L08", "L89", "L97", "L984"),
  inc_scarf    = c("N31", "N393", "N394", "R15", "R32", "T835", "Z466"),
  nut_scarf    = c("E41", "E43", "E44", "E46", "E53", "E55", "E66", "E83",
                   "E87", "R53", "R628", "R63", "R64", "X53"),
  cog_scarf    = c("F00", "F01", "F02", "F03", "F04", "F05", "F067", "G30",
                   "G31", "R41", "R54"),
  mood_scarf   = c("F2", "F3", "F41", "R44", "R45"),
  anae_scarf   = c("D50", "D51", "D52", "D53", "D63", "D64"),
  arth_scarf   = c("M05", "M06", "M07", "M09", "M10", "M11", "M12", "M13",
                   "M15", "M16", "M17", "M18", "M19", "M315", "M32", "M33",
                   "M34", "M35", "M36"),
  af_scarf     = c("I44", "I48", "I49", "Z450", "Z950"),
  cvd_scarf    = c("G45", "G46", "I6"),
  ren_scarf    = c("I12", "I13", "N01", "N03", "N05", "N07", "N08", "N18",
                   "N19", "N25", "I770", "Z49", "Z940", "Z992"),
  dm_scarf     = c("E109", "E119", "E129", "E139", "E149"),
  dmcomp_scarf = c("E100", "E101", "E102", "E103", "E104", "E105", "E106",
                   "E107", "E108", "E111", "E112", "E113", "E114", "E115",
                   "E116", "E117", "E118", "E120", "E121", "E122", "E123",
                   "E124", "E125", "E126", "E127", "E128", "E130", "E131",
                   "E132", "E133", "E134", "E135", "E136", "E137", "E138",
                   "E140", "E141", "E142", "E143", "E144", "E145", "E146",
                   "E147", "E148", "G590", "G632", "H360", "M142", "M146",
                   "N083"),
  hf_scarf     = c("I11", "I13", "I260", "I27", "I42", "I43", "I50", "I51"),
  valv_scarf   = c("I05", "I06", "I07", "I08", "I34", "I35", "I36", "I37",
                   "I390", "I391", "I392", "I393", "I394", "Z952", "Z953",
                   "Z954"),
  htn_scarf    = c("I10", "I11", "I12", "I13", "H350"),
  hypo_scarf   = c("I95", "R42", "R55", "E86", "E87"),
  ihd_scarf    = c("I20", "I21", "I22", "I23", "I24", "I250", "I251", "I252",
                   "I253", "I254", "I255", "I256", "I258", "I259"),
  foot_scarf   = c("B353", "G575", "G576", "L60", "M201", "M202", "M203",
                   "M204", "M205", "M206", "M213", "M214", "M215", "M216",
                   "M722", "M766", "M773", "M775", "S90", "S91", "S92", "S93",
                   "S94", "S97", "S99", "Q66"),
  frac_scarf   = c("M80", "M484", "S22", "S32", "S33", "S42", "S43", "S62",
                   "S72", "S73"),
  op_scarf     = c("M80", "M81", "M82"),
  pud_scarf    = c("K21", "K25", "K26", "K27", "K29", "R12"),
  pvd_scarf    = c("I65", "I70", "I71", "I72", "I73", "I771", "K551", "K558",
                   "K559", "R02", "Z958", "Z959"),
  resp_scarf   = c("J46", "J45", "J40", "J41", "J42", "J43", "J44", "J47",
                   "J60", "J61", "J62", "J63", "J64", "J65", "J684", "J70",
                   "J13", "J14", "J15", "J16", "J18", "J20", "J22", "J90",
                   "J961", "J980", "R06"),
  thy_scarf    = c("E03", "E04", "E05", "E06", "E079"),
  uri_scarf    = c("N30", "N34", "N390", "N398", "N399", "R31", "R33", "T835")
)

scarf_deficits <- paste0(
  c("cog", "mood", "soc", "dis", "fall", "nut", "anae", "af", "ihd", "valv",
    "hf", "htn", "hypo", "cvd", "pvd", "dm", "dmcomp", "ren", "resp", "arth",
    "op", "foot", "pd", "pud", "ulc", "thy", "vis", "hear", "inc", "frac",
    "care", "uri"),
  "_scarf"
)

# HES admin codes indicating care home admission/discharge (for SCARF care deficit)
scarf_carehome_admisorc <- c(54, 65, 85)
scarf_carehome_disdest  <- c(54, 65, 85)
