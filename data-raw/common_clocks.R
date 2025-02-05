## code to prepare `common_clocks` dataset goes here
cpg_extended <- c(
  "DNAmGait_noAge",
  "DNAmGrip_noAge",
  "DNAmVO2max",
  "DNAmGait_wAge",
  "DNAmGrip_wAge",
  "DNAmFEV1_noAge",
  "tnsc2",
  "irS2",
  "DNAmFitAge",
  "DNAmGrimAge2BasedOnPredictedAge",
  "DNAmGrimAgeBasedOnPredictedAge"
)

cpg_core <- c(
  "DNAmTL",
  "DNAmAge",
  "DNAmAgeHannum",
  "DNAmPhenoAge",
  "DNAmAgeSkinBloodClock",
  "DNAmGrimAge2BasedOnRealAge",
  "DNAmGrimAgeBasedOnRealAge",
  "DunedinPACE"
)

SystemsAge <- c(
  "Blood",
  "Brain",
  "Inflammation",
  "Heart",
  "Hormone",
  "Immune",
  "Kidney",
  "Liver",
  "Metabolic",
  "Lung",
  "MusculoSkeletal",
  "Age_prediction",
  "SystemsAge"
)

pc_clock <- c(
  "PCHorvath1",
  "PCHorvath2",
  "PCHannum",
  "PCPhenoAge",
  "PCDNAmTL",
  "PCGrimAge"
)

cpg_grim <- c(
  "DNAmADM",
  "DNAmB2M",
  "DNAmCystatinC",
  "DNAmGDF15",
  "DNAmLeptin",
  "DNAmPACKYRS",
  "DNAmPAI1",
  "DNAmTIMP1",
  "DNAmGDF_15",
  "DNAmCystatin_C",
  "DNAmTIMP_1",
  "DNAmadm",
  "DNAmpai_1",
  "DNAmleptin"
)

pc_grim <- c(
  "PCPACKYRS",
  "PCADM",
  "PCB2M",
  "PCCystatinC",
  "PCGDF15",
  "PCLeptin",
  "PCPAI1",
  "PCTIMP1"
)

common_clocks <- list(
  cpg_core = cpg_core,
  cpg_extended = cpg_extended,
  SystemsAge = SystemsAge,
  pc_clock = pc_clock,
  cpg_grim = cpg_grim,
  pc_grim = pc_grim
)

usethis::use_data(common_clocks, overwrite = TRUE)
