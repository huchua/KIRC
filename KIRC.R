install.packages(c('ggplot2', 'plyr', 'glmnet', 'survival','dplyr', 'plot3D', 'readr'))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c('ComplexHeatmap', 'survcomp'))

library(plyr)
library(ggplot2)
library(glmnet)
library(survival)
library(dplyr)
library(ComplexHeatmap)
library(plot3D)
library(gmodels)

cancer <- "KIRC"
num <- 0.6
color_type = c('red', 'blue', 'green', 'yellow')

###function
GetExpSurCli <- function(rt_T, rt_cli, rt_sur){
  print('the colnames of expression table should contain "Ensembl_ID", 
        you will get a list containing expression table, clinical data and survival data
        all of which contain the same sample id and the same order, 
        the expression data contain only tumor samples')
  
  ###survival data
  sam_sur <- matrix(unlist(strsplit(rt_sur$sample, "-")), ncol = 4, byrow = TRUE)
  rt_sur$sample <- paste(sam_sur[, 1], sam_sur[, 2], sam_sur[, 3], sep = "-")
  
  #clinical data
  print("the first column must be samples id which was combined wiht '-'")
  sam_cli <- matrix(unlist(strsplit(rt_cli[, 1], "-")), ncol = 4, byrow = TRUE)
  rt_cli$submitter_id.samples <- paste(sam_cli[, 1], sam_cli[, 2], sam_cli[, 3], sep = "-")
  
  ###intersection sampples
  sam_m <- intersect(colnames(rt_T),intersect(rt_sur$sample, rt_cli$submitter_id.samples))
  
  pair_T <- match(sam_m, colnames(rt_T), nomatch = NA)
  pair_S <- match(sam_m, rt_sur$sample , nomatch = NA)
  pair_C <- match(sam_m, rt_cli$submitter_id.samples , nomatch = NA)
  
  rt_T_m <- rt_T[, pair_T]
  rt_S_m <- rt_sur[pair_S, ]
  rt_C_m <- rt_cli[pair_C, ]
  
  esc_list <- list(rt_T_m,  rt_C_m, rt_S_m)
  if(all(colnames(rt_T_m) == rt_S_m$sample)&all(rt_S_m$sample == rt_C_m$submitter_id.samples)){
    return(esc_list)
  } else {
    stop("the sample id is diff among expression data, survival data and clinical data ")
  }
}


GEOExpCliM <- function(rt_exp, rt_cli_sur){
  int_id <- intersect(colnames(rt_exp), rt_cli_sur$Accession)
  rt_exp_m <- rt_exp[, match(int_id, colnames(rt_exp), nomatch = 0)]
  rt_cli_sur_m <- rt_cli_sur[match(int_id, rt_cli_sur$Accession, nomatch = 0), ]
  esc_list <- list(rt_exp_m,  rt_cli_sur_m)
  return(esc_list)
}


TCGAExpSurM <- function(rt_exp, rt_sur){
  
  
  ###survival data
  sam_sur <- matrix(unlist(strsplit(rt_sur$sample, "-")), ncol = 4, byrow = TRUE)
  rt_sur$sample <- paste(sam_sur[, 1], sam_sur[, 2], sam_sur[, 3], sep = "-")
  
  ###intersection sampples
  int_id <- intersect(colnames(rt_exp), rt_sur$sample)
  
  rt_exp_m <- rt_exp[, match(int_id, colnames(rt_exp), nomatch = 0)]
  rt_sur_m <- rt_sur[match(int_id, rt_sur$sample , nomatch = 0), ]
  
  esc_list <- list(rt_exp_m,  rt_sur_m)
  if(all(colnames(rt_exp_m) == rt_sur_m$sample)){
    return(esc_list)
  } else {
    stop("the sample id is diff among expression data, survival data and clinical data ")
  }
}

MergeCliOth <- function(rt_cli, rt_other){
  sam_int <- intersect(rt_cli$sample_id, row.names(rt_other))
  rt_cli_m <- rt_cli[match(sam_int, rt_cli$sample_id, nomatch = 0), ]
  rt_other_m <- rt_other[match(sam_int, row.names(rt_other), nomatch = 0), ]
  rt_other_cli <- data.frame(cbind(rt_other_m, rt_cli_m), stringsAsFactors = FALSE)
  return(rt_other_cli)
}

split_tcga_tn <-function(rt_tcga, sam_type = c("tumor", "normal"), split_type = c("[.]", "-", "_")){
  #the colnames is samples id and row.name is gene id 
  #the format of sample id must be "TCGA.76.4927.01A" 
  #sam_type must be set
  get_sam_nam <- function(rt_tcga, split_type){
    if(split_type == "[.]"){
      sam_G <- matrix(unlist(strsplit(colnames(rt_tcga), "[.]")), ncol = 4, byrow = TRUE)
    }else if(split_type == "-"){
      sam_G <- matrix(unlist(strsplit(colnames(rt_tcga), "-")), ncol = 4, byrow = TRUE)
    } else if(split_type == "_") {
      sam_G <- matrix(unlist(strsplit(colnames(rt_tcga), "_")), ncol = 4, byrow = TRUE)
    } else {
      stop("you must set split_type ")
    }
    sam_nam <- paste(sam_G[, 1], sam_G[, 2], sam_G[, 3], sep = "-")
    sam_list <- list(sam_G, sam_nam)
    return(sam_list)
  }
  
  sam_list <- get_sam_nam(rt_tcga, "[.]")
  sam_G <- sam_list[[1]]
  sam_nam <- sam_list[[2]]
  
  if(sam_type == "tumor"){
    tcga_pos_T <- c(grep("01", sam_G[, 4]), grep("02", sam_G[, 4]), grep("03", sam_G[, 4]), grep("04", sam_G[, 4]), grep("05", sam_G[, 4]), grep("06", sam_G[, 4]), grep("07", sam_G[, 4]), grep("08", sam_G[, 4]), grep("09", sam_G[, 4]))
    rt_T <- rt_tcga[, tcga_pos_T]
    colnames(rt_T) <- sam_nam[tcga_pos_T]
    data <- rt_T
  } else if(sam_type == "normal"){
    tcga_pos_N <- c(grep("11", sam_G[, 4]))
    if(length(tcga_pos_N) > 0) {
      rt_N <- rt_tcga[, tcga_pos_N]
      colnames(rt_N) <- sam_nam[tcga_pos_N]
      data <- rt_N
    } else {
      warning("this is no normal samples in this cancer")
      next
    }
  } else {
    stop("you must set sam_type to point which type samples you want to get")
  }
}

CliSort <- function(rt_cli, cancer){
  if(cancer == "LUAD"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic", "number_pack_years_smoked", "race.demographic",
                         "new_neoplasm_event_type", "pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "smoke", "race", "recurrence", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))#female male
    rt_esc$race <- factor(rt_esc$race, labels=c(1:5))#american indian or alaska native asian black or african american not reported white
    rt_esc$recurrence <- factor(rt_cli_m$new_tumor_event_after_initial_treatment, labels=c(3, 2, 1))#''  NO YES
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 4, 5))#T1 T1a T1b T2 T2a T2b T3 T4 TX
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 2, 3))#''  N0 N1 N2 N3 NX
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 2, 3))#''  M0 M1 M1a M1b MX
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 1, 1, 1, 1, 2, 2, 2))#not reported stage i stage ia stage ib stage ii stage iia stage iib stage iiia stage iiib stage iv
  }
  if(cancer == "UVM"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses",  "new_tumor_event_after_initial_treatment")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage", "recurrence")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3))
    rt_esc$N <- factor(rt_esc$N, labels=c(2, 1, 2))
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 3))
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 2, 2, 2, 2))
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 1, 2))
  }
  if(cancer == "BRCA"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5))
    rt_esc$N <- factor(rt_esc$N, labels=c(rep(1, 4), rep(2, 11), 3))
    rt_esc$M <- factor(rt_esc$M, labels=c(1, 2, 3, 1))
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3))
  }
  if(cancer == "ESCA"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id",  "age", "gender", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(6, 1, 2, 3, 4, 5, 5))
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 2, 3))
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 3))
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2))
  }
  if(cancer == "HNSC"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(6, 1, 2, 3, 4, 5, 5, 5, 6))
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 2, 2, 2, 2, 3))
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 3))
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 2, 2, 2, 2))
  }
  if(cancer == "KIRP"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic", "tumor_type", "pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses", "new_tumor_event_after_initial_treatment", "additional_therapy")]
    colnames(rt_esc) <- c("sample_id", "age", "gender","tumor_type", "T", "N", "M", "stage", "recurrence", "additional_therapy")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))# female male
    rt_esc$tumor_type <- factor(rt_esc$tumor_type, labels=c(3, 1, 2))# "" Type 1 Type 2
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5))# T1 T1a T1b T2 T2a T2b T3 T3a T3b T3c T4 TX
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 3))#"" N0 N1 N2 NX
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 3))#"" M0 M1 MX
    rt_esc$stage <- factor(rt_esc$stage, labels = c(3, 1, 1, 2, 2))# not reported stage i stage ii stage iii stage iv
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels = c(3, 2, 1))#  NO YES
    rt_esc$additional_therapy <- factor(rt_esc$additional_therapy, labels=c(2, 1)) #"NO" "YES"
  }
  if(cancer == "KIRC"){
    rt_esc <-  rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses", "new_tumor_event_after_initial_treatment")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage", "recurrence")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4))
    rt_esc$N <- factor(rt_esc$N, labels=c(1, 2, 3))
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 3))
    rt_esc$stage <- factor(rt_esc$stage, labels = c(3, 1, 1, 2, 2))
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels = c(3, 2, 1))
  }
  if(cancer == "LAML"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_diagnosis.diagnoses", "gender.demographic", "lab_procedure_hemoglobin_result_specified_value",
                         "lab_procedure_bone_marrow_blast_cell_outcome_percent_value", "lab_procedure_leukocyte_result_unspecified_value", "platelet_result_count", 
                         "lab_procedure_blast_cell_outcome_percentage_value", "leukemia_french_american_british_morphology_code", 
                         "acute_myeloid_leukemia_calgb_cytogenetics_risk_category", "cytogenetic_abnormality")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "hemoglobin", "blast_cell_BM", "leukocyte", "platelet", "blast_cell_PB", 
                          "AML_FAB_subtype", "cytogenetic_risk_group", "cytogenetic_abnormality")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels = c(1, 2))#female male
    rt_esc$AML_FAB_subtype <- factor(rt_esc$AML_FAB_subtype, labels = c(1, 2, 3, 4, 5, 6, 7, 8)) #"M0 Undifferentiated" M1 M2 M3 M4 M5 M6 M7 Not Classified
    rt_esc$cytogenetic_risk_group <- factor(rt_esc$cytogenetic_risk_group, labels = c(3, 1, 2, 2)) #"" "Favorable" "Intermediate/Normal" "Poor"
    rt_esc$cytogenetic_abnormality <- factor(rt_esc$cytogenetic_abnormality, labels = c(3, 2, 2, 1, 2, 2, 2, 2, 2, 2)) #" " "+8" "Complex - >= 3 distinct abnormalities" "Normal" "del(5q) / 5q-" "del(7q) / 7q-" "inv(16)" "t(15;17) and  variants" "t(8;21)" "t(9;11)"
    return(rt_esc)
  } 
  if(cancer == "LGG"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic", "neoplasm_histologic_grade")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "grade")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$grade <- factor(rt_esc$grade, labels=c(3, 1, 2))
  } 
  if(cancer == "OV"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "clinical_stage")]
    colnames(rt_esc) <- c("sample_id", "age", "clinical_stage")
    
    ###convert clinical element to numeric
    rt_esc$clinical_stage <- factor(rt_esc$clinical_stage, labels=c(3, 1, 1, 1, 1, 2, 2, 2, 2))
  } 
  if(cancer == "SKCM"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(5, 1, 2, 2, 2, 3, 3, 4, 4, 4))
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, rep(2, 8), 3))
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 2))
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 3, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2))
  } 
  if(cancer == "LUSC"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", 
                         "pathologic_N", "pathologic_M", "tumor_stage.diagnoses", 'new_tumor_event_after_initial_treatment')]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage", 'recurrence')
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2)) #female male
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 4)) # T1 T1a T1b T2 T2a T2b T3 T4
    rt_esc$N <- factor(rt_esc$N, labels=c(1, rep(2, 3), 3)) #N0 N1 N2 N3 NX
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 2, 3)) #  M0 M1 M1a M1b MX
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2)) # not reported stage i stage ia stage ib stage ii stage iia stage iib stage iii stage iiia stage iiib stage iv
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 2, 1)) #  NO YES
  } 
  if(cancer == "LIHC"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", 
                         "pathologic_N", "pathologic_M", "tumor_stage.diagnoses", 'new_tumor_event_after_initial_treatment')]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage", 'recurrence')
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2)) #female male
    rt_esc$T <- factor(rt_esc$T, labels=c(5, 1, 2, 2, 2, 3, 3, 3, 4, 5)) #'' T1 T2 T2a T2b T3 T3a T3b T4 TX
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 3)) #'' N0 N1 NX
    rt_esc$M <- factor(rt_esc$M, labels=c(1, 2, 3)) # M0 M1 MX
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 2, 2, 2, 2, 2, 2, 2)) #  not reported stage i stage ii stage iii stage iiia stage iiib stage iiic stage iv stage iva stage ivb
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 2, 1)) #'' NO YES
  }
  if(cancer == "READ"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", 
                         "pathologic_N", "pathologic_M", "tumor_stage.diagnoses", 'new_tumor_event_after_initial_treatment')]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage", 'recurrence')
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2)) #female male
    rt_esc$T <- factor(rt_esc$T, labels=c(5, 1, 2, 3, 4, 4, 4)) #'' T1 T2 T3 T4 T4a T4b
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, rep(2, 7), 3)) #'' N0 N1 N1a N1b N1c N2 N2a N2b NX
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 3)) #'' M0 M1 M1a MX
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, rep(1, 5), rep(2, 6))) #not reported stage i stage ii stage iia stage iib stage iic stage iii stage iiia stage iiib stage iiic stage iv stage iva
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 2, 1)) #'' NO YES
  }
  if(cancer == "SARC"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic", 'new_tumor_event_after_initial_treatment', 'tumor_depth')]
    colnames(rt_esc) <- c("sample_id", "age", "gender", 'recurrence', 'tumor_depth')
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2)) #female male
    rt_esc$tumor_depth <- factor(rt_cli_m$tumor_depth, labels=c(3, 2, 1))#'' Deep Superficial
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 2, 1)) #'' NO YES
  }
  rt_esc_c <- data.frame(apply(rt_esc, 2, as.character), stringsAsFactors = FALSE)
  colnames(rt_esc_c) <- colnames(rt_esc)
  write.table(rt_esc_c, file = paste(cancer, "_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  return(rt_esc_c)
}

GLCoxMain <- function(k, x, y, alpha = 1, s = "lambda.min", family = family, width, height){
  ScreenKeyGS <- function(k, family = family){
    #x is the expression table, the row.names is samples id and colname is the gene id
    #y id the survival data which contain time and status
    #k is the frequency times that you want to train in the model
    #width is the boxplot width and height is the boxplot height
    #theme_coef is theme of boxplot
    #x_coef is expression matrix contained active covariates in the model
    #risk_genes with the active covariates
    #coef_value_d is the coeffecient value of the risk genes
    if(family == "cox"){
      fit = glmnet(x, y, family = "cox", alpha = alpha)
      cvfit = cv.glmnet(x, y, family = "cox", alpha = alpha)
    } else if (family == "binomial"){
      fit = glmnet(x, y, family = "binomial", alpha = alpha)
      cvfit = cv.glmnet(x, y, family = "binomial", alpha = alpha)
    }
    
    pdf(file = 'coef_L1.pdf', 10, 6)
    plot(cvfit$glmnet.fit, label = TRUE, xvar = "lambda")
    dev.off()
    
    pdf(file = 'dev_lab.pdf', 10, 6)
    plot(cvfit)
    dev.off()
    coef_min = coef(cvfit, s = s)
    lambda = as.numeric(cvfit[s])
    x_coef <- x[, which(coef_min[which(row.names(coef_min) != "(Intercept)"), 1] != 0)]
    risk_genes <- colnames(x_coef)
    coef_value_d <- coef_min[which(coef_min[, 1] != 0)]
    risk_genes_lam <- paste(risk_genes, collapse = ";")
    coef_list <- list(risk_genes_lam, x_coef, coef_value_d, lambda)
    return(coef_list)
  }
  
  CountGS <- function(t, risk_gens_vec, UG_list){
    pos_genes <- which(risk_gens_vec == UG_list[[t]])[1]
    leng_t <- length(which(risk_gens_vec == UG_list[[t]]))
    type_t <- paste(length(unlist(strsplit(UG_list[[t]], ';'))), 'genes', sep = '_')
    genes <- UG_list[[t]]
    GS_out <- c(type_t, genes, leng_t, pos_genes)
    return(GS_out)
  }
  coef_list <- sapply(1: k, ScreenKeyGS, family = family)
  risk_genes_list <- coef_list[seq(1, by = 4, length = k)]
  risk_gens_vec <- do.call(rbind, risk_genes_list)
  UG_list <- unique(risk_genes_list)
  
  t = 1: length(unique(risk_genes_list))
  GS_list <- lapply(t, CountGS, risk_gens_vec, UG_list)
  GS_mat <- data.frame(do.call(rbind, GS_list), stringsAsFactors = FALSE)
  colnames(GS_mat) <- c('type', 'genes', 'freuency', 'pos')
  GS_mat$type <- paste(GS_mat$pos, GS_mat$type, sep = '_')
  GS_mat$freuency <- as.numeric(GS_mat$freuency)
  GS_mat$pos <- as.numeric(GS_mat$pos)
  write.table(GS_mat, file = 'modle_gene_counts.txt', row.names = TRUE, col.names = TRUE, sep = "\t")
  #get risk genes, risk genes expression table, risk genes' coefficient value
  
  GS_mat_d <- GS_mat[!(grep('genes', GS_mat$type) %in% grep('_0_genes', GS_mat$type)), ]
  print(GS_mat_d$freuency)
  risk_genes <- risk_genes_list[as.numeric(GS_mat_d$pos[unique(which(GS_mat_d$freuency == max(GS_mat_d$freuency)))])]
  x_coef_lsit <- coef_list[seq(2, by = 4, length = k)]
  x_coef <- x_coef_lsit[as.numeric(GS_mat_d$pos[which(GS_mat_d$freuency == max(GS_mat_d$freuency))])]
  coef_value_d_list <- coef_list[seq(3, by = 4, length = k)]
  coef_value_d <- coef_value_d_list[as.numeric(GS_mat_d$pos[which(GS_mat_d$freuency == max(GS_mat_d$freuency))])]
  lambda_list <- coef_list[seq(4, by = 4, length = k)]
  lambda_d <- lambda_list[as.numeric(GS_mat_d$pos[which(GS_mat_d$freuency == max(GS_mat_d$freuency))])]
  
  #do boxplot
  pdf(file = paste("fren_boxplot.pdf", sep = ""), width = width, height = height)
  p <- ggplot(GS_mat, aes(x = reorder(type, freuency) , y = freuency, fill = type)) +
    geom_bar(stat = "identity")  + labs(title = 'the frequency of models') + ylab("frequency")  + xlab('type') +
    geom_text(aes(label = freuency),vjust = 1, hjust = 0.5, color = 'white', size = 8)
  print(p)
  dev.off()
  
  x_coef_value_list <- list(risk_genes, x_coef, coef_value_d, GS_mat, lambda_d)
  return(x_coef_value_list)
}

uvm_count <- function(unam, rt_uvm, feature = c('OS', 'DFS')){
  #rt_uvm contain OS_Time, OS_Status and other element 
  if(feature == 'OS'){
    t <- rt_uvm[, "OS_Time"]
    d <- rt_uvm[, "OS_Status"]
  } else if (feature == 'DFS'){
    t <- rt_uvm[, "DFS_Time"]
    d <- rt_uvm[, "DFS_Status"]
  }
  y <- Surv(t, d)
  
  fc = coxph(y ~ rt_uvm[, unam], data = rt_uvm)
  tt <- summary(fc)
  uvm_out <- data.frame(cbind(round(tt$coefficients, digits = 3),  round(tt$conf.int, digits = 3)), stringsAsFactors = FALSE)[, c("coef", "exp.coef.", "lower..95", "upper..95", "Pr...z.." )]
  colnames(uvm_out) <- c("coef", "exp_coef", "lower_95%CI", "upper_95%CI", "pvalue")
  uvm_out <- data.frame(cbind(unam, uvm_out), stringsAsFactors = FALSE)
  colnames(uvm_out)[1] <- 'variable'
  return(uvm_out)
}

mvm_count <- function(con_nam, log_nam, tnam, rt_esc){
  t <- rt_esc[, "OS_Time"]
  d <- rt_esc[, "OS_Status"]
  y <- Surv(t, d)
  
  CountCoxph <- function(con_nam, log_nam, tnam, rt_esc){
    if(is.na(log_nam)){
      xc_con = rt_esc[, con_nam]
      ### can rescale you parameter, when you data have NA or NULL, you can use it
      fc = coxph(y ~ ., data = xc_con)
      return(fc)
    } else if (is.na(con_nam)){
      xc_con = rt_esc[, log_nam]
      fc = coxph(y ~ ., data = xc_con) 
      return(fc)
    } else {
      if(length(con_nam) > 1){
        xc_con = data.frame(apply(rt_esc[, con_nam], 2, as.numeric), stringsAsFactors = FALSE)
      } else {
        xc_con = as.numeric(rt_esc[, con_nam])
      }
      #xc_con = covariates(x.con = xc_con, con.rescale = T, fill.missing = T)
      xc_log = rt_esc[, log_nam]
      #xc_log = covariates(x.cat = xc_log, con.rescale = F, fill.missing = T)
      
      ### can rescale you parameter, when you data have NA or NULL, you can use it
      xc = cbind(xc_con, xc_log)
      fc = coxph(y ~ ., data = xc)
      return(fc)
    }
  }
  
  fc <- CountCoxph(con_nam, log_nam, tnam, rt_esc)
  tt <- summary(fc)
  
  mvm_out <- data.frame(cbind(round(tt$coefficients, digits = 3), round(tt$conf.int, digits = 3)), stringsAsFactors = FALSE)[, c("coef", "exp.coef.", "lower..95", "upper..95", "Pr...z.." )]
  colnames(mvm_out) <- c("coef", "HR", "lower_95%CI", "upper_95%CI", "pvalue")
  mvm_out <-data.frame(cbind(row.names(mvm_out), mvm_out), stringsAsFactors = FALSE)
  colnames(mvm_out)[1] <- 'variable'  
  return(mvm_out)
  
  if(all(!is.na(fc$coefficients))){
    pdf(file = paste(tnam, "_cox.pdf",  sep = ""), width = 10, height = 10)
    par(cex = 3)
    plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE,
            cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
    dev.off()
  } else {
    print('the results of some elements is na')
  }
  
}

MakSum <- function(i, mat_coef, coef_value){
  #i is the length of coef_value
  #the colnames of mat_coef is gene names and the rownames is samples id 
  #coef_value is the coefficient value of each  genes
  
  coef_exp <- mat_coef[, 1]*coef_value[1]
  for(i in 2 : i){
    coef_expi <- mat_coef[, i]*coef_value[i]
    coef_exp <- coef_exp + coef_expi
  }
  return(coef_exp)
}

DrawSurminer <- function(Gnam, Rt_Exp_Cli, DatType = c("ConType", "LogType"), theme_surv, len_a = len_a, len_b = len_b, s_cutoff = NULL, 
                         palette = palette, risk.table = risk.table, fun = NULL, rt_height = 0.3){

  print("your survival data title must be OS_Status and OS_Time")
  print("len_a < median")
  Gnam <<- Gnam
  Rt_Exp_Cli <<- Rt_Exp_Cli
  
  PreSurList <- function(Gnam, Rt_Exp_Cli){

    
    Rt_Exp_Cli$OS_Time <- as.numeric(as.character(Rt_Exp_Cli$OS_Time))
    Rt_Exp_Cli$OS_Status <- as.numeric(as.character(Rt_Exp_Cli$OS_Status))
    diff <- survdiff(Surv(OS_Time, OS_Status) ~ Rt_Exp_Cli[, Gnam], data = Rt_Exp_Cli)
    pValue <- 1-pchisq(diff$chisq, df=1)
    fit <- survfit(Surv(OS_Time, OS_Status) ~ Rt_Exp_Cli[, Gnam], data = Rt_Exp_Cli)
    pM_List <- list(pValue = pValue, fit = fit)
    return(pM_List)
  }
  
  MakSurList <- function(Gnam, Rt_Exp_Cli, DatType = c("ConType", "LogType"), len_a = len_a, len_b = len_b, s_cutoff = NULL){
    
    if(length(which((c("DFS_Time", "DFS_Status", "OS_Time", "OS_Status") %in% colnames(Rt_Exp_Cli))  == TRUE)) >1){
      if(DatType == "ConType" & length(unique(Rt_Exp_Cli[, Gnam])) < 3){
        stop("if yout data type is continuous variable, your unique elements should be over 2")
      }
      if(DatType == "ConType" & length(unique(Rt_Exp_Cli[, Gnam])) > 3){
        if(is.null(s_cutoff)){
          Rt_Exp_Cli[, Gnam][Rt_Exp_Cli[, Gnam] > median(Rt_Exp_Cli[, Gnam])] <- len_b
          Rt_Exp_Cli[, Gnam][which(Rt_Exp_Cli[, Gnam] != len_b)] <- len_a
          pM_List <- PreSurList(Gnam, Rt_Exp_Cli)
          return(pM_List)  
        } else if (!is.null(s_cutoff)){
          Rt_Exp_Cli[, Gnam][Rt_Exp_Cli[, Gnam] > s_cutoff]  <- len_b
          Rt_Exp_Cli[, Gnam][which(Rt_Exp_Cli[, Gnam] != len_b)] <- len_a
          pM_List <- PreSurList(Gnam, Rt_Exp_Cli)
          return(pM_List)  
        } 
      }
      if(DatType == "LogType" & length(unique(Rt_Exp_Cli[, Gnam])) < 2){
        stop("if yout data type is continuous variable, your unique elements should be over  2")
      } else if(DatType == "LogType" & length(unique(Rt_Exp_Cli[, Gnam])) >= 2){
        pM_List <- PreSurList(Gnam, Rt_Exp_Cli)
        return(pM_List)
      }
    }else{
      stop("please change your colnames which should contain DFS_Time and DFS_Status or OS_Time and OS_Status")
    }
  }
  
  pM_List <- MakSurList(Gnam, Rt_Exp_Cli, DatType, len_a = len_a, len_b = len_b, s_cutoff = s_cutoff)
  fit <- pM_List$fit
  print(fit)
  pValueR <- pM_List$pValue
  print(names(fit$strata))
  names(fit$strata) <- gsub("Rt_Exp_Cli., Gnam]=", "", names(fit$strata))
  
  if(is.null(fun)){
    pdf(file = paste(Gnam, ".pdf", sep = ""), 11, 11, onefile = FALSE)
    dt <- ggsurvplot(fit, title = Gnam, pval = paste('P = ', format(pValueR, digits = 4, scientific = TRUE), sep = ""), 
                     conf.int = FALSE, pval.size = 10, pval.coord = c(1, 0), xlab = "Time in days", 
                     size = 3, #change line size
                     risk.table = risk.table, # Add risk table
                     palette = palette, # custom color palette
                     risk.table.col = "strata", # Change risk table color by groups
                     ggtheme = theme_surv, fontsize = 8, risk.table.height = rt_height, fun = fun)
    print(dt)
    dev.off()
    
  } else {
    pdf(file = paste(Gnam, ".pdf", sep = ""), 11, 11, onefile = FALSE)
    dt <- ggsurvplot(fit,  title = Gnam, pval = paste('P = ', format(pValueR, digits = 4, scientific = TRUE), sep = ""), 
                     conf.int = FALSE, pval.size = 10, pval.coord = c(1, 1), xlab = "Time in days", 
                     size = 3, #change line size
                     risk.table = risk.table, # Add risk table
                     palette = palette, # custom color palette
                     risk.table.col = "strata", # Change risk table color by groups
                     ggtheme = theme_surv, fontsize = 8, risk.table.height = rt_height, fun = fun)
    print(dt)
    dev.off()
  }
  
  rm(Gnam, Rt_Exp_Cli, pos = ".GlobalEnv")
}

ConSamMat <- function(rt_exp, grp_nam, len_a = len_a, len_b = len_b, col_a = col_a, col_b = col_b, g_cutoff = NULL){
  if(is.null(g_cutoff)){
    rt_exp$group[rt_exp[, grp_nam] > median(rt_exp[, grp_nam])] <- len_b
    rt_exp$group[rt_exp[, grp_nam] <= median(rt_exp[, grp_nam])] <- len_a
  } else {
    rt_exp$group[rt_exp[, grp_nam] > g_cutoff] <- len_b
    rt_exp$group[rt_exp[, grp_nam] <= g_cutoff] <- len_a
  }
  
  rt_exp$color <- NA
  rt_exp$color[grep(len_a, rt_exp$group)] <- col_a
  rt_exp$color[grep(len_b, rt_exp$group)] <- col_b
  
  rt_sam <- data.frame(cbind(row.names(rt_exp), rt_exp[, c('group', 'color')]), stringsAsFactors = FALSE)
  colnames(rt_sam) <- c('samples_id', 'group', 'color')
  
  return(rt_sam)
}

PlotCliComHeat <- function(exp_mat, rt_cli_com, cancer, Fnam){
  if(cancer == "KIRP"){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            stage = rt_cli_com$stage, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       stage = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 8)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(10.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
  if(cancer == 'LAML'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(gender = rt_cli_com$gender, AML_FAB_subtype = rt_cli_com$AML_FAB_subtype, cytogenetic_risk_group = rt_cli_com$cytogenetic_risk_group,
                            cytogenetic_abnormality = rt_cli_com$cytogenetic_abnormality, hemoglobin = rt_cli_com$hemoglobin,  blast_cell_BM = rt_cli_com$blast_cell_BM,
                            leukocyte = rt_cli_com$leukocyte, platelet = rt_cli_com$platelet, blast_cell_PB = rt_cli_com$blast_cell_PB, 
                            col = list(gender = structure(names = c("1", "2"), c( 'black', 'grey')), 
                                       AML_FAB_subtype = structure(names = c("1", "2", "3", "4", "5", "6", "7", "8"), c("green4", "gold3", "orange1", "firebrick2", "plum3", "gold2", "darkred",  "yellow2")), 
                                       cytogenetic_risk_group = structure(names = c("1", "2", "3"), c("grey", "black", "white")),
                                       cytogenetic_abnormality = structure(names = c("1", "2", "3"), c("grey", "black", "white")), 
                                       hemoglobin = colorRamp2(c(0, max(rt_cli_com$hemoglobin)), c("white", "red")), 
                                       blast_cell_BM = colorRamp2(c(0, max(rt_cli_com$blast_cell_BM)), c("white", "red")), 
                                       leukocyte = colorRamp2(c(0, max(rt_cli_com$leukocyte)), c("white", "red")), 
                                       platelet = colorRamp2(c(0, max(rt_cli_com$platelet)), c("white", "red")), 
                                       blast_cell_PB = colorRamp2(c(0, max(rt_cli_com$blast_cell_PB)), c("white", "red"))),
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 5, 5, 5, 5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 10)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(13.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
  if(cancer == 'LUAD'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            stage = rt_cli_com$stage, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       stage = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 8)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(10.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
  if(cancer == 'LIHC'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            stage = rt_cli_com$stage, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       stage = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(2.5, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(2.5, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(2.5, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 25, 25, 25), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 9)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)/1.2), c(" blue", "white", "firebrick2")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(9, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
    
  }
  if(cancer == 'READ'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            stage = rt_cli_com$stage, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       stage = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 8)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(10.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
  if(cancer == 'SARC'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            tumor_depth = rt_cli_com$tumor_depth, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       tumor_depth = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 8)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(10.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
  if(cancer == 'KIRC'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            stage = rt_cli_com$stage, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       stage = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 8)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(10.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
}

MakPCA <- function(rt_pca, rt_sam = rt_sam, target_genes = NULL, nam = nam, phi = 40, pch = 18, len_a = len_a, len_b = len_b,
                   len_position = 'bottom', width = NULL, height = NULL){
  if(is.null(target_genes)){
    rt_pca <- rt_pca
  }else if (!is.null(target_genes)){
    rt_pca <- rt_pca[, target_genes]
  }
  
  data.a <- as.matrix(rt_pca)
  data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
  a <- summary(data.pca)
  tmp <- a$importance  # a include 4 sections which contain importance
  pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100
  pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100
  pro3 <- as.numeric(sprintf("%.3f", tmp[2,3]))*100# fetch the proportion of PC1 and PC2
  pc <- as.data.frame(a$x)  # convert to data.frame
  num_pos <- match(row.names(pc), rt_sam$samples_id, nomatch = NA)
  
  pc$group <-  rt_sam$group[num_pos]
  pc$color <-  rt_sam$color[num_pos]
  pc$group <- factor(pc$group, levels = unique(pc$group))
  
  xlab <- paste("PC1(", pro1, "%)", sep = "")
  ylab <- paste("PC2(", pro2, "%)", sep = "")
  zlab <- paste("PC3(", pro3, "%)", sep = "")
  
  pdf(file = paste(nam, "PCA.pdf", sep = ""), width, height)
  PCA <- scatter3D(pc$PC1, pc$PC2, pc$PC3,  xlab = xlab, ylab = ylab, zlab = zlab, bty = 'g', cex = 1.5,
                   pch  = pch, phi = phi, col = as.character(unique(pc$color)),  colvar = as.integer(pc$group), colkey = FALSE)
  legend(len_position, legend = levels(pc$group), col = as.character(unique(pc$color)), pch = 18, inset = -0.17, xpd = TRUE, horiz = TRUE, cex = 1.5)
  print(PCA)
  dev.off()
}

MakBarPlot <- function(mat_ht = mat_ht, Vnum = Vnum, Hnum = Hnum, Vnam = NULL, Hnam = NULL, Tnam = NULL, heigth = NULL, width = NULL,
                       fill_col = "steelblue", themeP = NULL, yaesmin = yaesmin, yaesmax = yaesmax, add_text = FALSE){
  if(add_text == TRUE){
    pdf(file = paste(Hnum, ".pdf", sep = ""), height = heigth, width = width)
    p = ggplot(mat_ht, aes(x = get(Vnum), y = get(Hnum))) +
      geom_bar(stat = "identity", fill = fill_col) + coord_cartesian(ylim = c(yaesmin, yaesmax)) +
      labs(x = Vnam, y = Hnam, title = Tnam) +  geom_text(aes(label = get(Hnum)), color = "black") +
      themeP
    print(p)
    dev.off()
  } else {
    pdf(file = paste(Hnum, ".pdf", sep = ""), height = heigth, width = width)
    p = ggplot(mat_ht, aes(x = get(Vnum), y = get(Hnum))) +
      geom_bar(stat = "identity", fill = fill_col) + coord_cartesian(ylim = c(yaesmin, yaesmax)) +
      labs(x = Vnam, y = Hnam, title = Tnam) +
      themeP
    print(p)
    dev.off()
  }
}

MakBoxPlot <- function(gnam, gggroup, rt_box_mat, width = 4, height = 6, theme_VB, ggtype = c('boxplot', 'violin'), ggcolor = ggcolor, ggylab = ggylab){
  gnam <<- gnam
  gggroup <<- gggroup
  rt_box <- rt_box_mat[, c(gnam, gggroup)]
  rt_box[, gnam] <- as.numeric(rt_box[, gnam])
  colnames(rt_box) <- c(gnam, "group")
  
  if(length(unique(rt_box$group)) > 2){
    print("ANOVA will be used becuase your group > 2")
    res_aov <- aov(rt_box[, gnam] ~ group, data = rt_box)
    Pval <- round(unlist(summary(res_aov))[9], digits = 4)
    
  } else if (length(unique(rt_box$group)) == 2){
    #do t.test
    print("t.test will be used becuase your group = 2")
    GPx <- rt_box[, gnam][grep(unique(rt_box$group)[1], rt_box$group)]
    GPy <- rt_box[, gnam][grep(unique(rt_box$group)[2], rt_box$group)]
    Obj <- t.test(GPx, GPy, var.equal = FALSE, paired = FALSE)
    Pval <- round(Obj$p.value, digits = 4)
  } else {
    stop("your group < 2")
  }
  ###do boxplot
  if(ggtype == "boxplot"){
    pdf(file = paste(gnam, '_', gggroup,"_boxplot_point.pdf", sep = ""), width = width, height = height)
    p <- ggplot(rt_box, aes(x = group, y = rt_box[, gnam], fill = group)) + geom_boxplot() +
      scale_fill_manual(values= ggcolor) + annotate("text", x = -Inf, y = -Inf, hjust = 0, vjust = 0, label= paste("P = ", Pval, sep = "")) +
      geom_jitter(shape=16, position = position_jitter(0.2)) +  labs(title = gnam) + xlab(gggroup) + ylab(ggylab) +
      theme_VB
    print(p)
    dev.off()
    
    #no point
    pdf(file = paste(gnam, '_', gggroup, "_boxplot.pdf", sep = ""), width = width, height = height)
    p <- ggplot(rt_box, aes(x = group, y = rt_box[, gnam], fill = group)) + geom_boxplot() +
      scale_fill_manual(values= ggcolor) + annotate("text",  x = -Inf, y = -Inf, hjust = 0, vjust = 0, label= paste("P = ", Pval, sep = ""))  +
      labs(title = gnam) + xlab(gggroup) + ylab(ggylab) + theme_VB
    print(p)
    dev.off()
    rm(gnam, gggroup, pos = ".GlobalEnv")
  } else if (ggtype == "violin"){
    pdf(file = paste(gnam, '_', gggroup, "_violin.pdf", sep = ""), width = width, height = width)
    p <- ggplot(rt_box, aes(x = group, y = rt_box[, gnam])) +  geom_violin(aes(fill = group), trim = FALSE) + geom_boxplot(width = 0.1) +
      scale_fill_manual(values = ggcolor)+ annotate("text", x = -Inf, y = -Inf, hjust = 0, vjust = 0, label = paste("P = ", Pval, sep = "")) +
      labs(title = gnam)  + xlab(gggroup) + ylab(ggylab) + theme_VB
    print(p)
    dev.off()
    rm(gnam, gggroup, pos = ".GlobalEnv")
  } else {
    stop('you must set ggtype')
  }
}

theme_E <- theme_bw()+
  theme(plot.title = element_text(size = 20, face="bold", hjust = 0.5),
        axis.title.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 90),
        axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        legend.title = element_text(colour = 'black', angle = 1, size = 15, face = 'bold'),
        legend.text = element_text(colour = 'black', angle = 1, size = 15, face = 'bold'),
        legend.position = "top",
        strip.text.x = element_text(size = 15, face="bold"),
        strip.text.y = element_text(size = 15, face="bold"),
        panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 1.5),
        panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 1.5))

StratifySur <- function(rt_cli_mat, sur_nam = sur_nam, pheno_ele = pheno_ele, strat_ele = strat_ele, title_nam = title_nam, theme_surv, 
                        DatType = c('ConType', 'LogType'), len_a = len_a, len_b = len_b, palette = c('red', 'blue'), risk.table = TRUE, fun = NULL){
  
  print("pheno_ele like stage/TMN, sur_nam like gene or risk score, title_nam like risk_score_M0")
  
  rt_cli_mat_strat <- rt_cli_mat[which(rt_cli_mat[, pheno_ele] == strat_ele), ]
  colnames(rt_cli_mat_strat)[which(colnames(rt_cli_mat_strat) == sur_nam)] <- title_nam
  
  str_sur <- DrawSurminer(title_nam, rt_cli_mat_strat, DatType = DatType, theme_surv, len_a = len_a, len_b = len_b, 
                          s_cutoff = NULL,  palette = palette, risk.table = risk.table, fun = fun)
}

theme_surv_a <- theme_survminer(base_size = 30, font.main = c(35, "bold", "black"), font.submain = c(35, "bold", "black"), 
                                font.x = c(30,"bold", "black"), font.y = c(30, "bold", "black"), font.caption = c(30, "bold", "black"), 
                                font.tickslab = c(30, "bold", "black"), 
                                legend ="top", font.legend = c(25, "bold", "black"))


rt_exp <- read_tsv(file = "https://gdc.xenahubs.net/download/TCGA-KIRC.htseq_fpkm.tsv.gz")
rt_cli <- read_tsv(file = "https://gdc.xenahubs.net/download/TCGA-KIRC.GDC_phenotype.tsv.gz")
rt_sur <- read_tsv(file = "https://gdc.xenahubs.net/download/TCGA-KIRC.survival.tsv.gz")
rt_T <- split_tcga_tn(rt_exp[, -1], sam_type = "tumor")
cancer_list <- GetExpSurCli(rt_T, rt_cli, rt_sur)
rt_T_m <- cancer_list[[1]]
rt_cli_m <- cancer_list[[2]]
rt_sur_m <- cancer_list[[3]]
cancer_cli_sort <- CliSort(rt_cli_m, cancer)
rt_cli_m <- cancer_cli_sort

#
exp_d_pos <- apply(rt_T_m, 1, function(x){median(x) > 0})
unam <- rt_exp$Ensembl_ID[exp_d_pos]
rt_cancer <- data.frame(cbind(t(rt_T_m), rt_sur_m[, c("X_OS_IND", "X_OS")]), stringsAsFactors = FALSE)
colnames(rt_cancer) <- c(unam, "OS_Status", "OS_Time")

uvm_list <- lapply(unam, uvm_count, rt_cancer)
data_uvm <- do.call(rbind, uvm_list)
sig_genes <-  row.names(data_uvm)[data_uvm$pvalue < 0.05]
rt_mvm <- rt_T_m[match(sig_genes, rt_exp$Ensembl_ID, nomatch = 0), ]
row.names(rt_mvm) <- sig_genes
row.names(rt_mvm) <- matrix(unlist(strsplit(row.names(rt_mvm), '[.]')), ncol = 2, byrow = TRUE)[, 1]

rt_mvm_sur <- data.frame(cbind(t(rt_mvm), rt_sur_m[, c("X_OS_IND", "X_OS")]), stringsAsFactors = FALSE)
rt_immune <- read.table(file = "GOappend.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
eg_imm <-bitr(rt_immune$Symbol, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = 'org.Hs.eg.db')
rt_mvm_sur_imm <- rt_mvm_sur[, match(c(eg_imm$ENSEMBL, "X_OS_IND", "X_OS"), colnames(rt_mvm_sur), nomatch = 0)]
set.seed(1234)
ind <- sample(2, nrow(rt_mvm_sur_imm), replace = TRUE, prob = c(0.6, 0.4))
ind <- ind_mat$x
rt_mvm_sur_train <- rt_mvm_sur_imm[ind == 1, ] 
rt_mvm_sur_test <- rt_mvm_sur_imm[ind == 2, ]

#
x <- data.matrix(rt_mvm_sur_train[, -c((length(colnames(rt_mvm_sur_train))-1):length(colnames(rt_mvm_sur_train)))])
y <- data.matrix(rt_mvm_sur_train[, c('X_OS', 'X_OS_IND')])
colnames(y) <- c("time", "status")
x_coef_value_list <- GLCoxMain(1000, x, y, s = 'lambda.min', family = "cox", 7, 6)

risk_genes <- strsplit(x_coef_value_list[1][[1]][[1]], ';')[[1]]
x_coef <- x_coef_value_list[2][[1]][[1]]
coef_value_d <- x_coef_value_list[3][[1]][[1]]
GC_mat <- x_coef_value_list[4][[1]]
lambda_d <- x_coef_value_list[5][[1]]

risk_genes_coef_mat <- data.frame(cbind(risk_genes, coef_value_d), stringsAsFactors = FALSE)
colnames(risk_genes_coef_mat) <- c('risk_gens', 'coef_value')
risk_genes <- risk_genes_coef_mat$risk_gens
coef_value_d <- risk_genes_coef_mat$coef_value
x_coef <- read.table(file = paste("/Users/stead/Desktop/subtype_analysis/signature/", cancer, "/", cancer, "_", num, "/", cancer, "_risk_genes_exp.txt", sep = ""),
                     sep = "\t", stringsAsFactors = FALSE)
eg_sym <- bitr(risk_genes, fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db', drop = FALSE)
risk_genes_sym <- eg_sym[, 'SYMBOL'][!duplicated(eg_sym$ENSEMBL)]

###
coef_score_train <- MakSum(length(coef_value_d), x_coef, coef_value_d)
rt_mvm_sur_train_min <- rt_mvm_sur_train[, c(risk_genes, 'X_OS', 'X_OS_IND')]
colnames(rt_mvm_sur_train_min) <- c(risk_genes_sym, 'X_OS', 'X_OS_IND')
rt_mvm_sur_train_min$score <- coef_score_train
colnames(rt_mvm_sur_train_min)[(length(colnames(rt_mvm_sur_train_min)) -2) : length(colnames(rt_mvm_sur_train_min))] <- c("OS_Time", "OS_Status", 'risk_score')

DrawSurminer('risk_score', rt_mvm_sur_train_min, DatType = "ConType", theme_surv_a, len_a = 'low_risk', len_b = 'high_risk',
             s_cutoff = NULL,  palette = c('red', 'blue'), risk.table = TRUE, fun = NULL)

###
rt_cli_m_train <-rt_cli_m[match(row.names(rt_mvm_sur_train_min), rt_cli_m$sample_id, nomatch = 0), ]
rt_mvm_sur_train_min$stage <- rt_cli_m_train$stage
rt_mvm_sur_train_min$gender <- rt_cli_m_train$gender
rt_mvm_sur_train_min$recurrence <- rt_cli_m_train$recurrence
rt_mvm_sur_train_min$age <- rt_cli_m_train$age

order_score_train <- order(rt_mvm_sur_train_min$risk_score)
rt_mvm_sur_train_order <- rt_mvm_sur_train_min[order_score_train, ]
rt_mvm_sur_train_order[, c('OS_Time', 'risk_score', 'age')] <- apply(rt_mvm_sur_train_order[, c('OS_Time', 'risk_score', 'age')], 2, as.numeric)
exp_train <- apply(rt_mvm_sur_train_order[, risk_genes_sym], 1, FUN = function(x){(x-mean(x))/sd(x)})
row.names(exp_train) <- risk_genes_sym
PlotCliComHeat(exp_train, rt_mvm_sur_train_order, 'KIRC', "train")

rt_sam <- ConSamMat(rt_mvm_sur_train_order, 'risk_score', len_a = 'low_risk', len_b = 'high_risk', col_a = 'blue', col_b = 'red', g_cutoff = NULL)
MakPCA(rt_mvm_sur_train_order, rt_sam = rt_sam, risk_genes_sym, nam = 'training', phi = 60, pch = 18, len_a = 'low_risk', len_b = 'high_risk',
       len_position = 'bottom', width = NULL, height = NULL)

###
i = length(coef_value_d)
x_coef_test <- rt_mvm_sur_test[, risk_genes]
coef_score_test <- MakSum(i, x_coef_test, coef_value_d)

rt_mvm_sur_test_min <- rt_mvm_sur_test[, c(risk_genes, 'X_OS', 'X_OS_IND')]
colnames(rt_mvm_sur_test_min) <- c(risk_genes_sym, 'X_OS', 'X_OS_IND')

rt_mvm_sur_test_min$score <- coef_score_test
colnames(rt_mvm_sur_test_min)[(length(colnames(rt_mvm_sur_test_min)) -2) : length(colnames(rt_mvm_sur_test_min))] <- c("OS_Time", "OS_Status", 'risk_score')

DrawSurminer('risk_score', rt_mvm_sur_test_min, DatType = "ConType", theme_surv_a, len_a = 'low_risk', len_b = 'high_risk',
             s_cutoff = NULL,  palette = c('red', 'blue'), risk.table = TRUE, fun = NULL)

cli_pos_test <- match(row.names(rt_mvm_sur_test_min), rt_cli_m$sample_id, nomatch = 0)
rt_mvm_sur_cli_test_min <- data.frame(cbind(rt_mvm_sur_test_min, rt_cli_m[cli_pos_test, -1]), stringsAsFactors = FALSE)

order_score_test <- order(rt_mvm_sur_cli_test_min$risk_score)
rt_mvm_sur_test_order <- rt_mvm_sur_cli_test_min[order_score_test, ]
rt_mvm_sur_test_order[, c('OS_Time', 'risk_score', 'age')] <- apply(rt_mvm_sur_test_order[, c('OS_Time', 'risk_score', 'age')], 2, as.numeric)
rt_mvm_sur_test_order$age[is.na(rt_mvm_sur_test_order$age)] <- 0
exp_test <- apply(rt_mvm_sur_test_order[, 1: length(risk_genes_sym)], 1, scale)
row.names(exp_test) <- risk_genes_sym
PlotCliComHeat(exp_test, rt_mvm_sur_test_order, 'KIRC', "test")

###
i = length(coef_value_d)
x_coef_all <- rt_mvm_sur[, risk_genes]
coef_score_all <- MakSum(i, x_coef_all, coef_value_d)

rt_mvm_sur_all_min <- rt_mvm_sur[, c(risk_genes, 'X_OS', 'X_OS_IND')]
colnames(rt_mvm_sur_all_min) <- c(risk_genes_sym, 'X_OS', 'X_OS_IND')

rt_mvm_sur_all_min$score <- coef_score_all
colnames(rt_mvm_sur_all_min)[(length(colnames(rt_mvm_sur_all_min)) -2) : length(colnames(rt_mvm_sur_all_min))] <- c("OS_Time", "OS_Status", 'risk_score')

DrawSurminer('risk_score', rt_mvm_sur_all_min, DatType = "ConType", theme_surv_a, len_a = 'low_risk', len_b = 'high_risk',
             s_cutoff = NULL,  palette = c('red', 'blue'), risk.table = TRUE, fun = NULL)

cli_pos_all <- match(row.names(rt_mvm_sur_all_min), rt_cli_m$sample_id, nomatch = 0)
rt_mvm_sur_cli_all_min <- data.frame(cbind(rt_mvm_sur_all_min, rt_cli_m[cli_pos_all, ]), stringsAsFactors = FALSE)

###do complex heatmap
order_score_all <- order(rt_mvm_sur_cli_all_min$risk_score)
rt_mvm_sur_all_order <- rt_mvm_sur_cli_all_min[order_score_all, ]
rt_mvm_sur_all_order[, c('OS_Time', 'risk_score', 'age')] <- apply(rt_mvm_sur_all_order[, c('OS_Time', 'risk_score', 'age')], 2, as.numeric)
exp_all <- apply(rt_mvm_sur_all_order[, risk_genes_sym], 1, scale)
row.names(exp_all) <- risk_genes_sym

PlotCliComHeat(exp_all, rt_mvm_sur_all_order, 'KIRC', "all")

rt_cli_box <- rt_mvm_sur_all_order[, c("gender", "T", "N",  "M", "stage" , "recurrence", "risk_score", "OS_Time", "OS_Status")]
for(gnam in colnames(rt_cli_box)[1:6]){
  rt_cli_box_g <- rt_cli_box[, c(gnam, 'risk_score')]
  colnames(rt_cli_box_g) <- c('group', gnam)
  MakBoxPlot(gnam, 'group', rt_cli_box_g, width = 4, height = 6, theme_E, ggtype = 'boxplot', ggcolor = color_type[1 : length(unique(rt_cli_box_g[, gnam]))], ggylab = 'risk_score')
}

