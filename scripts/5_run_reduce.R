#fifth script to run batchtools for imputation quality comparison

library(data.table)
library(batchtools)

source(file.path(source_dir, "result_source.R"))




#makes new registry for reducing data for chromosome 22
reg_reduce_22 <- makeRegistry(
  file.dir = file.path(dir_registry, "stahl_reduce_22"),
  make.default = FALSE,
  source = file.path(source_dir, "result_source.R"),
  conf.file = confic_dir
)
#makes new registry for reducing data for chromosome 19
reg_reduce_19 <- makeRegistry(
  file.dir = file.path(dir_registry, "stahl_reduce_19"),
  make.default = FALSE,
  source = file.path(source_dir,"result_source.R"),
  conf.file = confic_dir
)
#makes new registry for reducing data for chromosome 19 imputed with pbwt
reg_reduce_pbwt_19 <- makeRegistry(
  file.dir = file.path(dir_registry, "stahl_reduce_pbwt_19"),
  make.default = FALSE,
  source = file.path(source_dir,"result_source.R"),
  conf.file = confic_dir
)




#load registry of imputation results
reg<-loadRegistry(file.path(dir_registry,"stahl_phasing_imputation"), writeable=TRUE)
#load registry for imputation results for pbwt on chromosome 19
reg_pbwt<-loadRegistry(file.path(dir_registry,"stahl_imputation_pbwt_19"),
                       writeable=TRUE,make.default = FALSE)



#load registry for batchmapresults for chrom22
tmp_results_22<-loadRegistry(file.path(dir_registry,"stahl_batchmapresults_22"),
                             writeable=TRUE,make.default = FALSE)
#load registry for batchmapresults for chrom19
tmp_results_19<-loadRegistry(file.path(dir_registry,"stahl_batchmapresults_19"),
                             writeable=TRUE,make.default = FALSE)
#load_regi for batchmapresults for chrom19 pbwt
tmp_results_pbwt_19<-loadRegistry(file.path(dir_registry,"stahl_batchmapresults_19"),
                                  writeable=TRUE,make.default = FALSE)

#split data according to chromosome 19 and 22 jobs, pbwt on chromosome 19 separate
job_table <- batchtools::unwrap(getJobTable(reg = reg), cols = "prob.pars")
job_table_pbwt_19 <- batchtools::unwrap(getJobTable(reg = reg_pbwt), cols = "prob.pars")


job_table_results_22 <- batchtools::unwrap(getJobTable(reg = tmp_results_22))
job_table_22 <- job_table[job_table_results_22, on = c("job.id" = ".id")]

job_table_results_19 <- batchtools::unwrap(getJobTable(reg = tmp_results_19))
job_table_19 <- job_table[job_table_results_19, on = c("job.id" = ".id")]

job_table_results_pbwt_19 <- batchtools::unwrap(getJobTable(reg = tmp_results_pbwt_19))
job_table_pwbt_19 <- job_table_pbwt_19[job_table_results_pbwt_19, on = c("job.id" = ".id")]

#split data according to phasing/imputation combination for chromosome 19
split_list_job_tables_19 <-split(job_table_19, by=c("problem", "algorithm","ref.on"))
split_list_job_id_19<-lapply(split_list_job_tables_19, function(x) x$i.job.id)
#split data according to phasing/imputation combination for chromosome 22
split_list_job_tables_22 <-split(job_table_22, by=c("problem", "algorithm","ref.on"))
split_list_job_id_22<-lapply(split_list_job_tables_22, function(x) x$i.job.id)
#split data according to phasing for imputation of chromosome 19 with pbwt
split_list_job_tables_pbwt_19 <-split(job_table_pbwt_19, by="which.phasing")
split_list_job_id_pbwt_19<-lapply(split_list_job_tables_pbwt_19, function(x) x$job.id)





#define batchmap step with function to reduce data further to SNP basis and calculate different scores for chromosome 22,
#       individual participant data is no longer included after this step
ids_22 <- batchMap(
  fun = function(ids, reg.result) {
    res <- rbindlist(reduceResultsList(ids = ids, reg = reg.result))
    res[, DOSAGE := GP_1 + 2*GP_2]
    res[, SUPERDOSAGE := GP_1 + 2*GP_2]
    n <- res[, length(unique(SAMPLE))]
    res[, MAF_EST := sum(DOSAGE)/(2*n), by = .( POS, REF, ALT)]
    res[, .(
      MIN_H_SCORE = min(H_SCORE, na.rm = TRUE),
      Q25_H_SCORE = quantile(H_SCORE, probs = 0.25, na.rm = TRUE),
      MEDIAN_H_SCORE = median(H_SCORE, na.rm = TRUE),
      MEAN_H_SCORE = mean(H_SCORE, na.rm = TRUE),
      SD_H_SCORE = sd(H_SCORE, na.rm = TRUE),
      Q75_H_SCORE = quantile(H_SCORE, probs = 0.75, na.rm = TRUE), 
      MAX_H_SCORE = max(H_SCORE, na.rm = TRUE),
      MIN_SEN_SCORE = min(SEN_SCORE, na.rm = TRUE),
      Q25_SEN_SCORE = quantile(SEN_SCORE, probs = 0.25, na.rm = TRUE),
      MEDIAN_SEN_SCORE = median(SEN_SCORE, na.rm = TRUE),
      MEAN_SEN_SCORE = mean(SEN_SCORE, na.rm = TRUE),
      SD_SEN_SCORE = sd(SEN_SCORE, na.rm = TRUE),
      Q75_SEN_SCORE = quantile(SEN_SCORE, probs = 0.75, na.rm = TRUE), 
      MAX_SEN_SCORE = max(SEN_SCORE, na.rm = TRUE),
      IQS = iqs(.SD[!is.na(GT)]),
      MACHR2 = (sum(DOSAGE^2)/n-(sum(DOSAGE)/n)^2)/(2*MAF_EST*(1-MAF_EST)),
      BEAGLER2 = (sum(BEST_GUESS*DOSAGE)-1/n*sum(BEST_GUESS)*sum(DOSAGE))^2/((sum(SUPERDOSAGE)-1/n*sum(DOSAGE)^2)*(sum(BEST_GUESS^2)-1/n*sum(BEST_GUESS)^2)),
      IMPUTEINFO = 1-sum(SUPERDOSAGE-DOSAGE^2)/(2*n*MAF_EST*(1-MAF_EST))
    ),
    by = .(IMPUTED, PHASED, POS, REF, ALT, MAF_EST)]
  }, 
  split_list_job_id_22, more.args = list(reg.result = tmp_results_22), 
  reg = reg_reduce_22
)
ids_22[, chunk := 1]


#run batchmap jobs for chromosome 22
submitJobs(
  ids = ids_22,
  reg = reg_reduce_22,
  resources = list(
    ntasks = 1, ncpus = 1, memory = 250000, 
    walltime = 0L,
    partition = "prio",
    account = "dzhkomics",
    chunks.as.arrayjobs = TRUE
  )
)


#convert results of batchmap fpr chromosome 22 to data.table as an R object
reduced_tables_22 <- rbindlist(reduceResultsList(reg = reg_reduce_22), idcol = "job.id")

#add SNP specific info from original data set
variant_info_22 <- data.table::fread(dir_input_data_true_annotated_variant_infos_22, na.strings = ".")
setnames(variant_info_22, c("CHROM", "POS", "REF", "ALT", "VQSLOD", "QD", "InbreedingCoeff", "DP", "MAF"))
variant_info_22[, MAF := 0.5 - abs(MAF - 0.5)]
reduced_tables_22 <- reduced_tables_22[variant_info_22, on = .(POS, REF, ALT), nomatch = 0L]


#save data
#saveRDS(reduced_tables_22, file = dir_out_table_22)




#define batchmap step with function to reduce data further to SNP basis and calculate different scores for chromosome 19
ids_19 <- batchMap(
  fun = function(ids, reg.result) {
    res <- rbindlist(reduceResultsList(ids = ids, reg = reg.result))
    res[, DOSAGE := GP_1 + 2*GP_2]
    res[, SUPERDOSAGE := GP_1 + 2*GP_2]
    n <- res[, length(unique(SAMPLE))]
    res[, MAF_EST := sum(DOSAGE)/(2*n), by = .( POS, REF, ALT)]
    res[, .(
      MIN_H_SCORE = min(H_SCORE, na.rm = TRUE),
      Q25_H_SCORE = quantile(H_SCORE, probs = 0.25, na.rm = TRUE),
      MEDIAN_H_SCORE = median(H_SCORE, na.rm = TRUE),
      MEAN_H_SCORE = mean(H_SCORE, na.rm = TRUE),
      SD_H_SCORE = sd(H_SCORE, na.rm = TRUE),
      Q75_H_SCORE = quantile(H_SCORE, probs = 0.75, na.rm = TRUE), 
      MAX_H_SCORE = max(H_SCORE, na.rm = TRUE),
      MIN_SEN_SCORE = min(SEN_SCORE, na.rm = TRUE),
      Q25_SEN_SCORE = quantile(SEN_SCORE, probs = 0.25, na.rm = TRUE),
      MEDIAN_SEN_SCORE = median(SEN_SCORE, na.rm = TRUE),
      MEAN_SEN_SCORE = mean(SEN_SCORE, na.rm = TRUE),
      SD_SEN_SCORE = sd(SEN_SCORE, na.rm = TRUE),
      Q75_SEN_SCORE = quantile(SEN_SCORE, probs = 0.75, na.rm = TRUE), 
      MAX_SEN_SCORE = max(SEN_SCORE, na.rm = TRUE),
      IQS = iqs(.SD[!is.na(GT)]),
      MACHR2 = (sum(DOSAGE^2)/n-(sum(DOSAGE)/n)^2)/(2*MAF_EST*(1-MAF_EST)),
      BEAGLER2 = (sum(BEST_GUESS*DOSAGE)-1/n*sum(BEST_GUESS)*sum(DOSAGE))^2/((sum(SUPERDOSAGE)-1/n*sum(DOSAGE)^2)*(sum(BEST_GUESS^2)-1/n*sum(BEST_GUESS)^2)),
      IMPUTEINFO = 1-sum(SUPERDOSAGE-DOSAGE^2)/(2*n*MAF_EST*(1-MAF_EST))
    ),
    by = .(IMPUTED, PHASED, POS, REF, ALT, MAF_EST)]
  }, 
  split_list_job_id_19, more.args = list(reg.result = tmp_results_19), 
  reg = reg_reduce_19
)
ids_19[, chunk := 1]


#run batchmap jobs for chromosome 19
submitJobs(
  ids = ids_19,
  reg = reg_reduce_19,
  resources = list(
    ntasks = 1, ncpus = 1, memory = 250000, 
    walltime = 0L,
    partition = "prio",
    account = "dzhkomics",
    chunks.as.arrayjobs = TRUE
  )
)

#convert results of batchmap fpr chromosome 22 to data.table as an R object
reduced_tables_19 <- rbindlist(reduceResultsList(reg = reg_reduce_19), idcol = "job.id")

#add SNP specific info from original data set
variant_info_19 <- data.table::fread(dir_input_data_true_annotated_variant_infos_19, na.strings = ".")
setnames(variant_info_19, c("CHROM", "POS", "REF", "ALT", "VQSLOD", "QD", "InbreedingCoeff", "DP", "MAF"))
variant_info_19[, MAF := 0.5 - abs(MAF - 0.5)]
reduced_tables_19 <- reduced_tables_19[variant_info_19, on = .(POS, REF, ALT), nomatch = 0L]

#save data
#saveRDS(reduced_tables_19, file = dir_out_table_19)



ids_pbwt_19 <- batchMap(
  fun = function(ids, reg.result) {
    res <- rbindlist(reduceResultsList(ids = ids, reg = reg.result))
    res[, DOSAGE := GP_1 + 2*GP_2]
    res[, SUPERDOSAGE := GP_1 + 2*GP_2]
    n <- res[, length(unique(SAMPLE))]
    res[, MAF_EST := sum(DOSAGE)/(2*n), by = .( POS, REF, ALT)]
    res[, .(
      MIN_H_SCORE = min(H_SCORE, na.rm = TRUE),
      Q25_H_SCORE = quantile(H_SCORE, probs = 0.25, na.rm = TRUE),
      MEDIAN_H_SCORE = median(H_SCORE, na.rm = TRUE),
      MEAN_H_SCORE = mean(H_SCORE, na.rm = TRUE),
      SD_H_SCORE = sd(H_SCORE, na.rm = TRUE),
      Q75_H_SCORE = quantile(H_SCORE, probs = 0.75, na.rm = TRUE), 
      MAX_H_SCORE = max(H_SCORE, na.rm = TRUE),
      MIN_SEN_SCORE = min(SEN_SCORE, na.rm = TRUE),
      Q25_SEN_SCORE = quantile(SEN_SCORE, probs = 0.25, na.rm = TRUE),
      MEDIAN_SEN_SCORE = median(SEN_SCORE, na.rm = TRUE),
      MEAN_SEN_SCORE = mean(SEN_SCORE, na.rm = TRUE),
      SD_SEN_SCORE = sd(SEN_SCORE, na.rm = TRUE),
      Q75_SEN_SCORE = quantile(SEN_SCORE, probs = 0.75, na.rm = TRUE), 
      MAX_SEN_SCORE = max(SEN_SCORE, na.rm = TRUE),
      IQS = iqs(.SD[!is.na(GT)]),
      MACHR2 = (sum(DOSAGE^2)/n-(sum(DOSAGE)/n)^2)/(2*MAF_EST*(1-MAF_EST)),
      BEAGLER2 = (sum(BEST_GUESS*DOSAGE)-1/n*sum(BEST_GUESS)*sum(DOSAGE))^2/((sum(SUPERDOSAGE)-1/n*sum(DOSAGE)^2)*(sum(BEST_GUESS^2)-1/n*sum(BEST_GUESS)^2)),
      IMPUTEINFO = 1-sum(SUPERDOSAGE-DOSAGE^2)/(2*n*MAF_EST*(1-MAF_EST))
    ),
    by = .(IMPUTED, PHASED, POS, REF, ALT, MAF_EST)]
  }, 
  split_list_job_id_pbwt_19, more.args = list(reg.result = tmp_results_pbwt_19), 
  reg = reg_reduce_pbwt_19
)
ids_pbwt_19[, chunk := 1]


#run batchmap jobs for chromosome 19
submitJobs(
  ids = ids_pbwt_19,
  reg = reg_reduce_pbwt_19,
  resources = list(
    ntasks = 1, ncpus = 1, memory = 250000, 
    walltime = 0L,
    partition = "prio",
    account = "dzhkomics",
    chunks.as.arrayjobs = TRUE
  )
)

reduced_tables_pbwt_19 <- rbindlist(reduceResultsList(reg = reg_reduce_19), idcol = "job.id")

#add SNP specific info from original data set
reduced_tables_pbwt_19 <- reduced_tables_pbwt_19[variant_info_19, on = .(POS, REF, ALT), nomatch = 0L]





#get another table for concordance information



reduce_conc_22 <- makeRegistry(
  file.dir = file.path(dir_registry, "stahl_reduce_conc_22"),
  make.default = FALSE,
  source = file.path("source_scripts","result_source.R")
)


reduce_conc_19 <- makeRegistry(
  file.dir = file.path(dir_registry, "stahl_reduce_conc_19"),
  make.default = FALSE,
  source = file.path("source_scripts","result_source.R")
)


reduce_conc_pbwt_19 <- makeRegistry(
  file.dir = file.path(dir_registry, "stahl_reduce_conc_pbwt_19"),
  make.default = FALSE,
  source = file.path("source_scripts","result_source.R")
)







ids_conc_22 <- batchMap(
  fun = function(ids, reg.result) {
    res <- rbindlist(reduceResultsList(ids = ids, reg = reg.result))
    res[, .(
      CONCORDANCE_95 = concordance(.SD[!is.na(GT) & MAX_GP>=0.95]),
      CONCORDANCE_90 = concordance(.SD[!is.na(GT) & MAX_GP>=0.90 & MAX_GP<0.95]),
      CONCORDANCE_85 = concordance(.SD[!is.na(GT) & MAX_GP>=0.85 & MAX_GP<0.90]),
      CONCORDANCE_80 = concordance(.SD[!is.na(GT) & MAX_GP>=0.80 & MAX_GP<0.85]),
      CONCORDANCE_75 = concordance(.SD[!is.na(GT) & MAX_GP>=0.75 & MAX_GP<0.80]),
      CONCORDANCE_70 = concordance(.SD[!is.na(GT) & MAX_GP>=0.70 & MAX_GP<0.75]),
      CONCORDANCE_65 = concordance(.SD[!is.na(GT) & MAX_GP>=0.65 & MAX_GP<0.70]),
      CONCORDANCE_60 = concordance(.SD[!is.na(GT) & MAX_GP>=0.60 & MAX_GP<0.65]),
      CONCORDANCE_55 = concordance(.SD[!is.na(GT) & MAX_GP>=0.55 & MAX_GP<0.60]),
      CONCORDANCE_50 = concordance(.SD[!is.na(GT) & MAX_GP>=0.50 & MAX_GP<0.55]),
      C_F_99 = concordance(.SD[!is.na(GT) & MAX_GP>=0.99]),
      C_F_98 = concordance(.SD[!is.na(GT) & MAX_GP>=0.98 & MAX_GP<0.99]),
      C_F_97 = concordance(.SD[!is.na(GT) & MAX_GP>=0.97 & MAX_GP<0.98]),
      C_F_96 = concordance(.SD[!is.na(GT) & MAX_GP>=0.96 & MAX_GP<0.97]),
      C_F_95 = concordance(.SD[!is.na(GT) & MAX_GP>=0.95 & MAX_GP<0.96]),
      C_F_94 = concordance(.SD[!is.na(GT) & MAX_GP>=0.94 & MAX_GP<0.95]),
      C_F_93 = concordance(.SD[!is.na(GT) & MAX_GP>=0.93 & MAX_GP<0.94]),
      C_F_92 = concordance(.SD[!is.na(GT) & MAX_GP>=0.92 & MAX_GP<0.93]),
      C_F_91 = concordance(.SD[!is.na(GT) & MAX_GP>=0.91 & MAX_GP<0.92]),
      C_F_90 = concordance(.SD[!is.na(GT) & MAX_GP>=0.90 & MAX_GP<0.91]),
      CONCORDANCE_MA = concordance(.SD[!is.na(GT) & GT!=0]),
      CONCORDANCE = concordance(.SD[!is.na(GT)])),
      by = .(POS,REF,ALT, IMPUTED, PHASED)]
  }, 
  split_list_job_id_22, more.args = list(reg.result = tmp_results_22), 
  reg = reduce_conc_22
)


ids_conc_22[, chunk := 1]

submitJobs(
  ids = ids_conc_22,
  reg = reduce_conc_22,
  resources = list(
    ntasks = 1, ncpus = 1, memory = 250000, 
    walltime = 0L,
    partition = "prio",
    account = "dzhkomics",
    chunks.as.arrayjobs = TRUE
  )
)


ids_conc_19 <- batchMap(
  fun = function(ids, reg.result) {
    res <- rbindlist(reduceResultsList(ids = ids, reg = reg.result))
    res[, .(
      CONCORDANCE_95 = concordance(.SD[!is.na(GT) & MAX_GP>=0.95]),
      CONCORDANCE_90 = concordance(.SD[!is.na(GT) & MAX_GP>=0.90 & MAX_GP<0.95]),
      CONCORDANCE_85 = concordance(.SD[!is.na(GT) & MAX_GP>=0.85 & MAX_GP<0.90]),
      CONCORDANCE_80 = concordance(.SD[!is.na(GT) & MAX_GP>=0.80 & MAX_GP<0.85]),
      CONCORDANCE_75 = concordance(.SD[!is.na(GT) & MAX_GP>=0.75 & MAX_GP<0.80]),
      CONCORDANCE_70 = concordance(.SD[!is.na(GT) & MAX_GP>=0.70 & MAX_GP<0.75]),
      CONCORDANCE_65 = concordance(.SD[!is.na(GT) & MAX_GP>=0.65 & MAX_GP<0.70]),
      CONCORDANCE_60 = concordance(.SD[!is.na(GT) & MAX_GP>=0.60 & MAX_GP<0.65]),
      CONCORDANCE_55 = concordance(.SD[!is.na(GT) & MAX_GP>=0.55 & MAX_GP<0.60]),
      CONCORDANCE_50 = concordance(.SD[!is.na(GT) & MAX_GP>=0.50 & MAX_GP<0.55]),
      C_F_99 = concordance(.SD[!is.na(GT) & MAX_GP>=0.99]),
      C_F_98 = concordance(.SD[!is.na(GT) & MAX_GP>=0.98 & MAX_GP<0.99]),
      C_F_97 = concordance(.SD[!is.na(GT) & MAX_GP>=0.97 & MAX_GP<0.98]),
      C_F_96 = concordance(.SD[!is.na(GT) & MAX_GP>=0.96 & MAX_GP<0.97]),
      C_F_95 = concordance(.SD[!is.na(GT) & MAX_GP>=0.95 & MAX_GP<0.96]),
      C_F_94 = concordance(.SD[!is.na(GT) & MAX_GP>=0.94 & MAX_GP<0.95]),
      C_F_93 = concordance(.SD[!is.na(GT) & MAX_GP>=0.93 & MAX_GP<0.94]),
      C_F_92 = concordance(.SD[!is.na(GT) & MAX_GP>=0.92 & MAX_GP<0.93]),
      C_F_91 = concordance(.SD[!is.na(GT) & MAX_GP>=0.91 & MAX_GP<0.92]),
      C_F_90 = concordance(.SD[!is.na(GT) & MAX_GP>=0.90 & MAX_GP<0.91]),
      CONCORDANCE_MA = concordance(.SD[!is.na(GT) & GT!=0]),
      CONCORDANCE = concordance(.SD[!is.na(GT)])),
      by = .(POS, REF, ALT, IMPUTED, PHASED)]
  }, 
  split_list_job_id_19, more.args = list(reg.result = tmp_results_19), 
  reg = reduce_conc_19
)


ids_conc_19[, chunk := 1]

submitJobs(
  ids = ids_conc_19,
  reg = reduce_conc_19,
  resources = list(
    ntasks = 1, ncpus = 1, memory = 250000, 
    walltime = 0L,
    partition = "prio",
    account = "dzhkomics",
    chunks.as.arrayjobs = TRUE
  )
)




ids_conc_pbwt_19 <- batchMap(
  fun = function(ids, reg.result) {
    res <- rbindlist(reduceResultsList(ids = ids, reg = reg.result))
    res[, .(
      CONCORDANCE_95 = concordance(.SD[!is.na(GT) & MAX_GP>=0.95]),
      CONCORDANCE_90 = concordance(.SD[!is.na(GT) & MAX_GP>=0.90 & MAX_GP<0.95]),
      CONCORDANCE_85 = concordance(.SD[!is.na(GT) & MAX_GP>=0.85 & MAX_GP<0.90]),
      CONCORDANCE_80 = concordance(.SD[!is.na(GT) & MAX_GP>=0.80 & MAX_GP<0.85]),
      CONCORDANCE_75 = concordance(.SD[!is.na(GT) & MAX_GP>=0.75 & MAX_GP<0.80]),
      CONCORDANCE_70 = concordance(.SD[!is.na(GT) & MAX_GP>=0.70 & MAX_GP<0.75]),
      CONCORDANCE_65 = concordance(.SD[!is.na(GT) & MAX_GP>=0.65 & MAX_GP<0.70]),
      CONCORDANCE_60 = concordance(.SD[!is.na(GT) & MAX_GP>=0.60 & MAX_GP<0.65]),
      CONCORDANCE_55 = concordance(.SD[!is.na(GT) & MAX_GP>=0.55 & MAX_GP<0.60]),
      CONCORDANCE_50 = concordance(.SD[!is.na(GT) & MAX_GP>=0.50 & MAX_GP<0.55]),
      C_F_99 = concordance(.SD[!is.na(GT) & MAX_GP>=0.99]),
      C_F_98 = concordance(.SD[!is.na(GT) & MAX_GP>=0.98 & MAX_GP<0.99]),
      C_F_97 = concordance(.SD[!is.na(GT) & MAX_GP>=0.97 & MAX_GP<0.98]),
      C_F_96 = concordance(.SD[!is.na(GT) & MAX_GP>=0.96 & MAX_GP<0.97]),
      C_F_95 = concordance(.SD[!is.na(GT) & MAX_GP>=0.95 & MAX_GP<0.96]),
      C_F_94 = concordance(.SD[!is.na(GT) & MAX_GP>=0.94 & MAX_GP<0.95]),
      C_F_93 = concordance(.SD[!is.na(GT) & MAX_GP>=0.93 & MAX_GP<0.94]),
      C_F_92 = concordance(.SD[!is.na(GT) & MAX_GP>=0.92 & MAX_GP<0.93]),
      C_F_91 = concordance(.SD[!is.na(GT) & MAX_GP>=0.91 & MAX_GP<0.92]),
      C_F_90 = concordance(.SD[!is.na(GT) & MAX_GP>=0.90 & MAX_GP<0.91]),
      CONCORDANCE_MA = concordance(.SD[!is.na(GT) & GT!=0]),
      CONCORDANCE = concordance(.SD[!is.na(GT)])),
      by = .(POS, REF, ALT, IMPUTED, PHASED)]
  }, 
  split_list_job_id_pbwt_19, more.args = list(reg.result = tmp_results_pbwt_19), 
  reg = reduce_conc_pbwt_19
)

ids_conc_pbwt_19[, chunk := 1]

submitJobs(
  ids = ids_conc_pbwt_19,
  reg = reduce_conc_pbwt_19,
  resources = list(
    ntasks = 1, ncpus = 1, memory = 250000, 
    walltime = 0L,
    partition = "prio",
    account = "dzhkomics",
    chunks.as.arrayjobs = TRUE
  )
)



reduced_tables_conc_22 <- rbindlist(reduceResultsList(reg = reduce_conc_22), idcol = "job.id")
reduced_tables_conc_19 <- rbindlist(reduceResultsList(reg = reduce_conc_19), idcol = "job.id")
reduced_tables_conc_pbwt_19 <- rbindlist(reduceResultsList(reg = reduce_conc_pbwt_19), idcol = "job.id")


#save data
#saveRDS(reduced_tables_conc_22, file = dir_out_table_conc_22)
#saveRDS(reduced_tables_conc_19, file = dir_out_table_conc_19)
#saveRDS(reduced_tables_conc_pbwt_19, file = dir_out_table_conc_pbwt_19)