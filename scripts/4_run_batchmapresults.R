#fourth to run batchtools for imputation quality comparison


source(file.path(source_dir,"result_source.R"))




#split results of jobs from batchtools into 2 groups according to chromosomes 
reg<-loadRegistry(file.path(dir_registry,"stahl_phasing_imputation"), writeable=TRUE)


job_table <- getJobTable(reg = reg)
job_table <- job_table[, .(job.id, problem, algorithm)]

jobs_19<-job_table[problem=="phase_shapeit_19" |
                     problem=="phase_beagle5_19"|
                     problem=="phase_eagle_19"]

jobs_22<-job_table[problem=="phase_shapeit_22" |
                     problem=="phase_beagle5_22"|
                     problem=="phase_eagle_22"]




#make new registry for further operations on the batchtools results for chromosome 19
tmp_results_19 <- makeRegistry(
  file.dir = file.path(dir_registry, "stahl_batchmapresults_19"), 
  make.default = FALSE,
  source = c(file.path(source_dir, "directories.R"), 
             file.path(source_dir,"result_source.R")),
  conf.file = confic_dir
)

#define batchmaptresults call with the function "build_table_results2"  for chromosome 19, 
#       which adds the true genotype to the imputed genotypes and calculates first imputation quality scores
ids_19 <- batchMapResults(
  build_table_results2, 
  ids = merge(findDone(reg = reg),jobs_19), 
  target = tmp_results_19, 
  source = reg,
  more.args = list(
    true.data.annotated = dir_input_data_true_annotated_sample_genotypes_19
  )
)
ids_19[, chunk := 1]


#run batchmapresults to get a table with both imputed and true genotypes for chromosome 19
submitJobs(
  ids = ids_19,
  reg = tmp_results_19,
  resources = list(
    ntasks = 1, ncpus = 1, memory = 99000, 
    walltime = 0L,
    partition = "batch",
    chunks.as.arrayjobs = TRUE
  )
)







#make new registry for further operations on the batchtools results for chromosome 22
tmp_results_22 <- makeRegistry(
  file.dir = file.path(dir_registry, "stahl_batchmapresults_22"), 
  make.default = FALSE,
  source = c(file.path(source_dir, "directories.R"), 
             file.path(source_dir,"result_source.R")),
  conf.file = confic_dir
)

#define batchmaptresults call with the function "build_table_results2"  for chromosome 22, 
#       which adds the true genotype to the imputed genotypes and calculates first imputation quality scores
ids_22 <- batchMapResults(
  build_table_results2, 
  ids = merge(findDone(reg = reg),jobs_22), 
  target = tmp_results_22, 
  source = reg,
  more.args = list(
    true.data.annotated = dir_input_data_true_annotated_sample_genotypes_22
  )
)
ids_22[, chunk := 1]

#run batchmapresults to get a table with both imputed and true genotypes for chromosome 22
submitJobs(
  ids = ids_22,
  reg = tmp_results_22,
  resources = list(
    ntasks = 1, ncpus = 1, memory = 99000, 
    walltime = 0L,
    partition = "batch",
    chunks.as.arrayjobs = TRUE
  )
)

#start jobs again with higher memory, if they run out





#batchmapresults for PBWT imputation on chromosome 19
reg_pbwt_19<-loadRegistry(file.path(dir_registry,"stahl_imputation_pbwt_19"),make.default=FALSE, writeable=TRUE)


#make new registry for further operations on the batchtools results for pbwt imputation on chromosome 19
tmp_results_pbwt_19 <- makeRegistry(
  file.dir = file.path(dir_registry, "stahl_bmresults_pbwt_19"), 
  make.default = FALSE,
  source = c(file.path(source_dir, "directories.R"), 
             file.path(source_dir, "result_source.R")),
  conf.file = confic_dir
)



ids_19_pbwt <- batchMapResults(
  build_table_results2, 
  ids = findDone(reg = reg_pbwt_19), 
  target = tmp_results_pbwt_19, 
  source = reg_pbwt_19,
  more.args = list(
    true.data.annotated = dir_input_data_true_annotated_sample_genotypes_19
  )
)
ids_19_pbwt[, chunk := 1]

submitJobs(
  ids = ids_19_pbwt,
  reg = tmp_results_pbwt_19,
  resources = list(
    ntasks = 1, ncpus = 1, memory = 250000, 
    walltime = 0L,
    partition = "prio",
    account="dzhkomics",
    chunks.as.arrayjobs = TRUE
  )
)
