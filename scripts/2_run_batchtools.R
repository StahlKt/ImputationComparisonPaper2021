# run batchtools
#this is the second file to execute to run batchtools
#Please refer to the batchtools package for details on the structure of the specific functions




source_dir<-file.path(home_dir,"source_scripts")
source(file.path(source_dir,"directories.R"))
# checks output folders, if creates them if they don't exist
sapply(mget(ls(pattern = "out")), check_and_create_folders)

library("batchtools")


#make experiment registry for phasing and Imputation
reg <- makeExperimentRegistry(
    file.dir = file.path(dir_registry, "stahl_phasing_imputation"), 
    make.default = TRUE, 
    source = file.path(source_dir,"directories.R"), 
    conf.file = confic_dir
)


#define phasing as "problem", so all phasing results will be imputed by all imputation programs
#chromosome 22
addProblem("phase_shapeit_22", data = dir_input_data_annotated_22, fun = shapeit_wrap, seed = 20190509, cache = TRUE, reg = reg)
addProblem("phase_eagle_22", data = dir_input_data_annotated_22, fun = eagle_wrap, seed = 20190509, cache = TRUE, reg = reg)
addProblem("phase_beagle5_22", data = dir_input_data_annotated_22, fun = beagle5_phasing_wrap, seed = 20190509, cache = TRUE, reg = reg)


#define arguments for phasing functions for chromosome 22
arguments_phasing_22 <- list(
    phase_shapeit_22 = data.table(
        exec = dir_exec_shapeit, 
        outpref = file.path(dir_out_shapeit, "phased_shapeit_22"),
        map.no.chr = dir_input_map_ohne_chr_22, 
        reference.hap = dir_input_reference_hap_22, 
        ref.on = c(1, 0),
        chr.flag = 22
    ),
    phase_eagle_22 = data.table(
        exec = dir_exec_eagle, 
        outpref = file.path(dir_out_eagle, "phased_eagle_22"), 
        map.with.chr = dir_input_map_mit_chr_22, 
        reference.vcf = dir_input_reference_vcf_annotated_22, 
        ref.on = c(1, 0),
        chr.flag = 22
    ), 
    phase_beagle5_22 = data.table(
        exec = dir_exec_beagle5,
        outpref = file.path(dir_out_beagle5_phasing, "phased_beagle5_22"), 
        plink.map = dir_input_map_plink_22,
        chr.flag = 22)
)

#add phasing of chromosome 19 as separate problems
addProblem("phase_shapeit_19", data = dir_input_data_annotated_19, fun = shapeit_wrap, seed = 20190509, cache = TRUE, reg = reg)
addProblem("phase_eagle_19", data = dir_input_data_annotated_19, fun = eagle_wrap, seed = 20190509, cache = TRUE, reg = reg)
addProblem("phase_beagle5_19", data = dir_input_data_annotated_19, fun = beagle5_phasing_wrap, seed = 20190509, cache = TRUE, reg = reg)

#define arguments for phasing functions for chromosome 19
arguments_phasing_19 <- list(
    phase_shapeit_19 = data.table(
        exec = dir_exec_shapeit, 
        outpref = file.path(dir_out_shapeit, "phased_shapeit_19"),
        map.no.chr = dir_input_map_ohne_chr_19, 
        reference.hap = dir_input_reference_hap_19, 
        ref.on = c(1, 0),
        chr.flag=19
    ),
    phase_eagle_19 = data.table(
        exec = dir_exec_eagle, 
        outpref = file.path(dir_out_eagle, "phased_eagle_19"), 
        map.with.chr = dir_input_map_mit_chr_19, 
        reference.vcf = dir_input_reference_vcf_annotated_19, 
        ref.on = c(1, 0),
        chr.flag = 19
    ), 
    phase_beagle5_19 = data.table(
        exec = dir_exec_beagle5,
        outpref = file.path(dir_out_beagle5_phasing, "phased_beagle5_19"), 
        plink.map = dir_input_map_plink_19,
        chr.flag = 19)
)


#define imputation as "algorithms" for batchtools, in this step no distinction between chromosomes is needed
#       because it is included in the arguments and the right data set will be inherited from the "problem"
addAlgorithm("impute_impute2", fun = impute2_wrap, reg = reg)
addAlgorithm("impute_impute4", fun = impute4_wrap, reg = reg)
addAlgorithm("impute_beagle4", fun = beagle4_wrap, reg = reg)
addAlgorithm("impute_beagle5", fun = beagle5_impute_wrap, reg = reg)
addAlgorithm("impute_minimac3", fun = minimac3_wrap, reg = reg)
addAlgorithm("impute_minimac4", fun = minimac4_wrap, reg = reg)
addAlgorithm("impute_pbwt", fun = pbwt_wrap, reg = reg)

#define arguments for imputation functions for chromosome 22
arguments_imputation_22 <- list(
    impute_impute2 = data.table(
        exec = dir_exec_impute2, 
        outpref = file.path(dir_out_impute2, "imputed_impute2"),
        map.no.chr = dir_input_map_ohne_chr_22, 
        reference.hap = dir_input_reference_hap_22,
        int = list_of_chunks_22,
        chr.flag = 22
    ), 
    impute_impute4 = data.table(
        exec = dir_exec_impute4, 
        outpref = file.path(dir_out_impute4, "imputed_impute4"), 
        map.no.chr = dir_input_map_ohne_chr_22,
        reference.hap = dir_input_reference_hap_22,
        int = list_of_chunks_22,
        chr.flag = 22
    ), 
    impute_beagle4 = data.table(
        exec = dir_exec_beagle4, 
        outpref = file.path(dir_out_beagle4, "imputed_beagle4"),
        plink.map = dir_input_map_plink_22,
        reference.vcf = dir_input_reference_vcf_annotated_22,
        int = list_of_chunks_22,
        chr.flag = 22
    ), 
    impute_beagle5 = data.table(
        exec = dir_exec_beagle5, 
        outpref = file.path(dir_out_beagle5_imputation, "imputed_beagle5"),
        plink.map = dir_input_map_plink_22, 
        reference.bref3 = dir_input_reference_bref3_22, 
        int = list_of_chunks_22,
        chr.flag = 22
    ), 
    impute_minimac3 = data.table(
        exec = dir_exec_minimac3,
        outpref = file.path(dir_out_minimac3, "imputed_minimac3"), 
        reference.vcf = dir_input_reference_vcf_annotated_22,
        int = list_of_chunks_22,
        chr.flag = 22
    ), 
    impute_minimac4 = data.table(
        exec = dir_exec_minimac4, 
        outpref = file.path(dir_out_minimac4, "imputed_minimac4"), 
        map.with.chr = dir_input_map_mit_chr_22, 
        reference.m3vcf = dir_input_reference_m3vcf_22,
        int = list_of_chunks_22,
        chr.flag = 22
    ),
    impute_pbwt = data.table(
      exec = dir_exec_pbwt, 
      outpref = file.path(dir_out_pbwt, "imputed_pbwt"),
      reference.pbwt = dir_input_reference_pbwt_22,
      chr.flag = 22
    )
)



#define arguments for imputation functions for chromosome 19,
#           PBWT excluded, because it needs a different chunking, see run_batchtools_pbwt_19
arguments_imputation_19 <- list(
    impute_impute2 = data.table(
        exec = dir_exec_impute2, 
        outpref = file.path(dir_out_impute2, "imputed_impute2"),
        map.no.chr = dir_input_map_ohne_chr_19, 
        reference.hap = dir_input_reference_hap_19,
        int = list_of_chunks_19,
        chr.flag = 19
    ), 
    impute_impute4 = data.table(
        exec = dir_exec_impute4, 
        outpref = file.path(dir_out_impute4, "imputed_impute4"), 
        map.no.chr = dir_input_map_ohne_chr_19,
        reference.hap = dir_input_reference_hap_19,
        int = list_of_chunks_19,
        chr.flag = 19
    ), 
    impute_beagle4 = data.table(
        exec = dir_exec_beagle4, 
        outpref = file.path(dir_out_beagle4, "imputed_beagle4"),
        plink.map = dir_input_map_plink_19,
        reference.vcf = dir_input_reference_vcf_annotated_19,
        int = list_of_chunks_19,
        chr.flag = 19
    ), 
    impute_beagle5 = data.table(
        exec = dir_exec_beagle5, 
        outpref = file.path(dir_out_beagle5_imputation, "imputed_beagle5"),
        plink.map = dir_input_map_plink_19, 
        reference.bref3 = dir_input_reference_bref3_19, 
        int = list_of_chunks_19,
        chr.flag = 19
    ), 
    impute_minimac3 = data.table(
        exec = dir_exec_minimac3,
        outpref = file.path(dir_out_minimac3, "imputed_minimac3"), 
        reference.vcf = dir_input_reference_vcf_annotated_19,
        int = list_of_chunks_19,
        chr.flag = 19
    ), 
    impute_minimac4 = data.table(
        exec = dir_exec_minimac4, 
        outpref = file.path(dir_out_minimac4, "imputed_minimac4"), 
        map.with.chr = dir_input_map_mit_chr_19, 
        reference.m3vcf = dir_input_reference_m3vcf_19,
        int = list_of_chunks_19,
        chr.flag = 19
    ),
    impute_pbwt = data.table(
      exec = dir_exec_pbwt, 
      outpref = file.path(dir_out_pbwt, "imputed_pbwt"),
      reference.pbwt = dir_input_reference_pbwt_22,
      chr.flag = 22
    )
)



#add phasing/imputation pipelines as "experiments" for both chromosomes
addExperiments(prob.designs = arguments_phasing_22, algo.designs = arguments_imputation_22, reg = reg)
addExperiments(prob.designs = arguments_phasing_19, algo.designs = arguments_imputation_19, reg = reg)

#check, if all experiments are added and have the right combinations of phasing and imputation
summarizeExperiments(reg = reg)


#first run only one set of each phasing, so the phasing step is only done once and not everytime anew
#         for each imputation chunk

#extract job table
job_table <- getJobTable(reg = reg)
job_table <- unwrap(job_table[, .(job.id, problem, prob.pars)])

#extract of each phasing one job to run, use beagle5.1 as imputation for speed reasons 
first_run_ids <- job_table[, .(job.id = .SD[1, job.id]), by = .(problem, ref.on, algorithm)]
first_run_ids <- job_table[algorithm=="impute_beagle5", .(job.id = .SD[1, job.id]), by = .(problem, ref.on)]
#set chunk for batchtools on one, because chunking is already implicit in the algorithm design
first_run_ids[, chunk := 1]

#run first 5 jobs, one for each phasing
submitJobs(
    ids = first_run_ids, 
    reg = reg,
    resources = list(
        ntasks = 1, ncpus = 1, memory = 80000, 
        walltime = 0L,
        partition = "batch",
        chunks.as.arrayjobs = TRUE
    )
)

waitForJobs(reg = reg)




#find jobs that are not done to submit individually because of different memory needs
impute2_ids <- findExperiments(algo.name = "impute_impute2")
impute2_ids[, chunk := 1]

impute4_ids <- findExperiments(algo.name = "impute_impute4")
impute4_ids[, chunk := 1]


beagle4_ids <- findExperiments(algo.name = "impute_beagle4")
beagle4_ids[, chunk := 1]

beagle5_ids <- findExperiments(algo.name = "impute_beagle5")
beagle5_ids <- beagle5_ids[(-beagle5_ids[first_run_ids, on = "job.id", which = TRUE, nomatch = 0])]
beagle5_ids[, chunk := 1]

minimac3_ids <- findExperiments(algo.name = "impute_minimac3")
minimac3_ids[, chunk := 1]


minimac4_ids <- findExperiments(algo.name = "impute_minimac4")
minimac4_ids[, chunk := 1]


pbwt_ids <- findExperiments(algo.name = "impute_pbwt")
pbwt_ids[, chunk := 1]





#run PBWT imputation
submitJobs(
  ids = pbwt_ids,
  reg = reg,
  resources = list(
    ntasks = 1, ncpus = 1, memory = 250000, 
    walltime = 0L,
    partition = "prio",
    account= "dzhkomics",
    chunks.as.arrayjobs = TRUE
  )
)



#run IMPUTE2 imputation
submitJobs(
    ids = impute2_ids,
    reg = reg,
    resources = list(
        ntasks = 1, ncpus = 1, memory = 80000, 
        walltime = 0L,
        partition = "batch",
        chunks.as.arrayjobs = TRUE
    )
)

#run IMPUTE4 imputation
submitJobs(
    ids = impute4_ids,
    reg = reg,
    resources = list(
        ntasks = 1, ncpus = 1, memory = 20000, 
        walltime = 0L,
        partition = "batch",
        chunks.as.arrayjobs = TRUE
    )
)





#run Beagle4.1 imputation
submitJobs(
    ids = beagle4_ids,
    reg = reg,
    resources = list(
        ntasks = 1, ncpus = 1, memory = 20000, 
        walltime = 0L,
        partition = "batch",
        chunks.as.arrayjobs = TRUE
    )
)

#run Beagle5.1 imputation
submitJobs(
    ids = beagle5_ids,
    reg = reg,
    resources = list(
        ntasks = 1, ncpus = 1, memory = 40000, 
        walltime = 0L,
        partition = "batch",
        chunks.as.arrayjobs = TRUE
    )
)




#run minimac3 imputation
submitJobs(
    ids = minimac3_ids,
    reg = reg,
    resources = list(
        ntasks = 1, ncpus = 1, memory = 40000, 
        walltime = 0L,
        partition = "batch",
        chunks.as.arrayjobs = TRUE
    )
)


#run minimac4 imputation
submitJobs(
    ids = minimac4_ids,
    reg = reg,
    resources = list(
        ntasks = 1, ncpus = 1, memory = 40000, 
        walltime = 0L,
        partition = "batch",
        chunks.as.arrayjobs = TRUEnot
    )
)

waitForJobs(reg = reg)







