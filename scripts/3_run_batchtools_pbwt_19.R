#imputation with PBWT
#this is the third file to execute to run batchtools


#set up new registry for pbwt
reg_pbwt_19<- makeExperimentRegistry(
  file.dir = file.path(dir_registry, "stahl_imputation_pbwt_19"),
  make.default = FALSE,
  source = file.path(source_dir,"pbwt_19_source.R"),
  conf.file = confic_dir
)

#in case you need to load the already defined registry at a later point
#reg_pbwt_19<-loadRegistry(file.path(dir_registry,"stahl_imputation_pbwt_19"),make.default=FALSE, writeable=TRUE)

#add a "Problem", which divides the already phased data into two sets of participants
addProblem("divide_phased_data", data=dir_phasing_results, fun=divide_by_individual, seed = 20190509, cache = TRUE, reg = reg_pbwt_19 )


#define arguments for division function for the first half of individuals
args_problem_pbwt_19_first<-list(
  divide_phased_data=data.table(
    exec=dir_exec_bcftools,
    which.phasing=c("Eagle2/phased_eagle_19_mit_ref",
                    "Eagle2/phased_eagle_19_ohne_ref",
                    "Shapeit2/phased_shapeit_19_mit_ref",
                    "Shapeit2/phased_shapeit_19_ohne_ref",
                    "Beagle5/phased_beagle5_19"),
  file.individuals=dir_input_individuals_first_half
                     
  )
)  


#define arguments for division function for the second half of individuals
args_problem_pbwt_19_last<-list(
  divide_phased_data=data.table(
    exec=dir_exec_bcftools,
    which.phasing=c("Eagle2/phased_eagle_19_mit_ref",
                    "Eagle2/phased_eagle_19_ohne_ref",
                    "Shapeit2/phased_shapeit_19_mit_ref",
                    "Shapeit2/phased_shapeit_19_ohne_ref",
                    "Beagle5/phased_beagle5_19"),
    file.individuals=dir_input_individuals_second_half
  )
)  


#add imputation with pbwt on chromosome 19
addAlgorithm("impute_pbwt_19", fun = pbwt_wrap, reg = reg_pbwt_19)

args_algo_pbwt_19<-list(
    impute_pbwt_19 = data.table(
    exec = dir_exec_pbwt, 
    outpref = file.path(dir_out_pbwt, "imputed_pbwt"),
    reference.pbwt = dir_input_reference_pbwt_19,
    chr.flag = 19
    )
)



#define division of datasets imputation with PBWT as an experiment for the according halfs
addExperiments(prob.designs = args_problem_pbwt_19_first, algo.designs = args_algo_pbwt_19, reg = reg_pbwt_19)
addExperiments(prob.designs = args_problem_pbwt_19_last, algo.designs = args_algo_pbwt_19, reg = reg_pbwt_19)

#check, if the problems and algorithms are matched correctly 
summarizeExperiments(reg = reg_pbwt_19)

#find job ids to run and ensure that each job is processed as only one chunk
ids<-findNotDone(reg=reg_pbwt_19)
ids[,chunk:=1]


#run imputation with PBWT on chromosome 19
submitJobs(
  ids = ids, 
  reg = reg_pbwt_19,
  resources = list(
    ntasks = 1, ncpus = 1, memory = 250000, 
    walltime = 0L,
    partition = "prio",
    account = "dzhkomics",
    chunks.as.arrayjobs = TRUE
  )
)



