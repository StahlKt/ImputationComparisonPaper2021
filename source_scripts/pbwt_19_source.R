#functions and directories for pbwt_19


home_dir<-"/home/stahl/MA"
setwd(home_dir)
library(batchtools)
source("source_scripts/directories.R")


dir_phasing_results<-file.path(dir_output,"phasing")

dir_input_individuals_first_half<-file.path(dir_phasing_results,"samples_19_first.txt")

dir_input_individuals_second_half<-file.path(dir_phasing_results,"samples_19_last.txt")







divide_by_individual<- function(data, 
                                job,
                                exec,
                                which.phasing,
                                file.individuals){
  
  if(grepl(file.individuals, pattern="first")){
    flag_individuals<-"first_half"
  }
  else if(grepl(file.individuals, pattern="last")){
    flag_individuals<-"second_half"
  }



full_name<-paste(which.phasing,"_pbwt.vcf.gz",sep="")
input<-file.path(data,full_name)
input_no_suffix<-file.path(data,which.phasing)
output<-paste(input_no_suffix,"_",flag_individuals,"_pbwt.vcf.gz",sep="")
out_no_suffix<-paste(input_no_suffix,"_",flag_individuals,sep="")


system2(exec,c("view", 
               "-S", file.individuals,
               "-O z",
               "-o", output,
               input))
        
        
        
return<-c(out_no_suffix,paste(get_phasing_program(which.phasing),"_",flag_individuals, sep=""))
}




get_phasing_program<-function(datadir){
  
  if(grepl(datadir, pattern="beagle5_19")){
    phasing<-"phased_beagle5"
  }
  else if(grepl(datadir, pattern="shapeit_19_mit_ref")){
    phasing<-"phased_shapeit_ref_on"
  }
  else if(grepl(datadir, pattern="shapeit_19_ohne_ref")){
    phasing<-"phased_shapeit_ref_off"
  }
  else if(grepl(datadir, pattern="eagle_19_ohne_ref")){
    phasing<-"phased_eagle_ref_off"
  }
  else if(grepl(datadir, pattern="eagle_19_mit_ref")){
      phasing<-"phased_eagle_ref_on"
  }
    
  return(phasing)
}
  
  
  
  
