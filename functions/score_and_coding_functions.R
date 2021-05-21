#functions to calculate different imputation quality scores and to recode the dataset for consistency and handling


#calculate the hellinger score
hellinger_score<-function(imputed.prob,true.genotype.coded){
  score<-1-sqrt(1-sqrt(imputed.prob))
  return(score)
}


#calculate SEN-score
sen_score<-function(imputed.probs_0,imputed.probs_1,imputed.probs_2, true.genotype.coded){
  imputed.probs<-cbind(imputed.probs_0,imputed.probs_1,imputed.probs_2)
  imputed_dosage<- 2-(imputed.probs[,2]+2*imputed.probs[,1])
  score<-1-((true.genotype.coded-imputed_dosage)^2)/4
  return(score)
}

#get MAF 
get_maf<-function(genotypes){
  maf<-sum(genotypes)/(2*length(genotypes))  
  if(maf>0.5){
    maf<-1-maf
  }
  return(maf)
}

#determine, which programs were used according to the file path/file name
get_used_programs<-function(datadir){
  
  if(grepl(datadir, pattern="phased_beagle5")){
    phasing<-"beagle5"
  }
  else if(grepl(datadir, pattern="shapeit_ref_off")){
    phasing<-"shapeit no ref"
  }
  else if(grepl(datadir, pattern="shapeit_ref_on")){
    phasing<-"shapeit with ref"
  }
  else if(grepl(datadir, pattern="eagle_ref_off")){
    phasing<-"eagle no ref"
  }
  else if(grepl(datadir, pattern="eagle_ref_on")){
    phasing<-"eagle with ref"
  }
  if(grepl(datadir, pattern="imputed_beagle5")){
    imputation<-"beagle5"
  }
  else if(grepl(datadir, pattern="beagle4")){
    imputation<-"beagle4"
  }
  else if(grepl(datadir, pattern="impute2")){
    imputation<-"impute2"
  }
  else if(grepl(datadir, pattern="impute4")){
    imputation<-"impute4"
  }
  else if(grepl(datadir, pattern="minimac4")){
    imputation<-"minimac4"
  }
  else if(grepl(datadir, pattern="minimac3")){
    imputation<-"minimac3"
  }
  else if(grepl(datadir, pattern="pbwt")){
    imputation<-"pbwt"
  }
  return(c(phasing, imputation))
}



#calculate best guess genotype out of imputation probabilities
best_guess_genotype<-function(imputed.probs){
  return(best_guess<-which.max(imputed.probs)-1)  
}


#extract genotypes as 0 1 2 from vcf-type output format
code_genotypes<-function(true.genotype){
  #expect "1/1" "0/1" or "0/0" as true.genotype
  if(true.genotype=="1/1"){
    genotype_coded<-2
  }
  else if(true.genotype=="0/0"){
    genotype_coded<-0
  }  
  else{
    genotype_coded<-1
  }
  return(genotype_coded)
}


#calculate IQS
#expects data.table with one column best_guess and one true_genotype

iqs<-function(in.table){
  imputed_11<-as.numeric(in.table[BEST_GUESS==0, .N])
  imputed_12<-as.numeric(in.table[BEST_GUESS==1, .N])
  imputed_22<-as.numeric(in.table[BEST_GUESS==2, .N])
  
  true_11<-as.numeric(in.table[GT==0, .N])
  true_12<-as.numeric(in.table[GT==1, .N])
  true_22<-as.numeric(in.table[GT==2, .N])
  n<-as.numeric(length(in.table$GT))
  
  
  #observed proportion of agreement
  p_observed<-in.table[GT==BEST_GUESS, .N]/n
  
  #proportion of chance agreement
  p_chance<-(true_11*imputed_11+true_12*imputed_12+true_22*imputed_22)/n^2
  
  
  score<-(p_observed-p_chance)/(1-p_chance)
  
  return(score)
}

    
#calculate concordance rate
#expects data.table with one column best_guess and one true_genotype

concordance<-function(in.table){
  
  n<-as.numeric(length(in.table$GT))
  
  
  #observed proportion of agreement
  conc<-in.table[GT==BEST_GUESS, .N]/n
  return(conc)
}
   
    
    
    

