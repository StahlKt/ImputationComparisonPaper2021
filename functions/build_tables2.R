#essential function for batchmapresults



build_table_results2<-function(x, true.data.annotated){
  
  #read in true data for comparison
  if(file.exists(true.data.annotated)){
    true_table <- fread(input = true.data.annotated, drop=2)
    
    #set names
    setnames(true_table, c("SAMPLE","POS","REF","ALT", "GT"))
    
    #swap GT for coded GT
    true_table[, GT := as.numeric(factor(GT, levels = c("0/0", "0/1", "1/0", "1/1", "0|0", "0|1", "1|0", "1|1"), labels = c(0, 1, 1, 2, 0, 1, 1, 2)))-1]
    
    #set key to merge for later
    setkey(true_table,  POS, REF, ALT, SAMPLE)
  } else {
    stop("Annotated true data not found!")
  }
  
  #read in table with samples,chrom pos, ref, alt, GP(genotype probabilities in 3 columns)
  combined_table <- fread(
    cmd = sprintf(
      fmt = "%s %s %s",
      dir_exec_bcftools, "query -f \'[%SAMPLE,%POS,%REF,%ALT,%GP\\n]\'",
      x), 
    sep = ","
  )
  #set names for compability
  setnames(combined_table,c("SAMPLE","POS","REF","ALT","0","1", "2"))
  
  #add best_guess und max_GP
  combined_table[, BEST_GUESS := as.integer(colnames(.SD)[max.col(.SD, ties.method = "first")]), .SDcols = c("0", "1", "2")]
  combined_table[, MAX_GP := get(colnames(.SD)[max.col(.SD, ties.method = "first")]), .SDcols = c("0", "1", "2")]
  
  #set names for compability
  setnames(combined_table,c("SAMPLE","POS","REF","ALT","GP_0", "GP_1", "GP_2", "BEST_GUESS", "MAX_GP"))
  

  
  
  #get which programs were used
  col_length <- combined_table[, .N]
  programs <- get_used_programs(x)
  combined_table[, IMPUTED := programs[2]]
  combined_table[, PHASED := programs[1]]
  
  #get true genotypes and mafs from true_table
  setkey(combined_table, SAMPLE,  POS, REF, ALT)
  combined_table <- merge(combined_table, true_table, by=c("SAMPLE","POS","REF","ALT"), all = FALSE)
  
  #remove table with true genotypes to free up space
  rm(true_table)

  #calculate Hellinger score and SEN score
  combined_table[!is.na(GT), H_SCORE := hellinger_score(.SD[[sprintf("GP_%d", GT)]], GT), by = .(GT)]
  combined_table[, SEN_SCORE := sen_score(GP_0, GP_1, GP_2, GT)]
  
  #return table
  return(combined_table)
}



#get_max

get_max<-function (x, ...){
  max(unlist(x))
}
