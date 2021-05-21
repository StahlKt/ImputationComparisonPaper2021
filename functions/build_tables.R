#x ist directory

build_table_results<-function(x, true.data.annotated){
  ##aggr kann im init einfach auf 0 gesetzt werden, also 
  #batchreduce(fun=build_table_wrap, ergebnisse, init=0, 
  #chuncks=chunk(ergebnisse, chunk.size=1),reg=tmp )
  
  
  if(file.exists(true.data.annotated)){
    #read in
    true_table <- fread(input = true.data.annotated)
   
    #set names
    setnames(true_table, c("SAMPLE","CHROM","POS","REF","ALT", "GT"))
    
    #swap GT for coded GT
    true_table[, GT := as.numeric(factor(GT, levels = c("0/0", "0/1", "1/0", "1/1", "0|0", "0|1", "1|0", "1|1"), labels = c(0, 1, 1, 2, 0, 1, 1, 2)))-1]

    #set key to merge for later
    setkey(true_table, CHROM, POS, REF, ALT, SAMPLE)
  } else {
    stop("Annotated true data not found!")
  }
  
  #liest table ein mit samples,chrom pos, ref, alt, GP(in3 spalten)
  combined_table <- fread(
    cmd = sprintf(
      fmt = "%s %s %s",
      dir_exec_bcftools, "query -f \'[%SAMPLE,%CHROM,%POS,%REF,%ALT,%GP\\n]\'",
      x), 
    sep = ","
    )
  setnames(combined_table,c("SAMPLE","CHROM","POS","REF","ALT","0","1", "2"))
  
  #add best_guess
  combined_table[, BEST_GUESS := as.integer(colnames(.SD)[max.col(.SD, ties.method = "first")]), .SDcols = c("0", "1", "2")]
  setnames(combined_table,c("SAMPLE","CHROM","POS","REF","ALT","GP_0", "GP_1", "GP_2", "BEST_GUESS"))
  
  #add_MAX_GEN_PROB
  combined_table[, MAX_GEN_PROB:= max(GP_0,GP_1, GP_2), by=.(GP_0,GP_1,GP_2)]

  
  #get which programs were used
  col_length <- combined_table[, .N]
  programs <- get_used_programs(x)
  combined_table[, IMPUTED := programs[2]]
  combined_table[, PHASED := programs[1]]

  #get true genotypes and mafs from true_table
  setkey(combined_table, SAMPLE, CHROM, POS, REF, ALT)
  combined_table <- merge(combined_table, true_table, by=c("SAMPLE","CHROM","POS","REF","ALT"), all = FALSE)
  
  #kannn ich dataframe mit einer liste an spalten einträge aus zeilen auswählen??
  #das hier wird ewig dauern
  #scores
  combined_table[!is.na(GT), H_SCORE := hellinger_score(.SD[[sprintf("GP_%d", GT)]], GT), by = .(GT)]

  combined_table[, SEN_SCORE := sen_score(GP_0, GP_1, GP_2, GT)]
  
  return(combined_table)
}



 build_plot_tables<-function(sub.table, group.name, plot.subject){

   
   data.table(IMPUTED=sub.table[, .SD[1,IMPUTED]],
              PHASED=sub.table[, .SD[1,PHASED]],
              GROUP=group.name,
              sub.table[, as.list(summary(get(plot.subject)))]) 
   
 }
 
