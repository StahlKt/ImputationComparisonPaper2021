####functions for annotation
#this get sources by source_scripts/directories.R, which contains all other needed variables


##annotate Reference Panel
#in.file for input path to unannotated reference panel
#out.file for path to not-yet-existing annotated reference panel
set_ref_Panel <- function(in.file, out.file) {
    
  #annotate to set ID in the right format
    system2(dir_exec_bcftools, c("annotate --set-id '%CHROM\\:%POS\\_%REF\\_%FIRST_ALT'", "-o", out.file, "-O z", in.file))
  #index new file 
   system2(dir_exec_bcftools, c("index -f", out.file))
    return(out.file)
}


##annotate test data to ready for imputation course
#in.file for input path to unannotated test file
#out.file for path to not-yet-existing annotated test file
#rename.chr.file contrains old name and new names for chromomes, refer to bcftools documentation for details
set_test_data <- function(in.file, out.file, rename.chr.file) {
    
  
  #set paths for intermediate files
    filtered <- file.path(dir_input, "data_filtered.vcf.gz")
    drop <- file.path(dir_input, "data_drop.vcf.gz")
    set_chr <- file.path(dir_input, "data_filtered_set_chr.vcf.gz")
    
    #filter out SNPs that does not "PASS" the quality criteria and have less than 10% missing entries
    system2(dir_exec_bcftools, c("filter -i \"FILTER='PASS' && TYPE='snp' && F_PASS(GT='./.')<0.1\"", "-o", filtered, "-O z", in.file))
    # drop anything but biallelic SNPs
    system2(dir_exec_bcftools, c("view -m2 -M2 -v snps", "-o", drop, "-O z", filtered))
    #rename chromosomes
    system2(dir_exec_bcftools, c("annotate --rename-chrs", rename.chr.file, "-o", set_chr, "-O z", drop))
    #annotate SNPsnames
    system2(dir_exec_bcftools, c("annotate --set-id '%CHROM\\:%POS\\_%REF\\_%FIRST_ALT'", "-o", out.file, "-O z", set_chr))
    
    #remove obsolete files
    file.remove(filtered)
    file.remove(drop)
    file.remove(set_chr)
    
    #index final file
    system2(dir_exec_bcftools, c("index -f", out.file))
    
    return(out.file)
}



#annotate true data, as in full data set without the SNPs contained in the test dataset, which were not imputed
#in.file for input path to unannotated file with true genotypes for comparison (complete data set)
#out.file for path to not-yet-existing annotated genotype file for comparison
#chr.flag to indicate which chromosome 
set_true_data <- function(in.file, out.file, chr.flag){
  
  #remove output suffix for handling reasons
  out.prefix<-file_path_sans_ext(out.file,compression =TRUE)
  #name files for intermediate steps
  filtered <-file.path(dir_input,"data_filtered.vcf.gz")
  set_chr<- file.path(dir_input,"data_set_chr.vcf.gz")
  set_id<-file.path(dir_input,"dataset_id.vcf.gz")
  
  #check for chromosomes to get right files for SNP-array genotypes and files to annotate the chromosomes
  if(chr.flag==19){
    rename_chr<-dir_input_rename_chr_19
    chip_data <-dir_input_data_annotated_19
  }
  
  if(chr.flag==22){
    rename_chr<-dir_input_rename_chr_22
    chip_data <-dir_input_data_annotated_22
  }
  
  
  
  system2(
    command = dir_exec_bcftools, 
    args = c(
      "view", "-m2", "-M2", "-v snps", "-Ou", in.file,
      "|",
      dir_exec_bcftools, "annotate", "-x", "^INFO/AF,^INFO/DP,^INFO/QD,^INFO/InbreedingCoeff,^INFO/VQSLOD,^FORMAT/GT", "-Ou",
      "|",
      dir_exec_bcftools, "filter", "-i \"FILTER=\'PASS\' && TYPE=\'snp\'\"", "-Ou",
      "|",
      dir_exec_bcftools, "annotate", "--rename-chrs", rename_chr, "-Ou",
      "|",
      dir_exec_bcftools, "annotate", "--set-id \'%CHROM\\:%POS\\_%REF\\_%FIRST_ALT\'", "-Oz", "-o", set_id
    )
  )
  
  
  #index file
  system2(
    command = dir_exec_bcftools,
    args = c(
      "index", "-f", set_id
    )
  )
  
  #filter out SPNs that were already included in the test data set, which were therefore not imputed
  system2(
    command = dir_exec_bcftools,
    args = c(
      "isec", "--complement", "-Oz", "-p", out.prefix, set_id, chip_data
    )
  )
  
  #rename output files for further steps to avoid mix-ups
  file.rename(file.path(out.prefix, "0000.vcf.gz"), out.file)
  file.rename(file.path(out.prefix, "0000.vcf.gz.tbi"), paste(out.file, "tbi", sep = "."))
  
  #delete obsolete files
  unlink(out.prefix, recursive = TRUE, force = TRUE)
  
  return(out.file)
}
