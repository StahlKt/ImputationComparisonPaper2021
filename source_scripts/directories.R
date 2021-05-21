## Directories

library("tools")
library("batchtools")



# home directory/working dir and other folders
home_dir <- "/imbs/home/stahl/MA"

#directory to functions folder
dir_functions <- file.path(home_dir, "functions")

#directories for input and output data
dir_data <- "/imbs/projects/students/stahl/Masterarbeit_Stahl/data"
dir_input <- file.path(dir_data, "input")
dir_output <- file.path(dir_data, "output_new")

#directory for the registries that will be created by batchtools
dir_registry <- file.path(dir_data, "registries")

#directory for folder with executables
dir_exec <- file.path("/software/bin")

#directory for source_scripts folder
source_dir<-file.path(home_dir,"source_scripts")

#directory for batchtools configuration file
confic_dir<-"./batchtools.conf.R"

### Directories for required data
#input data for imputation
dir_input_data_raw_19 <- file.path(dir_input, "dzhk_omics_infiniumomni5_chr19_neu.vcf.gz")
dir_input_data_raw_22 <- file.path(dir_input, "dzhk_omics_infiniumomni5_chr22_neu.vcf.gz")

# reference panels in vcf
dir_input_reference_vcf_raw_19 <- file.path(dir_input, "HRC.r1-1.EGA.GRCh37.chr19.haplotypes.vcf.gz")
dir_input_reference_vcf_raw_22 <- file.path(dir_input, "HRC.r1-1.EGA.GRCh37.chr22.haplotypes.vcf.gz")

# reference panels in hap (no suffixes for the files making up the .hap format, prefix must therefore be the same)
dir_input_reference_hap_19 <- file.path(dir_input, "HRC.r1-1.EGA.GRCh37.chr19")
dir_input_reference_hap_22 <- file.path(dir_input, "HRC.r1-1.EGA.GRCh37.chr22")

# map files with a column for chromosomes 
dir_input_map_mit_chr_19 <- file.path(dir_input, "genetic_map_chr19_combined_b37_added_chr.txt")
dir_input_map_mit_chr_22 <- file.path(dir_input, "genetic_map_chr22_combined_37_added_chr.txt")

#  map files without a column for chromosomes 
dir_input_map_ohne_chr_19 <- file.path(dir_input, "genetic_map_chr19_combined_b37.txt")
dir_input_map_ohne_chr_22 <- file.path(dir_input, "genetic_map_chr22_combined_b37.txt")

# map files in plink format
dir_input_map_plink_19 <- file.path(dir_input, "plink.chr19.GRCh37.map")
dir_input_map_plink_22 <- file.path(dir_input, "plink.chr22.GRCh37.map")

# file for renaming chromosomes and SNP names, used in annotation
dir_input_rename_chr_19 <- file.path(dir_input, "chr19.map")
dir_input_rename_chr_22 <- file.path(dir_input, "chr22.map")



#directories to the needed executables
# 
dir_exec_shapeit <- file.path(dir_exec, "shapeit.v2.r904.GLIBCv2.17")
# 
dir_exec_eagle <- file.path(dir_exec, "eagle_v2.4.1")
# 
dir_exec_impute2 <- file.path(dir_exec, "impute2.3.2")
# 
dir_exec_impute4 <- file.path(dir_exec, "impute4.1.2_r300.2")
# 
dir_exec_beagle4 <- "/software/src/beagle/beagle.11Mar19.69c.jar"
# 
dir_exec_beagle5 <- "/software/src/beagle/beagle.21Sep19.ec3.jar"
# 
dir_exec_minimac3 <- file.path(dir_exec, "minimac_v3.2.0.1")
# 
dir_exec_minimac4 <- file.path(dir_exec, "minimac_v4.1.0.0")
# 
dir_exec_pbwt <- file.path(dir_exec, "pbwt_3.0-64c4ffa")
# 
dir_exec_bref <- file.path(dir_exec, "bref_v3.20190516")
# 
dir_exec_bcftools <- file.path(dir_exec, "bcftools1.9")
dir_exec_bgzip <- file.path(dir_exec, "bgzip1.9")
dir_exec_java <- "/usr/bin/java"



## The following directories do not have to be changes necessarily, but they dictate
# what the newly created files are saved as
dir_input_data_annotated_19 <- file.path(dir_input, "data_annotated_19.vcf.gz")
dir_input_data_annotated_22 <- file.path(dir_input, "data_annotated_22.vcf.gz")

dir_input_reference_vcf_annotated_19 <- file.path(dir_input, "refPanel_annotated_19.vcf.gz")
dir_input_reference_vcf_annotated_22 <- file.path(dir_input, "refPanel_annotated_22.vcf.gz")

dir_input_exclude_shapeit_19 <- file.path(dir_input, "unaligned_snps_19")
dir_input_exclude_shapeit_22 <- file.path(dir_input, "unaligned_snps_22")

dir_input_samples_19 <- file.path(dir_input, "sample_list_pbwt_19.txt")
dir_input_samples_22 <- file.path(dir_input, "sample_list_pbwt_22.txt")

dir_input_reference_pbwt_19 <- file.path(dir_input, "refPanel_19.pbwt")
dir_input_reference_pbwt_22 <- file.path(dir_input, "refPanel_22.pbwt")

dir_input_reference_bref3_19 <- file.path(dir_input, "refPanel_19.bref3")
dir_input_reference_bref3_22 <- file.path(dir_input, "refPanel_22.bref3")

dir_input_reference_m3vcf_19 <- file.path(dir_input, "refPanel_19.m3vcf.gz")
dir_input_reference_m3vcf_22 <- file.path(dir_input, "refPanel_22.m3vcf.gz")




### directiries for output folders
#if these folders don't exist, they will be created

dir_out_shapeit <- file.path(dir_output, "phasing", "Shapeit2")
dir_out_eagle <- file.path(dir_output, "phasing", "Eagle2")
dir_out_beagle5_phasing <- file.path(dir_output, "phasing", "Beagle5")
dir_out_impute2 <- file.path(dir_output, "imputation", "Impute2")
dir_out_impute4 <- file.path(dir_output, "imputation", "Impute4")
dir_out_minimac3 <- file.path(dir_output, "imputation", "Minimac3")
dir_out_minimac4 <- file.path(dir_output, "imputation", "Minimac4")
dir_out_beagle4 <- file.path(dir_output, "imputation", "Beagle4.1")
dir_out_beagle5_imputation <- file.path(dir_output, "imputation", "Beagle5")
dir_out_pbwt <- file.path(dir_output, "imputation", "PBWT")



# global variables: chunks
chunk_22 <- 4406250
list_of_chunks_22 <- seq(16050000, 16050000 + 7 * chunk_22, chunk_22)

chunk_19 <- 4549231
list_of_chunks_19 <- seq(60000, 60000 + 12 * chunk_19, chunk_19)



#source functions for imputation
source(file.path(dir_functions,"annotate_functions.R"))
source(file.path(dir_functions,"beagle4.R"))
source(file.path(dir_functions,"beagle5.R"))
source(file.path(dir_functions,"check_folder_function.R"))
source(file.path(dir_functions,"convert_functions.R"))
source(file.path(dir_functions,"eagle.R"))
source(file.path(dir_functions,"impute2.R"))
source(file.path(dir_functions,"impute4.R"))
source(file.path(dir_functions,"minimac3.R"))
source(file.path(dir_functions,"minimac4.R"))
source(file.path(dir_functions,"pbwt.R"))
source(file.path(dir_functions,"shapeit.R"))
