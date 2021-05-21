# run prep
#this is the first file to execute to run batchtools


## Set working directory to folder that contains following folders
#"source_scripts" containint the R files from the github repository.
#'exec' with the executables
#'data' containing another folder 'input' with inputfiles. 
#the file paths in directories.R need to be set according to the position of files in your system 

install.packages('batchtools') 
install.packages('snow')
install.packages("tools")
install.packages("data.table")

# home directory/working dir and other folders
home_dir <- "/imbs/home/stahl/MA"
setwd(home_dir)

#source all needed directories and functions, needs to be adjusted
source_dir<-file.path(home_dir,"source_scripts")
source(file.path(source_dir,"directories.R"))
source(file.path(source_dir,"result_source.R"))

# checks output folders, if creates them if they don't exist
sapply(mget(ls(pattern = "out")), check_and_create_folders)

# checks if executables are in place
sapply(mget(ls(pattern = "exec")), file.exists)
# check if needed data is in place


# annotateting reference panel and test data
set_test_data(dir_input_data_raw_19, dir_input_data_annotated_19, dir_input_rename_chr_19)
set_test_data(dir_input_data_raw_22, dir_input_data_annotated_22, dir_input_rename_chr_22)

set_ref_Panel(dir_input_reference_vcf_raw_19, dir_input_reference_vcf_annotated_19)
set_ref_Panel(dir_input_reference_vcf_raw_22, dir_input_reference_vcf_annotated_22)


# create samplelist for pbwt_prep
pbwt_get_samples(dir_input_data_annotated_22, dir_input_samples_22)
pbwt_get_samples(dir_input_data_annotated_19, dir_input_samples_19)


# generates list of SNPs not well aligned with refPanel, to be excluded in Shapeit
shapeit_prep(dir_input_data_annotated_19, dir_exec_shapeit, dir_input_exclude_shapeit_19, dir_input_reference_hap_19)
shapeit_prep(dir_input_data_annotated_22, dir_exec_shapeit, dir_input_exclude_shapeit_22, dir_input_reference_hap_22)

#generate reference panels in all needed formats
vcf_pbwt(dir_input_reference_vcf_annotated_19, dir_exec_pbwt, dir_input_reference_pbwt_19)
vcf_bref3(dir_input_reference_vcf_annotated_19, dir_exec_bref, dir_input_reference_bref3_19)
vcf_m3vcf(dir_input_reference_vcf_annotated_19, dir_exec_minimac3, dir_input_reference_m3vcf_19)

vcf_pbwt(dir_input_reference_vcf_annotated_22, dir_exec_pbwt, dir_input_reference_pbwt_22)
vcf_bref3(dir_input_reference_vcf_annotated_22, dir_exec_bref, dir_input_reference_bref3_22)
vcf_m3vcf(dir_input_reference_vcf_annotated_22, dir_exec_minimac3, dir_input_reference_m3vcf_22)


# check again if needed data is in place
#length(ls(pattern = "input")) == 15
#ls(pattern = "input")
#sapply(mget(ls(pattern = "input")), file.exists)





#annotate sequenced data sets
set_true_data(dir_input_data_true_raw_19,dir_input_data_true_annotated_19)
set_true_data(dir_input_data_true_raw_22,dir_input_data_true_annotated_22,22)

#divide datasets into one for genotypes only for imputation quality comparison and one for information on the SNPs
system2(
    command = dir_exec_bcftools, 
    args = c("query", 
             "-f \'[%SAMPLE %CHROM %POS %REF %ALT %GT\\n]\'", 
             dir_input_data_true_annotated_19,
             ">", dir_input_data_true_annotated_sample_genotypes_19)
)

system2(
    command = dir_exec_bcftools, 
    args = c("query", 
             "-f \'%CHROM %POS %REF %ALT %INFO/VQSLOD %INFO/QD %INFO/InbreedingCoeff %INFO/DP %INFO/AF\\n\'", 
             dir_input_data_true_annotated_19,
             ">", dir_input_data_true_annotated_variant_infos_19)
)
 system2(
    command = dir_exec_bcftools, 
    args = c("query", 
             "-f \'[%SAMPLE %CHROM %POS %REF %ALT %GT\\n]\'", 
             dir_input_data_true_annotated_22,
             ">", dir_input_data_true_annotated_sample_genotypes_22)
)
system2(
    command = dir_exec_bcftools, 
    args = c("query", 
             "-f \'%CHROM %POS %REF %ALT %INFO/VQSLOD %INFO/QD %INFO/InbreedingCoeff %INFO/DP %INFO/AF\\n\'", 
             dir_input_data_true_annotated_22,
             ">", dir_input_data_true_annotated_variant_infos_22)
)





