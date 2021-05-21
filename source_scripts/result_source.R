##source for results

library("R.utils")
library("batchtools")
library("tools")

home_dir <- "/imbs/home/stahl/MA"

source(file.path(home_dir, "source_scripts/directories.R"))
source(file.path(home_dir, "functions/score_and_coding_functions.R"))
source(file.path(home_dir, "functions/build_tables2.R"))
source(file.path(home_dir, "functions/reduce_functions.R"))




### Directories for required data
#input data for comparison, complete true genotypes
dir_input_data_true_raw_19<-file.path(dir_input,"dzhk_omics_chr19_neu.vcf.gz")
dir_input_data_true_raw_22<-file.path(dir_input,"dzhk_omics_chr22_neu.vcf.gz")


###
#directory where to save annotated complete "true" data
dir_input_data_true_annotated_19<-file.path(dir_input,"dzhk_omics_chr19_neu_annotated.vcf.gz")
dir_input_data_true_annotated_22<-file.path(dir_input,"dzhk_omics_chr22_neu_annotated.vcf.gz")

#split into sample and variant info for easier handling
dir_input_data_true_annotated_sample_genotypes_19<-file.path(dir_input,"dzhk_omics_chr19_neu_annotated_sample.genotypes")
dir_input_data_true_annotated_variant_infos_19<-file.path(dir_input,"dzhk_omics_chr19_neu_annotated_variant.infos")

dir_input_data_true_annotated_sample_genotypes_22<-file.path(dir_input,"dzhk_omics_chr22_neu_annotated_sample.genotypes")
dir_input_data_true_annotated_variant_infos_22<-file.path(dir_input,"dzhk_omics_chr22_neu_annotated_variant.infos")

#dir_out_plot_table<-file.path(dir_output,"plot_table")
#dir_out_score_table<-file.path(dir_output,"score_table")
#dir_out_error_table<-file.path(dir_output,"error_table")



chunk_22 <- 4406250
list_of_chunks_22 <- seq(16050000, 16050000 + 7 * chunk_22, chunk_22)

chunk_19 <- 4549231
list_of_chunks_19 <- seq(60000, 60000 + 12 * chunk_19, chunk_19)

#output directory for reduced data tables.
dir_out_table_22<-file.path(dir_output, "reduced_tables_22.rds")
dir_out_table_19<-file.path(dir_output, "reduced_tables_19.rds")


