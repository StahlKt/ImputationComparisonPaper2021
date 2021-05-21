#several functions that convert data set from one format to another
#Input: in.file expects the path to the to be converted file,
#       out.file expects the path to the not-yet-existing ouput file
#       exec.XXX expects the path to the respective executable 



# vcf to bref3
vcf_bref3 <- function(in.file, exec.bref, out.file) {
    system2(dir_exec_java, c("-jar", exec.bref, in.file, ">", out.file))
    return(out.file)
}

# vcf to m3vcf
vcf_m3vcf <- function(in.file, exec.minimac, out.file) {
    system2(exec.minimac, c("--refHaps", in.file, "--processReference --prefix", 
                            file_path_sans_ext(file_path_sans_ext(out.file))))
    return(out.file)
}

# convert vcf to pbwt
vcf_pbwt <- function(in.file, exec.pbwt, out.file) {
    system2(exec.pbwt, c("-readVcfGT", in.file, "-writeAll", file_path_sans_ext(out.file)))
    return(out.file)
}



# vcf to hap/sample, used after phasing, not in the beginning for reference panels
vcf_hap <- function(data) {
    # data only prefix, exec ist bcftools
    system2(dir_exec_bcftools, c("convert --hapsample", data, paste(data, ".vcf.gz", sep = "")))
    
}
