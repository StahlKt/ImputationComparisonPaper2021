# ImputationComparisonPaper2021
This is the Code used in the Paper "Assessment of Imputation Quality: Comparison of Phasing and Imputation Algorithms in Real Data" from 2021


In this paper we compare the imputation quality of 35 combinations of imputation and phasing programs regarding imputation quality and speed.


Since this is the protocol for the imputation and phasing process and the cumulating of results after, this code will generally not run on other machines without changes and the needed tools and data. 


To run this code on your machine, you would need to change the file paths to fit your data structure, obtain the test data and reference panel, download the different phasing and imputation programs, create a configuration file for the R package „Batchtools“ and change the “run” command according to your system. The test data and reference panel used in this paper is not open access, but may be requested at the DZHK Omics and the HCR respectively. You can find all needed file paths you might need to change in the files contained in the source-scripts folder. All run commands are contained in the code within the folder “scripts”. For the phasing and imputation programs, please follow the links provided below. For the configuration file for batchtools, please see the batchtools documentation.


The code used to run the process are contained in the folder scripts according to the number. This is not code meant to be sourced blindly, since the single jobs may be terminated because of compatibility issues and in general the jobs of one script need to be finished successfully before the code of the next script can be run.  

1_run_prep.R: checks the existence, annotates  the used data sets and converts the reference panels before phasing and imputation takes place. 

2_run_batchtools.R: The code defines the phasing as a set of “Problems” whose results are then fed into the imputation as “Algorithms” according to batchtools.

3_run_batchtools_pbwt_19.R: PBWT for chromosome 19 is imputed as a different set of problems and algorithms, since we used different chunking because of computations reasons. 

4_run_batchmapresults.R: The created VCF files containing the imputed genotypes get extracted into large data tables and recoded for further comparison. 

5_run_reduce.R: The large tables from the last step get reduced, so that there is only one row per imputed SNP, phasing and imputation protocol with imputation quality measures and information on the SNPs. 


Beagle 4.1

https://faculty.washington.edu/browning/beagle/b4_1.html

Beagle 5.1

https://faculty.washington.edu/browning/beagle/old.beagle.html

Eagle 2.4.1

https://alkesgroup.broadinstitute.org/Eagle/

IMPUTE2 

https://mathgen.stats.ox.ac.uk/impute/impute_v2.html

IMPUTE4

https://jmarchini.org/software/

Minimac3

https://genome.sph.umich.edu/wiki/Minimac3

Minimac4

https://genome.sph.umich.edu/wiki/Minimac4

PBWT

https://github.com/richarddurbin/pbwt

SHAPEIT

https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html



Further Tools: 

BCFtools

http://samtools.github.io/bcftools/bcftools.html

bgzip

http://www.htslib.org/doc/bgzip.html

JAVA

https://java.com/
