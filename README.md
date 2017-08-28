# HLA-check

##Summary

HLA-Check: Evaluating HLA data from SNP information

This is the official repository for the HLA-Check tool. The accompanying paper is still under review.

##Compiling

The code only depends on OpenMP:
`g++ -O3 judge_sample.cpp -fopenmp -lgomp -std=c++11 -o hla-check`

##Usage

`./hla-check <HLA> <imputed_file> <hla_file> <files/ path>`

Where:
* HLA is one in A, B, C, DRB1, DPA1, DPB1, DQA1, DQB1.
* imputed_file is the output file from IMPUTE2
* hla_file is a file containing HLA typings for individuals in the IMPUTE2 .sample file format: columns should be Family_id, Individual_id, 0,0,2,-9, then hla alleles in 4-digit, using 0 for missing data: `700` to code for 2-digit "07" type and `701` to code for 07:01 2-filed typing. Every individual should have typings for all class 1 and class 2 alleles, ordered in the chromosome order: HLA-A1 HLA-C, HLA-B, HLA-DRB, HLA-DQA, HLA-DQB, HLA-DPA, HLA-DPB. 
* files/ is the path of the files/ folder of this repo (if you're running the program from this folder, just put "files/")

##Included files

The `files/` folder contains files used internally by hla-check.
* HLA alleles alignments, generated from the IPD-IMGT/HLA Database (hla_nuc)
* Manual alignments between the reference allele and the NM sequence
* List of annotated SNPs from RefSeq in HLA exons, with their NM and their position (downloaded from SNP-Nexus)

##Planned improvements

* Take alleles in the colon-separated fields format to support >99 hla alleles
* Add header management to ease file manipulation


