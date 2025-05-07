#ifndef RCPP_BCFTOOLS_MERGE_H
#define RCPP_BCFTOOLS_MERGE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>

// Constant for output format
#define FT_BCF        1
#define FT_VCF        2
#define FT_VCF_GZ     4
#define FT_BCF_GZ     8

// Constants for collapse mode
#define COLLAPSE_NONE       0
#define COLLAPSE_SNPS       1
#define COLLAPSE_INDELS     2
#define COLLAPSE_BOTH       (COLLAPSE_SNPS|COLLAPSE_INDELS)
#define COLLAPSE_ANY        4
#define COLLAPSE_SOME       8
#define COLLAPSE_SNP_INS_DEL (1<<10)

// Constants for filter logic
#define FLT_LOGIC_ADD    0
#define FLT_LOGIC_REMOVE 1

// Structure to hold vcfmerge parameters
typedef struct {
    int record_cmd_line;         // append version and command line to header
    int write_index;             // automatically index output files
    int output_type;             // output file type
    int clevel;                  // compression level
    int n_threads;               // number of worker threads
    int missing_to_ref;          // missing genotypes are treated as reference
    int force_samples;           // resolve duplicate sample names
    int header_only;             // only print the merged header
    int merge_by_id;             // merge sites by ID
    int collapse;                // merge sites with identical position
    int trim_star_allele;        // trim unobserved star alleles
    int local_alleles;           // when to use local alleles
    int keep_AC_AN;              // keep AC/AN fields
    int force_single;            // allow single file
    int filter_logic;            // filter logic for overlapping records
    int no_index;                // merge unindexed files
    const char* output_fname;    // output file name
    const char* header_fname;    // custom header file
    const char* regions_list;    // regions to restrict merge to
    const char* info_rules;      // rules for merging INFO fields
    const char* missing_rules;   // rules for replacing missing values
    const char* gvcf_fai_fname;  // reference fasta for gVCF
    const char** input_fnames;   // array of input file names
    int n_input_files;           // number of input files
} merge_params_t;

#ifdef __cplusplus
extern "C" {
#endif

// Function to run bcftools merge from command line arguments
int run_bcftools_merge(int argc, char** argv);

// Function to run bcftools merge directly with parameters struct
int run_bcftools_merge_direct(const merge_params_t* params);

#ifdef __cplusplus
}
#endif

#endif // RCPP_BCFTOOLS_MERGE_H