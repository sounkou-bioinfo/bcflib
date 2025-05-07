#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include "BcftoolsMerge.h"

// Include the declaration for the bcftools main_vcfmerge function
// but don't modify the original code
extern "C" {
    int main_vcfmerge(int argc, char *argv[]);
}

// Function to run bcftools merge from command line arguments
int run_bcftools_merge(int argc, char **argv) {
    return main_vcfmerge(argc, argv);
}

// Function to build and run the bcftools merge command from parameters struct
int run_bcftools_merge_direct(const merge_params_t* params) {
    if (!params) return -1;

    // Count the number of arguments needed
    int n_args = 2; // Program name + at least one file
    
    // Count options
    if (params->record_cmd_line == 0) n_args += 1; // --no-version
    if (params->write_index) n_args += 1; // --write-index
    if (params->output_type != FT_VCF || params->clevel >= 0) n_args += 2; // -O option
    if (params->n_threads > 0) n_args += 2; // --threads option
    if (params->missing_to_ref) n_args += 1; // -0 option
    if (params->force_samples) n_args += 1; // --force-samples
    if (params->header_only) n_args += 1; // --print-header
    if (params->header_fname) n_args += 2; // --use-header option
    if (params->regions_list) n_args += 2; // -r option
    if (params->info_rules) n_args += 2; // -i option
    if (params->missing_rules) n_args += 2; // -M option
    if (params->output_fname && strcmp(params->output_fname, "-") != 0) n_args += 2; // -o option
    if (params->gvcf_fai_fname) n_args += 2; // -g option
    if (params->local_alleles > 0) n_args += 2; // -L option
    if (params->merge_by_id || params->collapse != COLLAPSE_BOTH) n_args += 2; // -m option
    if (params->filter_logic) n_args += 2; // -F option
    if (params->force_single) n_args += 1; // --force-single
    
    // Add one argument for each input file
    n_args += params->n_input_files;

    // Allocate arguments array
    char** args = (char**)malloc(n_args * sizeof(char*));
    if (!args) return -1;
    
    // First argument is the program name
    args[0] = strdup("bcftools");
    
    // Initialize argument counter
    int arg_idx = 1;
    
    // Add options
    if (params->record_cmd_line == 0)
        args[arg_idx++] = strdup("--no-version");
    
    if (params->write_index)
        args[arg_idx++] = strdup("--write-index");
    
    if (params->output_type != FT_VCF || params->clevel >= 0) {
        args[arg_idx++] = strdup("-O");
        char output_type[8];
        
        switch (params->output_type) {
            case FT_BCF_GZ: output_type[0] = 'b'; break;
            case FT_BCF: output_type[0] = 'u'; break;
            case FT_VCF_GZ: output_type[0] = 'z'; break;
            default: output_type[0] = 'v'; break;
        }
        
        if (params->clevel >= 0) {
            output_type[1] = '0' + params->clevel;
            output_type[2] = '\0';
        } else {
            output_type[1] = '\0';
        }
        args[arg_idx++] = strdup(output_type);
    }
    
    if (params->n_threads > 0) {
        args[arg_idx++] = strdup("--threads");
        char threads_str[16];
        snprintf(threads_str, sizeof(threads_str), "%d", params->n_threads);
        args[arg_idx++] = strdup(threads_str);
    }
    
    if (params->missing_to_ref)
        args[arg_idx++] = strdup("-0");
    
    if (params->force_samples)
        args[arg_idx++] = strdup("--force-samples");
    
    if (params->header_only)
        args[arg_idx++] = strdup("--print-header");
    
    if (params->header_fname) {
        args[arg_idx++] = strdup("--use-header");
        args[arg_idx++] = strdup(params->header_fname);
    }
    
    if (params->regions_list) {
        args[arg_idx++] = strdup("-r");
        args[arg_idx++] = strdup(params->regions_list);
    }
    
    if (params->info_rules) {
        args[arg_idx++] = strdup("-i");
        args[arg_idx++] = strdup(params->info_rules);
    }
    
    if (params->missing_rules) {
        args[arg_idx++] = strdup("-M");
        args[arg_idx++] = strdup(params->missing_rules);
    }
    
    if (params->output_fname && strcmp(params->output_fname, "-") != 0) {
        args[arg_idx++] = strdup("-o");
        args[arg_idx++] = strdup(params->output_fname);
    }
    
    if (params->gvcf_fai_fname) {
        args[arg_idx++] = strdup("-g");
        args[arg_idx++] = strdup(params->gvcf_fai_fname);
    }
    
    if (params->local_alleles > 0) {
        args[arg_idx++] = strdup("-L");
        char local_str[16];
        snprintf(local_str, sizeof(local_str), "%d", params->local_alleles);
        args[arg_idx++] = strdup(local_str);
    }
    
    if (params->merge_by_id || params->collapse != COLLAPSE_BOTH) {
        args[arg_idx++] = strdup("-m");
        char* collapse_str;
        
        if (params->merge_by_id) {
            collapse_str = strdup("id");
        } else {
            switch (params->collapse) {
                case COLLAPSE_SNPS: collapse_str = strdup("snps"); break;
                case COLLAPSE_INDELS: collapse_str = strdup("indels"); break;
                case COLLAPSE_NONE: collapse_str = strdup("none"); break;
                case COLLAPSE_ANY: collapse_str = strdup("all"); break;
                case COLLAPSE_SNP_INS_DEL: collapse_str = strdup("snp-ins-del"); break;
                default: collapse_str = strdup("both"); break;
            }
        }
        
        if (params->trim_star_allele == 1) {
            // Append * to the collapse string
            char* new_str = (char*)malloc(strlen(collapse_str) + 2);
            strcpy(new_str, collapse_str);
            strcat(new_str, "*");
            free(collapse_str);
            collapse_str = new_str;
        } else if (params->trim_star_allele > 1) {
            // Append ** to the collapse string
            char* new_str = (char*)malloc(strlen(collapse_str) + 3);
            strcpy(new_str, collapse_str);
            strcat(new_str, "**");
            free(collapse_str);
            collapse_str = new_str;
        }
        
        args[arg_idx++] = collapse_str;
    }
    
    if (params->filter_logic) {
        args[arg_idx++] = strdup("-F");
        args[arg_idx++] = strdup(params->filter_logic == FLT_LOGIC_ADD ? "+" : "x");
    }
    
    if (params->force_single)
        args[arg_idx++] = strdup("--force-single");
    
    // Add input files
    for (int i = 0; i < params->n_input_files; i++) {
        args[arg_idx++] = strdup(params->input_fnames[i]);
    }
    
    // Make sure we've used the correct number of arguments
    if (arg_idx != n_args) {
        fprintf(stderr, "Error: argument count mismatch (expected %d, got %d)\n", n_args, arg_idx);
        for (int i = 0; i < arg_idx; i++)
            free(args[i]);
        free(args);
        return -1;
    }
    
    // Execute the main_vcfmerge function
    int ret = main_vcfmerge(arg_idx, args);
    
    // Free allocated memory
    for (int i = 0; i < arg_idx; i++)
        free(args[i]);
    free(args);
    
    return ret;
}