/* The MIT License

   Copyright (C) 2022-2024 Giulio Genovese

   Author: Giulio Genovese <giulio.genovese@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#ifndef BCFTOOLS_INDEX_H
#define BCFTOOLS_INDEX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>

// Constants for index type
#define IDX_CSI 1  // CSI index (.csi suffix)
#define IDX_TBI 2  // Tabix index (.tbi suffix)

// Structure to hold index parameters
typedef struct {
    int min_shift;     // min_shift for CSI index (0 for TBI)
    int n_threads;     // Number of threads to use
    int index_type;    // Type of index (CSI or TBI)
    int force;         // Overwrite existing index
    const char* fname; // File to index
} index_params_t;

// Function to create an index for a VCF/BCF file
int create_vcf_index(const index_params_t* params);

#endif // BCFTOOLS_INDEX_H