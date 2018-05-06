## Introduction

`homomcl` is a C++ implementation of the BLAST similarity graph construction algorithm used by [OrthoMCL](https://github.com/stajichlab/OrthoMCL). It differs from the OrthoMCL algorithm only in that it does not include the Reciprocal Best Hits reweighting step designed to identify orthologous pairwise relationships between sequences. Thus, `homomcl` produces a similarity graph (in [ABC format](https://micans.org/mcl/man/mcl.html#started)) where edges represent simple pairwise homology. This graph can be used with [MCL](https://micans.org/mcl/) to cluster sequences into homologous families.

## Usage

`homomcl <length-file> <e-value> <blast-file>`

The input arguments are

`<length-file>` - A tab-delimited file of all sequence ids and their lengths

`<e-value>` - The e-value cutoff used to prune edges from the graph

`<blast-file>` - A BLAST results file in [tabular format (-outfmt 6)](https://github.com/seqan/lambda/wiki/BLAST-Output-Formats)