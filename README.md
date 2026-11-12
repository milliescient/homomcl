## Introduction

`homomcl` is a C++ implementation of the BLAST similarity graph construction algorithm used by [OrthoMCL](https://github.com/stajichlab/OrthoMCL). It differs from the OrthoMCL algorithm only in that it does not include the Reciprocal Best Hits reweighting step designed to identify orthologous pairwise relationships between sequences. Thus, `homomcl` produces a similarity graph (in [ABC format](https://micans.org/mcl/man/mcl.html#started)) where edges represent simple pairwise homology. This graph can be used with [MCL](https://micans.org/mcl/) to cluster sequences into homologous families.

## Usage

`homomcl <length-file> <pml-cutoff> <blast-file>`

The input arguments are

`<length-file>` - A tab-delimited file of all sequence ids and the corresponding sequence lengths

`<pml-cutoff>` - The percent match length cutoff (0-100) used to prune edges from the graph. An edge is kept when the BLAST HSPs cover at least this percentage of the shorter sequence, in both directions. OrthoMCL's default is 50. This is not an e-value: apply any e-value cutoff when running BLAST (e.g. `-evalue 1e-5`).

`<blast-file>` - A BLAST results file in [tabular format (-outfmt 6)](https://github.com/seqan/lambda/wiki/BLAST-Output-Formats)
