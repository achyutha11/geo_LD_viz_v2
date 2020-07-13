# geo_LD_viz_v2
---------------

This is a project exploring geographic variation in LD-patterns across the
populations of the 1000 Genomes Project. 


# Environment Setup

TODO : detail steps to setup the environment via conda

# Analyses Steps

## 1. Low-Resolution Representation of LD Matrices

* Description of the K x P matrices
* Binning method for 8-bit storage
* Superpopulation vs population

### Snakemake Setup

In order to generate the banded LD matrices, we need to compute them via
`snakemake`. On midway you can use the command: 

```
./run_snakemake -s snakefiles/gen_ld_mats.snake split_collapsed_ld_mats
```

Currently this rule will generate matrices for `K=200` for chromosome 22 for
each superpopulation and population defined in the 1000 Genomes Project. The
rule automatically splits the chromosome into 10 pieces to parallelize across
so works very fast.

## 2. Fast Computation of Statistics on Low-Resolution Representation

* LD-scores
* Local LD-matrices

## 3. Plotting Statistics on a Geographic Map

* Develop plotting function via cartopy for single variants or local LD
  matrices
