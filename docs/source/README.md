# BASHer (BiAllellic SNP Haplotyper)

## What does this app do?

This script filters and classifies SNPs.

## What are typical use cases for this script?


## What inputs are required for this software to run?


## How does this script work?


## Known Limitations of this Software


### Sex checking of embryos

Different methods are used for calling informative SNPs in XL conditions depending on the sex of the embryo.  This relies on accurate sex checking of the embryos.  This should be done by the PGD team downstream of the BASHer and is outside the scope of this software - only the gene and a flanking region are loaded into BASHer and this is not enough info to draw inferences regarding the sex of the embryo.

### Sample mixup

If there has been a sample mixup and one of the trio is not related to the tested embryo’s then this could produce an incorrect diagnosis.

### Recombination

Recombination events can affect the diagnosis.  Single recombination events occur at 

### Consanguinity

Consanguineous AR cases are dealt with using a slightly modified version of the AR algorithm.  This is because they usually have less informative SNPs than non-consanguineous cases so there is a need to look at SNPs which aren't considered in other cases.  It is therefore important to identify consanguieneous cases either through the family history or via a lack of heterozygosity.

### Telomeric & Centromeric Regions

Genes are flanked with a 2mb window, for some genes in telomeric regions (FHSD1 - D4Z4 repeat, PKD1) this may mean that the window extends beyond the start/end of the chromosome. The plots produced by BASHer in these regions may look as though they have no informative SNPs over part of the window.
