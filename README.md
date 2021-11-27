# alleleIntegrator

An R package, `alleleIntegrator`, that integrates whole genome sequencing of tumour DNA and single cell RNA sequencing of a tumour to identify which single cells are cancer derived.  As described in [this paper](https://www.biorxiv.org/content/10.1101/2021.11.25.469995v1).

## Installation

The R package should be installable in the usual way with `devtools`

```R
devtools::install_github('constantAmateur/alleleIntegrator')
```

## Dependencies

alleleIntegrator depends on a few external tools that need to be installed for it to function fully.

[alleleCount](https://github.com/cancerit/alleleCount) - Used to count reads supporting genotype in DNA and RNA bam files.  There is a pre-compiled binary in `inst/bin` that should work on ubuntu, otherwise compile using instructions in the alleleCount repository.

bcftools - Used to call variants in DNA bam files.  On ubuntu run `sudo apt install bcftools`.

[ASCAT](https://github.com/VanLoo-lab/ascat) used for calling copy number state from whole genome sequencing.  Pretty much all of alleleIntegrator will work without this.  Install with `devtools::install_github('VanLoo-lab/ascat/ASCAT')`

## Usage

There are four basic steps to using alleleIntegrator to identify cancer transcriptomes:

### 1 - Specify copy number profile

You would usually obtain this from running a copy number caller, such as [ASCAT](https://github.com/VanLoo-lab/ascat).  alleleIntegrator needs a `GRanges` object with chromosome segments and columns `matNum`, `patNum` and `tumFrac = matNum/(matNum+patNum)`.

### 2 - Identify heterozygous SNPs

Use the `findHetSNPs` function on the whole genome sequencing BAM.  These can also be manually provided if you know them from some other source.

### 3 - Phase heterozygous SNPs

Use the `phaseSNPsFromCN` function on the whole genome sequencing of tumour BAM.  You will also need to provide the heterozygous SNPs and copy number segment definitions determined in the previous steps.

### 4 - Fit model

Finally, the 10X BAM files are interrogated to determine which allele is expressed at each phased SNP in each copy number segment using `getAllelicExpression`.  These counts should then be filtered using `filterCells` and the model fit by running `calcASE`, `calcOverDispersion` and `abbSegProb`.  

## Tutorial

Look at `exampleRun.R` for a demonstration how to use `alleleIntegrator` to identify cancer cell transcriptomes.  This example starts from whole genome sequencing of tumour/blood from a Neuroblastoma from [this paper](https://www.science.org/doi/10.1126/sciadv.abd3311) along with cellranger mapped BAMs containing 10X 3' expression from the tumour.  

This example should provide a template for how to apply this package.  However, it does have some eccentricities that are specific to this sample.  In particular, the DNA BAM files have been mapped against hg19, while the scRNA data is mapped against GRCh38.  As such, it is necessary to perform a liftover to convert phased, heterozygous SNP coordinates between the references.

As an optional first step of this example, the `matchBAMs` function is used to compare the genotype of all DNA and RNA BAM files.  This step is completely optional, but experience indicates that sample mixups are depressingly common and it is sensible to verify that the samples that you specify really are related in the way that you expect.

Also demonstrated in this example is how to aggregate cells into clusters and apply inference at the cluster level.  Usually there is sufficient information to identify the presence/absence of the cancer genotype in each cell.  But in some circumstances, it can be useful to assume that all cells in a cluster are of the same genotype and treat each cluster as if it were a single cell.  If the cells in a cluster are of various genotypes, this will no improve anything, but if this assumption holds than the power to detect allelic imbalances is boosted by aggregating counts across all cells in a cluster.

## Citation

If you use this package in your work, please cite: Trinh, M.K, et al., Precise identification of cancer cells from allelic imbalances in single cell transcriptomes, BiorXiv, 2021
