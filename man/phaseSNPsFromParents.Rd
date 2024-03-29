% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phaseSNPsFromParents.R
\name{phaseSNPsFromParents}
\alias{phaseSNPsFromParents}
\title{Phases heterozygous SNPs in a patient using parents genotype}
\usage{
phaseSNPsFromParents(
  hSNPs,
  refGenome,
  mBAM = NULL,
  pBAM = NULL,
  outBase = NULL,
  FDR = 0.05,
  nullAF = 0.05,
  notNullAF = 0.3,
  alleleCounterParams = list(f = 3, F = 3852, m = 20, q = 35),
  nParallel = 1,
  ...
)
}
\arguments{
\item{hSNPs}{The heterozygous SNP locations.  See \code{\link{findHetSNPs}}.}

\item{refGenome}{Reference genome used to map \code{mBAM} and \code{pBAM}.}

\item{mBAM}{BAM file containing maternal genotype.}

\item{pBAM}{BAM file containing paternal genotype.}

\item{outBase}{Base name used to construct files storing genotyping calculations.  If left at NULL, results are not saved and cannot be reloaded in subsequent runs.}

\item{FDR}{False discovery rate for declaring an allele not present in the parental data.}

\item{nullAF}{As well as passing the FDR criteria, an allele must have frequency less than this to be declared abscent in the parent.}

\item{notNullAF}{As well as passing the FDR criteria, an allele must have frequency greater than this to be declared present in the parent.}

\item{alleleCounterParams}{Parameters to be used when running allele counter.}

\item{nParallel}{How many threads.}

\item{...}{Does nothing.}
}
\value{
\code{hSNPs} but with extra columns indicating results of phasing.
}
\description{
Given heterozygous SNPs in an individual, as called using \code{\link{findHetSNPs}}, checks genotype in Maternal and/or Paternal BAM files at these locations.  This information is used to identify which allele is maternal and paternal at each SNP.
}
\details{
To be phased, one of the alleles must be completely abscent in one of the parents.  This is determined by a binomial test against heterozygosity with FDR of \code{FDR}, plus a requirement that the allele frequency not exceed \code{nullAF}.

Where both parents are available, additionaly sanity checks are performed to make sure the trio is consistent.  This tests for presence of one allele in each parent, as measured by the allele frequency \code{notNullAF}.

Note that all three BAMs (the one used to generate \code{hSNPs}, \code{mBAM}, and \code{pBAM}) must have been mapped to the same reference genome, \code{refGenome}.
}
