#' A simple example showing how to use alleleIntegrator
#'
#' Last tested on v0.7.3

#############
# Libraries #
#############

library(alleleIntegrator)
library(GenomicFeatures)
library(Seurat)

##########
# Params #
##########

outDir = 'Results/Genotyping/'
refGenome = 'Data/DNA/genomeDNA.fa'
refGenome10X = 'Data/scRNA/genomeRNA.fa'
liftChain = 'Data/hg19ToHg38_noChr.over.chain'
gtf = 'Data/gtf10X_GRCh38_120.gtf'
nParallel=12

###########
# PD46693 #
###########

#########################
# Sample specific params
tumourDNA = 'Data/DNA/PD46693a.bam'
patientDNA = 'Data/DNA/PD46693b.bam'
bams10X = list.files('Data/scRNA/PD46693',full.names=TRUE)
bams10X = setNames(file.path(bams10X,'possorted_genome_bam.bam'),basename(bams10X))
#Get annotated cells
srat = readRDS('Data/srat_inhouse.rds')
srat@meta.data$cellID = paste0('4602',rownames(srat@meta.data),'-1')
normCells = srat@meta.data$cellID[!grepl('tumour',srat@meta.data$idents_for_plot)]
#Define CN segments roughly
altChrs = c(1,4,5,6,9,10,11,12,13,15,16,20,21,22)
segs = GRanges(altChrs,IRanges(rep(1,length(altChrs)),1e9))
segs$matNum = c(2,2,2,2,2,2,2,2,2,2,2,2,2,2)
segs$patNum = c(1,0,1,1,1,0,0,1,1,1,1,1,0,1)
segs$tumFrac = segs$matNum/(segs$patNum+segs$matNum)
segs$normFrac = 0.5
names(segs) = altChrs
#############################
# Check genotype consistency
#Are all the BAMs you're going to use from the same individual?  Check before you start
genoCheck = matchBAMs(BAMs = c(norm=patientDNA,tum=tumourDNA,bams10X),
                    refGenomes = rep(c(refGenome,refGenome10X),c(2,3)),
                    outputs = file.path(outDir,paste0(c(basename(c(patientDNA,tumourDNA)),names(bams10X)),'_genotypeCheck.tsv')),
                    liftOvers=rep(c(NA,liftChain),c(2,3)),
                    is10X=rep(c(FALSE,TRUE),c(2,3)),
                    nParallel=nParallel)
#If anything is less than 0.8 and you should be concerned...
message(sprintf("The minimum similarity found was %g",min(genoCheck$ibs$ibs)))
######################
# Call and phase SNPs
hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,'PD46693_patient_hetSNPs.vcf'),nParallel=nParallel)
#Expectation is that we'll find ~ 3 million of them
message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
#Use tumour DNA to phase them.  Set up plot area so we can inspect the CN changes
par(mfrow=c(5,4))
phSNPs = phaseSNPsFromCN(hSNPs,segs,refGenome,tumourDNA,outPath=file.path(outDir,'PD46693_tumour_countAtHetSNPs.tsv'),nParallel=nParallel)
#Liftover to GRCh38
phSNPs38 = changeGenomeVersion(phSNPs,liftChain)
#Annotate SNPs using GTF
phSNPs38 = annotateSNPs(phSNPs38,gtf)
########################
# Integrate with 10X.  
#If the majority of the high coverage SNPs don't look heterozygous, something has gone wrong...
phCnts = getAllelicExpression(phSNPs38,refGenome10X,bams10X,outputs=file.path(outDir,paste0('PD46693_',names(bams10X),'_scRNA_alleleCounts.tsv')),nParallel=nParallel)
##############
# Filter data
#You don't **have** to provide cluster information here.  You could just tell filterCells which bacodes represent cells by specificying "passCellIDs".  But clustering information helps a lot with interpretation
clusterIDs = setNames(srat@meta.data$idents_for_plot,srat@meta.data$cellID)
#If you don't have much power at the individual cell level, you could consider collapsing all counts into small clusters and then treating each cluster as a cell.  The code below demonstrates how to do this.  After running aggregateByClusters you can proceed in the same way as if there had been no aggregation.
aggToClust=FALSE
if(aggToClust){
gCnts = aggregateByClusters(phCnts,clusterIDs)
gCnts = filterCells(gCnts,passCellIDs=levels(clusterIDs),normIDs=c('immune','vascular endothelium'))
}else{
gCnts = filterCells(phCnts,clusterIDs=clusterIDs,normIDs=c('immune','vascular endothelium'))
}
##################
# Calibrate model
#Specify the error rate
gCnts$errRate = c('Exonic'=0.01,'Intronic'=0.05,'Intergenic'=0.15)[gCnts$regionType]
#Detect allele specific expression
gCnts = calcASE(gCnts)
#Get over-dispersion
od = calcOverDispersion(gCnts)
############
# Inference
pp = abbSegProb(gCnts,od)
#############
# Validation
dat = plotRawData(gCnts,returnData=TRUE)
plotPosteriorHeatmap(pp,'nLL')
#Integrate global call with Seurat, dropping chr4 with sub-clone
pp = abbSegProb(gCnts,od,segs = gCnts@metadata$segs[-2],abbFrac = 'tumFrac',globalOnly=TRUE)
m = match(pp$cellID,srat@meta.data$cellID)
srat@meta.data$call = NA
srat@meta.data$call[m] = ifelse(pp$maxPostProb>0.8,pp$mostLikelyState,'Uncalled')
DimPlot(srat,group.by='call')
table(srat@meta.data$call,srat@meta.data$idents_for_plot)

