#Shut up about global variables.
utils::globalVariables(c("ExAcSNPs",'biasedGenes','MAF','covBins','tumFrac','SNPs1k'))

#This is how to process the imprinted genes data.table to impGenes
#x = read.table('~/trueHome/Projects/Common/Data/imprintedGenes.txt',sep='\t',header=TRUE)
#impGenes = unique(x[x$Species=='Homo sapiens' & x$Status=='Imprinted',]$Gene)
#impGenes = gsub(' V2$','',gsub('[@*]$','',impGenes))
#Full definition of biasedGenes
#biasedGenes = list(MHCI = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-K", "HLA-L", "HLA-N", "HLA-P", "HLA-S", "HLA-T", "HLA-U", "HLA-V", "HLA-W", "HLA-X", "HLA-Y", "HLA-Z"),
#                   MHCII = c("HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPA2", "HLA-DPA3", "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DQB3", "HLA-DRA", "HLA-DRB1", "HLA-DRB2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8", "HLA-DRB9"),
#                   HB = c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"),
#                   IG = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC"),
#                   MT = c("MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB"),
#                   Imprinted = c("ADTRP", "AIM1", "ANO1", "ATP10A", "ATP5F1EP2", "BLCAP", "CCDC71L", "CDKN1C", "CMTM1", "COPG2IT1", "CPA4", "DDC", "DGCR6", "DGCR6L", "DIO3", "DIO3OS", "DIRAS3", "DLGAP2", "DLK1", "DLX5", "DNMT1", "DSCAM", "ERAP2", "ESR2", "FAM50B", "GDAP1L1", "GLI3", "GLIS3", "GNAS", "GNASAS", "GPR1", "GRB10", "H19", "HECW1", "HNF1A", "HOXA4", "HYMAI", "IGF2", "IGF2AS", "INPP5F", "INS", "IRAIN", "KCNK9", "KCNQ1", "KCNQ1DN", "KCNQ1OT1", "KLF14", "L3MBTL1", "LIN28B", "LRRTM1", "MAGEL2", "MAGI2", "MCTS2", "MEG3", "MEG8", "MEST", "MESTIT1", "MIMT1", "MIR296", "MIR298", "MIR371A", "MKRN3", "NAA60", "NAP1L5", "NDN", "NLRP2", "NNAT", "NPAP1", "NTM", "OSBPL5", "PEG10", "PEG13", "PEG3-AS1", "PEG3", "PHLDA2", "PLAGL1", "PPP1R9A", "PSIMCT-1", "PWAR6", "PWCR1", "PXDC1", "RASGRF1", "RB1", "RBP5", "RHOBTB3", "RNU5D-1", "RTL1", "SANG", "SGCE", "SGK2", "SLC22A18", "SLC22A2", "SLC22A3", "SMOC1", "SNORD107", "SNORD108", "SNORD109A", "SNORD109B", "SNORD113-1", "SNORD114-1", "SNORD115-48", "SNORD115", "SNORD116", "SNORD64", "SNRPN", "SNURF", "ST8SIA1", "SVOPL", "TCEB3C", "TFPI2", "TP53", "TP73", "UBE3A", "VTRNA2-1", "WT1-AS", "WT1", "ZC3H12C", "ZDBF2", "ZFAT-AS1", "ZFAT", "ZFP90", "ZIM2", "ZNF396", "ZNF597")
#                   )
#biasedGenes = unique(unlist(biasedGenes,use.names=FALSE))

