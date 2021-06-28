#########################################################################################################################

#TranscriptomeWideAssocaitionStudy_GOEnrich.R
#script to conduct TWAS for metabolites and seed phenotypes for the diversity panel and GO enrichment

#########################################################################################################################

#### 1. Load libraries

library(rrBLUP)
library(statgenGWAS)
library(topGO)


##### 2. Load data and genotypes
diversity_panel_drBLUPs <- read.csv("./DiversityPanel_drBLUPs.csv")

diversity_panel_gbs <- read.csv("./DiversityPanel_GBS.csv")
rownames(diversity_panel_gbs) <- diversity_panel_gbs$X
diversity_panel_gbs <- diversity_panel_gbs[,-c(1)]

diversity_panel_geneExpr <- read.csv("./DiversityPanel_geneExprPEERresid.csv")
rownames(diversity_panel_geneExpr) <- diversity_panel_geneExpr$X
diversity_panel_geneExpr <- diversity_panel_geneExpr[,-c(1)]

###### 3. Calculate relationship matrices
divpanel_amat <- rrBLUP::A.mat(as.matrix(diversity_panel_gbs), min.MAF = 0.05)


###### 4. Conduct TWAS 
#follow #https://cran.r-project.org/web/packages/statgenGWAS/vignettes/GWAS.html

#4a 'marker matrix' (=GE) names of the markers in its column names and the genotypes in its row names. 
#diversity_panel_geneExpr

#4b marker map placeholder for TWAS ; chr pos are columns, 'makers' are row names
marker_map_forSGGWAS <- matrix(nrow=ncol(diversity_panel_geneExpr), ncol = 2)
rownames(marker_map_forSGGWAS) <- colnames(diversity_panel_geneExpr)
colnames(marker_map_forSGGWAS) <- c("chr", "pos")
marker_map_forSGGWAS[,1] <- c(rep(1:4, (nrow(marker_map_forSGGWAS)/4)))
marker_map_forSGGWAS[,2] <- seq(1:nrow(marker_map_forSGGWAS))
marker_map_forSGGWAS <- as.data.frame(marker_map_forSGGWAS)

#4c kinship 
divpanel_amat_for_gencorr <- divpanel_amat[which(rownames(divpanel_amat) %in% rownames(diversity_panel_geneExpr)), which(rownames(divpanel_amat) %in% rownames(diversity_panel_geneExpr))]
dim(divpanel_amat_for_gencorr) # 306 x 306

#4d - pheno first column of all elements of pheno should be genotype and all the other columns should represent different traits. 
pheno_forSGGWAS <- diversity_panel_drBLUPs[which(diversity_panel_drBLUPs$GID_1 %in% rownames(diversity_panel_geneExpr)), c(1,14:ncol(diversity_panel_drBLUPs))] #remove PCs, DTH, line name, and match to samples that have GBS data
colnames(pheno_forSGGWAS)[1] <- "genotype"

#4e covar This data.frame has genotypes in its row names and the covariates in the column names
cols_to_keep_pc <- c("GID_1", "PC1", "PC2","PC3","PC4","PC5")
covar_forSGGWAS <- diversity_panel_drBLUPs[which(diversity_panel_drBLUPs$GID_1 %in% rownames(diversity_panel_geneExpr)),which(colnames(diversity_panel_drBLUPs) %in% cols_to_keep_pc)]
colnames(covar_forSGGWAS)[1] <- "genotype"

dim(diversity_panel_geneExpr)
dim(marker_map_forSGGWAS)
dim(divpanel_amat_for_gencorr)
dim(pheno_forSGGWAS)
dim(covar_forSGGWAS)

# 4f conduct TWAS 
twas_res_diversity_panel <- createGData(geno = diversity_panel_geneExpr, 
                                        map = marker_map_forSGGWAS, 
                                        kin = divpanel_amat_for_gencorr, 
                                        pheno = pheno_forSGGWAS, 
                                        covar = covar_forSGGWAS)

TWAS_results_divpanel <- runSingleTraitGwas(gData = twas_res_diversity_panel)
summary(TWAS_results_divpanel)
#saveRDS(TWAS_results_divpanel, "./TWAS_results_divpanel.RDS")

# 4g apply fdr correction
TWAS_res_with_pdfr <- as.data.frame(matrix(ncol=11))
colnames(TWAS_res_with_pdfr) <- c(colnames(TWAS_results_divpanel$GWAResult$pheno_forSGGWAS), "pfdr")

traits <- unique(TWAS_results_divpanel$GWAResult$pheno_forSGGWAS$trait)
for (i in 1:length(traits)) {
  temp <- TWAS_results_divpanel$GWAResult$pheno_forSGGWAS[TWAS_results_divpanel$GWAResult$pheno_forSGGWAS$trait == traits[i],]
  temp$pfdr <- p.adjust(temp$pValue, method = "fdr")
  TWAS_res_with_pdfr <- rbind(TWAS_res_with_pdfr,temp)
}

TWAS_res_with_pdfr <- TWAS_res_with_pdfr[-c(1),]
dim(TWAS_res_with_pdfr)
#write.csv(TWAS_res_with_pdfr, "./TWAS_res_with_pdfr.csv", row.names = F)

TWAS_res_with_pdfr_05 <- TWAS_res_with_pdfr[TWAS_res_with_pdfr$pfdr < 0.05,]
#write.csv(TWAS_res_with_pdfr[TWAS_res_with_pdfr$pfdr < 0.25,], "./TWASpdfr_p25.csv", row.names = F)


#4g read in gene annotations
#gene_annot <- read.csv("./GeneAnnotations.csv")
#str(gene_annot)
#subet by transcripts of interest



###### 5. Conduct GO enrichment
library(topGO)
library(ALL)
data(ALL)
data(geneList)

#read in RDS of list with all genes and their associated GO terms
avenasat_geneID2GO <- readRDS("./avenasat_geneID2GO.rds")
str(avenasat_geneID2GO)

traits <- unique(TWAS_res_with_pdfr$trait)
for (i in 1:length(traits)) {
  go_temp  <- data.frame(TWAS_res_with_pdfr[TWAS_res_with_pdfr$trait == traits[i], c(2,11)], stringsAsFactors = F)
  colnames(go_temp) <- c("geneID", "padj"); rownames(go_temp) <- NULL; go_temp$padj <- as.numeric(go_temp$padj)
  go_temp_vector <- c(go_temp$padj); names(go_temp_vector) <- go_temp$geneID
  
  GOdata_BP <- new("topGOdata", ontology = "BP", 
                   allGenes = go_temp_vector, 
                   geneSel = topDiffGenes,
                   annot = annFUN.gene2GO, 
                   gene2GO = avenasat_geneID2GO)
  
  resultKS5 <- runTest(GOdata_BP, algorithm = "weight01", statistic = "ks")
  allRes_BP <- GenTable(GOdata_BP,  weKS = resultKS5, orderBy = "weKS",  topNodes = 4710)
  allRes_BP$padj <- p.adjust(allRes_BP$weKS, "fdr")
  
  outputfilename <- paste0("./GO_BP_", traits[i], ".csv")
  #write.csv(allRes_BP, outputfilename, row.names = F)
}








