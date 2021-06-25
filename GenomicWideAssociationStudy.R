#########################################################################################################################

#GenomeWideAssocaitionStudy.R
#script to conduct GWAS for metabolites and seed phenotypes for all germplasm panels

#########################################################################################################################

#### 1. Load libraries

library(rrBLUP)
library(statgenGWAS)


##### 2. Load data and genotypes

diversity_panel_drBLUPs <- read.csv("./DiversityPanel_drBLUPs.csv")
elite_panel_drBLUPs <- read.csv("./ElitePanel_drBLUPs.csv")

diversity_panel_gbs <- read.csv("./DiversityPanel_GBS.csv")
rownames(diversity_panel_gbs) <- diversity_panel_gbs$X
diversity_panel_gbs <- diversity_panel_gbs[,-c(1)]

elite_panel_gbs <- read.csv("./ElitePanel_GBS.csv")
rownames(elite_panel_gbs) <- elite_panel_gbs$X
elite_panel_gbs <- elite_panel_gbs[,-c(1)]

###### 3. Calculate relationship matrices

# 3a. diveristy panel
divpanel_amat <- rrBLUP::A.mat(as.matrix(diversity_panel_gbs), min.MAF = 0.05)
divpanel_amat[1:10,1:10]
str(divpanel_amat)

# 3b. elite panel by location
elite_panel_drBLUPs_MN <- elite_panel_drBLUPs[elite_panel_drBLUPs$env == "MN",]
elite_panel_drBLUPs_SD <- elite_panel_drBLUPs[elite_panel_drBLUPs$env == "SD",]
elite_panel_drBLUPs_WI <- elite_panel_drBLUPs[elite_panel_drBLUPs$env == "WI",]


elite_mn_amat <- rrBLUP::A.mat(as.matrix(elite_panel_gbs[which(rownames(elite_panel_gbs) %in% elite_panel_drBLUPs_MN$line),]), min.MAF = 0.05)
elite_mn_amat[1:10,1:10]

elite_sd_amat <- rrBLUP::A.mat(as.matrix(elite_panel_gbs[which(rownames(elite_panel_gbs) %in% elite_panel_drBLUPs_SD$line),]), min.MAF = 0.05)
elite_sd_amat[1:10,1:10]

elite_wi_amat <- rrBLUP::A.mat(as.matrix(elite_panel_gbs[which(rownames(elite_panel_gbs) %in% elite_panel_drBLUPs_WI$line),]), min.MAF = 0.05)
elite_wi_amat[1:10,1:10]


###### 4. Conduct GWAS 
#follow #https://cran.r-project.org/web/packages/statgenGWAS/vignettes/GWAS.html

#4a GWAS for diversity panel

# 4a.1 marker map = chromosome, postion are columns, makers are row names
marker_map_forSGGWAS <- read.csv("./marker_genome_position.csv")
rownames(marker_map_forSGGWAS) <- marker_map_forSGGWAS$Marker
marker_map_forSGGWAS <- marker_map_forSGGWAS[,c(3,2)]
colnames(marker_map_forSGGWAS) <- c("chr", "pos")

#4a.2 marker matrix names of the markers in its column names and the genotypes in its row names. 
marker_matrix_forSGGWAS <-  as.data.frame(diversity_panel_gbs)

#4a.3 - pheno first column of all elements of pheno should be genotype and all the other columns should represent different traits. 
pheno_forSGGWAS <- diversity_panel_drBLUPs[which(diversity_panel_drBLUPs$GID_1 %in% rownames(marker_matrix_forSGGWAS)), c(1,14:ncol(diversity_panel_drBLUPs))] #remove PCs, DTH, line name, and match to samples that have GBS data
colnames(pheno_forSGGWAS)[1] <- "genotype"

#4a.4 kinship - match to invididuals with genotypes, phenotyptes
divpanel_amat_for_gencorr <- divpanel_amat[which(rownames(divpanel_amat) %in% diversity_panel_drBLUPs$GID_1), which(rownames(divpanel_amat) %in% diversity_panel_drBLUPs$GID_1)]
dim(divpanel_amat_for_gencorr) # 343 x 343

#4a.5 covar This data.frame has genotypes in its row names and the covariates in the column names
cols_to_keep_pc <- c("GID_1", "PC1", "PC2","PC3","PC4","PC5")
covar_forSGGWAS <- diversity_panel_drBLUPs[,which(colnames(diversity_panel_drBLUPs) %in% cols_to_keep_pc)]
colnames(covar_forSGGWAS)[1] <- "genotype"

# 4a.6 run GWAS
gwas_res_diversity_panel <- createGData(geno = marker_matrix_forSGGWAS, 
                              map = marker_map_forSGGWAS, 
                              kin = divpanel_amat_for_gencorr, 
                              pheno = pheno_forSGGWAS, 
                              covar = covar_forSGGWAS)

GWAS_results_divpanel <- runSingleTraitGwas(gData = gwas_res_diversity_panel)
summary(GWAS_results_divpanel)



#4b GWAS for elite panel

#4b.1 marker matrix names of the markers in its column names and the genotypes in its row names. 
marker_matrix_forSGGWAS_elite <-  as.data.frame(elite_panel_gbs)
dim(marker_matrix_forSGGWAS_elite)

# 4b.2 marker map = chromosome, postion are columns, makers are row names
marker_map_forSGGWAS_elite <- marker_map_forSGGWAS[which(rownames(marker_map_forSGGWAS) %in% colnames(marker_matrix_forSGGWAS_elite)),]
dim(marker_map_forSGGWAS_elite)

#4a.3 - pheno first column of all elements of pheno should be genotype and all the other columns should represent different traits. 
pheno_forSGGWAS_MN <- elite_panel_drBLUPs_MN[which(elite_panel_drBLUPs_MN$line %in% rownames(marker_matrix_forSGGWAS_elite)), c(2,15:ncol(elite_panel_drBLUPs_MN))] #remove PCs, DTH, line name, and match to samples that have GBS data
pheno_forSGGWAS_SD <- elite_panel_drBLUPs_SD[which(elite_panel_drBLUPs_SD$line %in% rownames(marker_matrix_forSGGWAS_elite)), c(2,15:ncol(elite_panel_drBLUPs_SD))] #remove PCs, DTH, line name, and match to samples that have GBS data
pheno_forSGGWAS_WI <- elite_panel_drBLUPs_WI[which(elite_panel_drBLUPs_WI$line %in% rownames(marker_matrix_forSGGWAS_elite)), c(2,15:ncol(elite_panel_drBLUPs_WI))] #remove PCs, DTH, line name, and match to samples that have GBS data
colnames(pheno_forSGGWAS_MN)[1] <- "genotype"; colnames(pheno_forSGGWAS_SD)[1] <- "genotype"; colnames(pheno_forSGGWAS_WI)[1] <- "genotype"

#4a.4 kinship - match to invididuals with genotypes, phenotyptes
#use kinship calculated above

#4a.5 covar This data.frame has genotypes in its row names and the covariates in the column names
cols_to_keep_pc <- c("line", "PC1", "PC2","PC3","PC4")
covar_forSGGWAS_MN <- elite_panel_drBLUPs_MN[which(elite_panel_drBLUPs_MN$line %in% rownames(marker_matrix_forSGGWAS_elite)),which(colnames(elite_panel_drBLUPs_MN) %in% cols_to_keep_pc)]
covar_forSGGWAS_SD <- elite_panel_drBLUPs_SD[which(elite_panel_drBLUPs_SD$line %in% rownames(marker_matrix_forSGGWAS_elite)),which(colnames(elite_panel_drBLUPs_SD) %in% cols_to_keep_pc)]
covar_forSGGWAS_WI <- elite_panel_drBLUPs_WI[which(elite_panel_drBLUPs_WI$line %in% rownames(marker_matrix_forSGGWAS_elite)),which(colnames(elite_panel_drBLUPs_WI) %in% cols_to_keep_pc)]
colnames(covar_forSGGWAS_MN)[1] <- "genotype";colnames(covar_forSGGWAS_SD)[1] <- "genotype";colnames(covar_forSGGWAS_WI)[1] <- "genotype"


# 4a.6 run GWAS
gwas_res_elite_panel_MN <- createGData(geno = marker_matrix_forSGGWAS_elite, 
                                        map = marker_map_forSGGWAS_elite, 
                                        kin = elite_mn_amat, 
                                        pheno = pheno_forSGGWAS_MN, 
                                        covar = covar_forSGGWAS_MN)

GWAS_results_eliteMN <- runSingleTraitGwas(gData = gwas_res_elite_panel_MN)
summary(GWAS_results_eliteMN)


gwas_res_elite_panel_SD <- createGData(geno = marker_matrix_forSGGWAS_elite, 
                                       map = marker_map_forSGGWAS_elite, 
                                       kin = elite_sd_amat, 
                                       pheno = pheno_forSGGWAS_SD, 
                                       covar = covar_forSGGWAS_SD)

GWAS_results_eliteSD <- runSingleTraitGwas(gData = gwas_res_elite_panel_SD)
summary(GWAS_results_eliteSD)


gwas_res_elite_panel_WI <- createGData(geno = marker_matrix_forSGGWAS_elite, 
                                       map = marker_map_forSGGWAS_elite, 
                                       kin = elite_wi_amat, 
                                       pheno = pheno_forSGGWAS_WI, 
                                       covar = covar_forSGGWAS_WI)

GWAS_results_eliteWI <- runSingleTraitGwas(gData = gwas_res_elite_panel_WI)
summary(GWAS_results_eliteWI)



all_GWAS_results <- list(GWAS_results_divpanel = GWAS_results_divpanel,
                         GWAS_results_eliteMN = GWAS_results_eliteMN,
                         GWAS_results_eliteSD = GWAS_results_eliteSD,
                         GWAS_results_eliteWI = GWAS_results_eliteWI)

#saveRDS(all_GWAS_results, "./all_GWAS_results.RDS")
