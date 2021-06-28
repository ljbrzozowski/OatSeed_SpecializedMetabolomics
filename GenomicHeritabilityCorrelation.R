#########################################################################################################################

#GenomicHeritabilityCorrelation.R
#script to calculate genomic heritability and genomic correlation of metabolites and seed phenotypes for all germplasm panels

#########################################################################################################################

#### 1. Load libraries

library(rrBLUP)
library(sommer)


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


###### 4. Calculate heritability

# 4a. diversity panel

# subset data by individuals that have both genotypes and phenotypes
div_panel_drblups_final_for_gencor <- diversity_panel_drBLUPs[which(diversity_panel_drBLUPs$GID_1 %in% rownames(divpanel_amat)) ,]
str(div_panel_drblups_final_for_gencor) #343 x 23

#calculate genomic heritability
ans_avnA <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="LC.03.1204", K=divpanel_amat)
ans_avnA$Vg / (ans_avnA$Vg + ans_avnA$Ve) 
ans_avnB <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="LC.03.0125", K=divpanel_amat)
ans_avnB$Vg / (ans_avnB$Vg + ans_avnB$Ve) 
ans_aec1 <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="LC.03.0066", K=divpanel_amat)
ans_aec1$Vg / (ans_aec1$Vg + ans_aec1$Ve) 
ans_aec2 <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="LC.03.0609", K=divpanel_amat)
ans_aec2$Vg / (ans_aec2$Vg + ans_aec2$Ve) 
ans_aosA <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="LC.03.0001", K=divpanel_amat)
ans_aosA$Vg / (ans_aosA$Vg + ans_aosA$Ve) 
ans_aosdA <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="LC.03.0052", K=divpanel_amat)
ans_aosdA$Vg / (ans_aosdA$Vg + ans_aosdA$Ve) 
ans_aosB <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="LC.03.0061", K=divpanel_amat)
ans_aosB$Vg / (ans_aosB$Vg + ans_aosB$Ve) 
ans_sv <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="seedVol", K=divpanel_amat)
ans_sv$Vg / (ans_sv$Vg + ans_sv$Ve) 
ans_ssa <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="seedSA", K=divpanel_amat)
ans_ssa$Vg / (ans_ssa$Vg + ans_ssa$Ve) 
ans_ssavol <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="SA_vol_ratio", K=divpanel_amat)
ans_ssavol$Vg / (ans_ssavol$Vg + ans_ssavol$Ve) 
ans_sl <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="Seedlength", K=divpanel_amat)
ans_sl$Vg / (ans_sl$Vg + ans_sl$Ve) 
ans_sw <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="SeedWidth", K=divpanel_amat)
ans_sw$Vg / (ans_sw$Vg + ans_sw$Ve) 
ans_sh <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="Seedheight", K=divpanel_amat)
ans_sh$Vg / (ans_sh$Vg + ans_sh$Ve) 
ans_hkw <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="hundred_kernel_weight", K=divpanel_amat)
ans_hkw$Vg / (ans_hkw$Vg + ans_hkw$Ve) 
ans_hhw <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="hundred_hull_weight", K=divpanel_amat)
ans_hhw$Vg / (ans_hhw$Vg + ans_hhw$Ve) 
ans_gp <- rrBLUP::kin.blup(data=div_panel_drblups_final_for_gencor, geno="GID_1", pheno="groat_pct", K=divpanel_amat)
ans_gp$Vg / (ans_gp$Vg + ans_gp$Ve) 



# 4b. elite panel
# subset data by individuals that have both genotypes and phenotypes
elite_mn_panel_drblups_final_for_gencor <- elite_panel_drBLUPs_MN[which(elite_panel_drBLUPs_MN$line %in% rownames(elite_mn_amat)) ,]
str(elite_mn_panel_drblups_final_for_gencor) #192, 30
dim(elite_mn_amat)#192, 192

elite_sd_panel_drblups_final_for_gencor <- elite_panel_drBLUPs_SD[which(elite_panel_drBLUPs_SD$line %in% rownames(elite_sd_amat)) ,]
str(elite_sd_panel_drblups_final_for_gencor) #191, 13
dim(elite_sd_amat)#191, 191

elite_wi_panel_drblups_final_for_gencor <- elite_panel_drBLUPs_WI[which(elite_panel_drBLUPs_WI$line %in% rownames(elite_wi_amat)) ,]
str(elite_wi_panel_drblups_final_for_gencor) #191, 13
dim(elite_wi_amat) #191, 191

ans_avnA <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0325", K=elite_mn_amat)
ans_avnA$Vg / (ans_avnA$Vg + ans_avnA$Ve) 
ans_avnB <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0242", K=elite_mn_amat)
ans_avnB$Vg / (ans_avnB$Vg + ans_avnB$Ve) 
ans_aec1 <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0099", K=elite_mn_amat)
ans_aec1$Vg / (ans_aec1$Vg + ans_aec1$Ve) 
ans_aec2 <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0551", K=elite_mn_amat)
ans_aec2$Vg / (ans_aec2$Vg + ans_aec2$Ve) 
ans_aosA <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0012", K=elite_mn_amat)
ans_aosA$Vg / (ans_aosA$Vg + ans_aosA$Ve) 
ans_aosdA <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0076", K=elite_mn_amat)
ans_aosdA$Vg / (ans_aosdA$Vg + ans_aosdA$Ve) 
ans_aosB <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0441", K=elite_mn_amat)
ans_aosB$Vg / (ans_aosB$Vg + ans_aosB$Ve) 
ans_sv <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="seedVol", K=elite_mn_amat)
ans_sv$Vg / (ans_sv$Vg + ans_sv$Ve) 
ans_ssa <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="seedSA", K=elite_mn_amat)
ans_ssa$Vg / (ans_ssa$Vg + ans_ssa$Ve) 
ans_ssavol <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="SA_vol_ratio", K=elite_mn_amat)
ans_ssavol$Vg / (ans_ssavol$Vg + ans_ssavol$Ve) 
ans_sl <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="Seedlength", K=elite_mn_amat)
ans_sl$Vg / (ans_sl$Vg + ans_sl$Ve) 
ans_sw <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="Seedwidth", K=elite_mn_amat)
ans_sw$Vg / (ans_sw$Vg + ans_sw$Ve) 
ans_sh <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="Seedheight", K=elite_mn_amat)
ans_sh$Vg / (ans_sh$Vg + ans_sh$Ve) 
ans_hkw <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="hundred_kernel_weight", K=elite_mn_amat)
ans_hkw$Vg / (ans_hkw$Vg + ans_hkw$Ve) 
ans_hhw <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="hundred_hull_weight", K=elite_mn_amat)
ans_hhw$Vg / (ans_hhw$Vg + ans_hhw$Ve)  
ans_gp <- rrBLUP::kin.blup(data=elite_mn_panel_drblups_final_for_gencor, geno="line", pheno="groat_pct", K=elite_mn_amat)
ans_gp$Vg / (ans_gp$Vg + ans_gp$Ve) 

ans_avnA <- rrBLUP::kin.blup(data=elite_sd_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0325", K=elite_sd_amat)
ans_avnA$Vg / (ans_avnA$Vg + ans_avnA$Ve) 
ans_avnB <- rrBLUP::kin.blup(data=elite_sd_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0242", K=elite_sd_amat)
ans_avnB$Vg / (ans_avnB$Vg + ans_avnB$Ve) 
ans_aec1 <- rrBLUP::kin.blup(data=elite_sd_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0099", K=elite_sd_amat)
ans_aec1$Vg / (ans_aec1$Vg + ans_aec1$Ve) 
ans_aec2 <- rrBLUP::kin.blup(data=elite_sd_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0551", K=elite_sd_amat)
ans_aec2$Vg / (ans_aec2$Vg + ans_aec2$Ve) 
ans_aosA <- rrBLUP::kin.blup(data=elite_sd_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0012", K=elite_sd_amat)
ans_aosA$Vg / (ans_aosA$Vg + ans_aosA$Ve) 
ans_aosdA <- rrBLUP::kin.blup(data=elite_sd_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0076", K=elite_sd_amat)
ans_aosdA$Vg / (ans_aosdA$Vg + ans_aosdA$Ve) 
ans_aosB <- rrBLUP::kin.blup(data=elite_sd_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0441", K=elite_sd_amat)
ans_aosB$Vg / (ans_aosB$Vg + ans_aosB$Ve) 
ans_sh <- rrBLUP::kin.blup(data=elite_sd_panel_drblups_final_for_gencor, geno="line", pheno="Seedheight", K=elite_sd_amat)
ans_sh$Vg / (ans_sh$Vg + ans_sh$Ve) 
ans_hkw <- rrBLUP::kin.blup(data=elite_sd_panel_drblups_final_for_gencor, geno="line", pheno="hundred_kernel_weight", K=elite_sd_amat)
ans_hkw$Vg / (ans_hkw$Vg + ans_hkw$Ve) 
ans_hhw <- rrBLUP::kin.blup(data=elite_sd_panel_drblups_final_for_gencor, geno="line", pheno="hundred_hull_weight", K=elite_sd_amat)
ans_hhw$Vg / (ans_hhw$Vg + ans_hhw$Ve)  
ans_gp <- rrBLUP::kin.blup(data=elite_sd_panel_drblups_final_for_gencor, geno="line", pheno="groat_pct", K=elite_sd_amat)
ans_gp$Vg / (ans_gp$Vg + ans_gp$Ve) 

ans_avnA <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0325", K=elite_wi_amat)
ans_avnA$Vg / (ans_avnA$Vg + ans_avnA$Ve) 
ans_avnB <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0242", K=elite_wi_amat)
ans_avnB$Vg / (ans_avnB$Vg + ans_avnB$Ve) 
ans_aec1 <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0099", K=elite_wi_amat)
ans_aec1$Vg / (ans_aec1$Vg + ans_aec1$Ve) 
ans_aec2 <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="LC.02.0551", K=elite_wi_amat)
ans_aec2$Vg / (ans_aec2$Vg + ans_aec2$Ve) 
ans_sv <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="seedVol", K=elite_wi_amat)
ans_sv$Vg / (ans_sv$Vg + ans_sv$Ve) 
ans_ssa <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="seedSA", K=elite_wi_amat)
ans_ssa$Vg / (ans_ssa$Vg + ans_ssa$Ve) 
ans_ssavol <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="SA_vol_ratio", K=elite_wi_amat)
ans_ssavol$Vg / (ans_ssavol$Vg + ans_ssavol$Ve) 
ans_sl <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="Seedlength", K=elite_wi_amat)
ans_sl$Vg / (ans_sl$Vg + ans_sl$Ve) 
ans_sw <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="Seedwidth", K=elite_wi_amat)
ans_sw$Vg / (ans_sw$Vg + ans_sw$Ve) 
ans_sh <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="Seedheight", K=elite_wi_amat)
ans_sh$Vg / (ans_sh$Vg + ans_sh$Ve) 
ans_hkw <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="hundred_kernel_weight", K=elite_wi_amat)
ans_hkw$Vg / (ans_hkw$Vg + ans_hkw$Ve) 
ans_hhw <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="hundred_hull_weight", K=elite_wi_amat)
ans_hhw$Vg / (ans_hhw$Vg + ans_hhw$Ve)  
ans_gp <- rrBLUP::kin.blup(data=elite_wi_panel_drblups_final_for_gencor, geno="line", pheno="groat_pct", K=elite_wi_amat)
ans_gp$Vg / (ans_gp$Vg + ans_gp$Ve) 


###### 5. Calculate genetic correlations
# run in separate analyses so that converge
# only used traits with h2 > 0.05

divpanel_amat_for_gencorr <- divpanel_amat[which(rownames(divpanel_amat) %in% diversity_panel_drBLUPs$GID_1), which(rownames(divpanel_amat) %in% diversity_panel_drBLUPs$GID_1)]
dim(divpanel_amat_for_gencorr) # 343 x 343
div_panel_drblups_final_for_gencor$GID_1 <- as.character(div_panel_drblups_final_for_gencor$GID_1)

# 5a. diversity panel
ans.p2 <- mmer(cbind(LC.03.0001,LC.03.0052,LC.03.0061,LC.03.0066,LC.03.0609,LC.03.1204,LC.03.0125)~1,
               random=~ vs(GID_1, Gu=divpanel_amat_for_gencorr, Gtc=unsm(7)), rcov=~ vs(units, Gtc=unsm(7)), data=div_panel_drblups_final_for_gencor, tolparinv = 1e-2)
gencor_div2 <- cov2cor(ans.p2$sigma$`u:GID_1`)
gencor_div2 

ans.p2.seed <- mmer(cbind(LC.03.0001,LC.03.0052,LC.03.0061,LC.03.0066,LC.03.0609,LC.03.1204,LC.03.0125,seedVol,seedSA)~1,
                    random=~ vs(GID_1, Gu=divpanel_amat_for_gencorr, Gtc=unsm(9)), rcov=~ vs(units, Gtc=unsm(9)), data=div_panel_drblups_final_for_gencor, tolparinv = 1e-2)
gencor_div2.seed <- cov2cor(ans.p2.seed$sigma$`u:GID_1`)
gencor_div2.seed 

ans.p2.seedw <- mmer(cbind(LC.03.0001,LC.03.0052,LC.03.0061,LC.03.0066,LC.03.0609,LC.03.1204,LC.03.0125,hundred_kernel_weight,hundred_hull_weight,groat_pct)~1,
                     random=~ vs(GID_1, Gu=divpanel_amat_for_gencorr, Gtc=unsm(10)), rcov=~ vs(units, Gtc=unsm(10)), data=div_panel_drblups_final_for_gencor, tolparinv = 1e-2)
gencor_div2.seedw <- cov2cor(ans.p2.seedw$sigma$`u:GID_1`)
gencor_div2.seedw 


# 5b. Elite panel by location
elite_mn_panel_drblups_final_for_gencor$line <- as.character(elite_mn_panel_drblups_final_for_gencor$line)
elite_sd_panel_drblups_final_for_gencor$line <- as.character(elite_sd_panel_drblups_final_for_gencor$line)
elite_wi_panel_drblups_final_for_gencor$line <- as.character(elite_wi_panel_drblups_final_for_gencor$line)

ans_mn2_h2_05_seed <- mmer(cbind(LC.02.0325, LC.02.0242, LC.02.0099, LC.02.0551, LC.02.0012, LC.02.0076,seedVol,seedSA)~1,
                           random=~ vs(line, Gu=elite_mn_amat, Gtc=unsm(8)), rcov=~ vs(units, Gtc=unsm(8)), data=elite_mn_panel_drblups_final_for_gencor, tolparinv = 1e-2)
gencor_elitemn2_h205_seed <- cov2cor(ans_mn2_h2_05_seed$sigma$`u:line`)
gencor_elitemn2_h205_seed

ans_mn2 <- mmer(cbind(LC.02.0325, LC.02.0242,  LC.02.0099, LC.02.0551, LC.02.0012, LC.02.0076, hundred_kernel_weight, groat_pct)~1,
                random=~ vs(line, Gu=elite_mn_amat, Gtc=unsm(8)), rcov=~ vs(units, Gtc=unsm(8)), data=elite_mn_panel_drblups_final_for_gencor, tolparinv = 1e-2)
gencor_elitemn2 <- cov2cor(ans_mn2$sigma$`u:line`)
gencor_elitemn2

ans_sd2 <- mmer(cbind(LC.02.0325, LC.02.0242, LC.02.0099, LC.02.0551, LC.02.0012, hundred_kernel_weight, groat_pct)~1,
                random=~ vs(line, Gu=elite_sd_amat, Gtc=unsm(7)), rcov=~ vs(units, Gtc=unsm(7)), data=elite_sd_panel_drblups_final_for_gencor, tolparinv = 1e-2)
gencor_elitesd2 <- cov2cor(ans_sd2$sigma$`u:line`)
gencor_elitesd2

ans_wi2_seed <- mmer(cbind(LC.02.0325, LC.02.0242, LC.02.0099, LC.02.0551, seedVol, seedSA)~1,
                     random=~ vs(line, Gu=elite_wi_amat, Gtc=unsm(6)), rcov=~ vs(units, Gtc=unsm(6)), data=elite_wi_panel_drblups_final_for_gencor, tolparinv = 1e-2)
gencor_elitewi2_seed <- cov2cor(ans_wi2_seed$sigma$`u:line`)
gencor_elitewi2_seed

ans_wi2 <- mmer(cbind(LC.02.0325, LC.02.0242, LC.02.0099, LC.02.0551, hundred_kernel_weight, groat_pct)~1,
                random=~ vs(line, Gu=elite_wi_amat, Gtc=unsm(6)), rcov=~ vs(units, Gtc=unsm(6)), data=elite_wi_panel_drblups_final_for_gencor, tolparinv = 1e-2)
gencor_elitewi2 <- cov2cor(ans_wi2$sigma$`u:line`)
gencor_elitewi2







