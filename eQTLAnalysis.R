#########################################################################################################################

#eQTLAnalysis.R
#script to conduct eQTL analysis 

#########################################################################################################################

#### 1. Load libraries
library(MatrixEQTL)


##### 2. Load data and genotypes
diversity_panel_drBLUPs <- read.csv("./DiversityPanel_drBLUPs.csv")

diversity_panel_gbs <- read.csv("./DiversityPanel_GBS.csv")
rownames(diversity_panel_gbs) <- diversity_panel_gbs$X
diversity_panel_gbs <- diversity_panel_gbs[,-c(1)]

diversity_panel_geneExpr <- read.csv("./DiversityPanel_geneExprPEERresid.csv")
rownames(diversity_panel_geneExpr) <- diversity_panel_geneExpr$X
diversity_panel_geneExpr <- diversity_panel_geneExpr[,-c(1)]

# see TWAS script to generate the TWAS_res_with_pfdr_05 data frame
# this is list of TWAS transcripts of interest (shared with AVN_A, AVN_B)
transcript_list <- intersect(unique(TWAS_res_with_pdfr_05[TWAS_res_with_pdfr_05$trait == "LC.03.1204" , 2]), unique(TWAS_res_with_pdfr_05[TWAS_res_with_pdfr_05$trait == "LC.03.0125" , 2]))
transcript_list <- transcript_list[-c(1)] #remove NA




#### 3. Prepare files for eQTL analysis
# follow http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html

#gbs snps
snp_for_matrixeQTL <- diversity_panel_gbs[which(rownames(diversity_panel_gbs) %in% rownames(diversity_panel_geneExpr)),]
snp_for_matrixeQTL <- as.data.frame(t(snp_for_matrixeQTL))
snp_for_matrixeQTL <- as.data.frame(cbind(rownames(snp_for_matrixeQTL), snp_for_matrixeQTL))
colnames(snp_for_matrixeQTL)[1] <- "id"
str(snp_for_matrixeQTL[,1:10])
dim(snp_for_matrixeQTL) #54284 x  307
write.table(snp_for_matrixeQTL, "./snp_for_matrixeQTL.txt", row.names = F, sep = "\t")

#gene expression
ge_for_matrixeQTL_50rem_peer_twas <- as.data.frame(t(diversity_panel_geneExpr[,which(colnames(diversity_panel_geneExpr) %in% transcript_list)]))
ge_for_matrixeQTL_50rem_peer_twas <- as.data.frame(cbind(rownames(ge_for_matrixeQTL_50rem_peer_twas), ge_for_matrixeQTL_50rem_peer_twas))
colnames(ge_for_matrixeQTL_50rem_peer_twas)[1] <- "id"
str(ge_for_matrixeQTL_50rem_peer_twas) #51 x 307
write.table(ge_for_matrixeQTL_50rem_peer_twas, "./ge_for_matrixeQTL_twas_only.txt", row.names = F, sep = "\t")

# covariates (PCs)
cols_to_keep_pc <- c("GID_1", "PC1", "PC2","PC3","PC4","PC5")
cov_for_matrixeQTL <- diversity_panel_drBLUPs[which(diversity_panel_drBLUPs$GID_1 %in% rownames(diversity_panel_geneExpr)),which(colnames(diversity_panel_drBLUPs) %in% cols_to_keep_pc)]
rownames(cov_for_matrixeQTL) <- cov_for_matrixeQTL$GID_1
cov_for_matrixeQTL <- cov_for_matrixeQTL[,-c(1)]
head(cov_for_matrixeQTL)
cov_for_matrixeQTL <- as.data.frame(t(cov_for_matrixeQTL))
cov_for_matrixeQTL <- as.data.frame(cbind(rownames(cov_for_matrixeQTL), cov_for_matrixeQTL))
colnames(cov_for_matrixeQTL)[1] <- "id"
str(cov_for_matrixeQTL[,1:10])
dim(cov_for_matrixeQTL) #5x307 = 306 inds
write.table(cov_for_matrixeQTL, "./cov_for_matrixeQTL.txt", row.names = F, sep = "\t")


###### 4 run analysis 
## Settings
base.dir = find.package('MatrixEQTL');
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = "./snp_for_matrixeQTL.txt"
covariates_file_name = "./cov_for_matrixeQTL.txt"
ge_file_name= "./ge_for_matrixeQTL_twas_only.txt"

# Output file name
output_file_name= tempfile();

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(ge_file_name);
##

#Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis with stringent pvalue threshold 
pvOutputThreshold3 = 0.99 #chosen for evaluating TWAS results only; if evaluting larger data sets, choose higher cutoff

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold3,
  useModel = useModel,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

eqtl_res_twas_only <- as.data.frame(me$all$eqtls)
dim(eqtl_res_twas_only)
head(eqtl_res_twas_only)
length(unique(eqtl_res_twas_only$gene)) #51
length(unique(eqtl_res_twas_only$snps)) #54284


##### 5. calculate FDR by gene, and keep eQTL that pass a pfdr < 0.25 threshold

eqtl_res_twas_only_temp <- eqtl_res_twas_only[eqtl_res_twas_only$gene == transcript_list[1],]
head(eqtl_res_twas_only_temp)
hist(eqtl_res_twas_only_temp$pvalue)
eqtl_res_twas_only_temp$pfdr <- p.adjust(eqtl_res_twas_only_temp$pvalue, "fdr")
plot(eqtl_res_twas_only_temp$pfdr, eqtl_res_twas_only_temp$FDR)
dim(eqtl_res_twas_only_temp[eqtl_res_twas_only_temp$pfdr < 0.25,])
head(eqtl_res_twas_only_temp)

eqt_res_twas_fdr_p25 <- eqtl_res_twas_only_temp[eqtl_res_twas_only_temp$pfdr < 0.25,]

for (i in 2:length(transcript_list)) {
  eqtl_res_twas_only_temp <- eqtl_res_twas_only[eqtl_res_twas_only$gene == transcript_list[i],]
  eqtl_res_twas_only_temp$pfdr <- p.adjust(eqtl_res_twas_only_temp$pvalue, "fdr")
  eqt_res_twas_fdr_p25 <- rbind(eqt_res_twas_fdr_p25, eqtl_res_twas_only_temp[eqtl_res_twas_only_temp$pfdr < 0.25,])
}

#dim(eqt_res_twas_fdr_p25)
#write.csv(eqt_res_twas_fdr_p25, "./eqt_res_twas_fdr_p25.csv", row.names = F)
