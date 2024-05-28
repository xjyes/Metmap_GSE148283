setwd("/Users/xujingyi/Documents/szbl/metastatic cell genomic/code/met map")
rm(list = ls())

library(readxl)
library(limma)
library(edgeR)


organ <- c("brain","bone","kidney","liver","lung","HCC1954_BREAST","JIMT1_BREAST","HCC1806_BREAST")

# Separate gene sets
rm(list = ls())
load("data148283.Rdata")
load("information.Rdata")
brain_g1 <- data148283[,information$Sample_id %in% c("GSM4459610","GSM4459617","GSM4459622","GSM4459632","GSM4459614","GSM4459627") ]
bone_g1 <- data148283[,information[information$Sample_characteristics == "group 1 (g1)" & information$tissue_or_cell_type == "bone",]$Sample_id]
kidney_g1 <- data148283[,information[information$Sample_characteristics == "group 1 (g1)" & information$tissue_or_cell_type == "kidney",]$Sample_id]
liver_g1 <- data148283[,information[information$Sample_characteristics == "group 1 (g1)" & information$tissue_or_cell_type == "liver",]$Sample_id]
lung_g1 <- data148283[,information[information$Sample_characteristics == "group 1 (g1)" & information$tissue_or_cell_type == "lung",]$Sample_id]
JIMT_g1 <- data148283[,information[information$Sample_characteristics == "group 1 (g1)" & information$tissue_or_cell_type == "JIMT1_BREAST",]$Sample_id]
HCC_g1 <- data148283[,information[information$Sample_characteristics == "group 1 (g1)" & information$tissue_or_cell_type == "HCC1806_BREAST",]$Sample_id]
brainJIMT_g1 <- data148283[,information$Sample_id %in% c("GSM4459610","GSM4459617","GSM4459622","GSM4459632") ]
brainHCC_g1 <- data148283[, information$Sample_id %in% c("GSM4459614","GSM4459627")]
mammary_g1 <- data148283[, information$Sample_id %in% c("GSM4459636","GSM4459637","GSM4459638","GSM4459639")]

# Get group information
gl_bb_g1 <- c(rep("control",ncol(brain_g1)),rep("treat",ncol(bone_g1)))
gl_bk_g1 <- c(rep("control",ncol(brain_g1)),rep("treat",ncol(kidney_g1)))
gl_blv_g1 <- c(rep("control",ncol(brain_g1)),rep("treat",ncol(liver_g1)))
gl_blg_g1 <- c(rep("control",ncol(brain_g1)),rep("treat",ncol(lung_g1)))
gl_bJIMT_g1 <- c(rep("control",ncol(brainJIMT_g1)),rep("treat",ncol(JIMT_g1)))
gl_bHCC_g1 <- c(rep("control",ncol(brainHCC_g1)),rep("treat",ncol(HCC_g1)))
gl_bm_g1 <- c(rep("control",ncol(brainJIMT_g1)),rep("treat",ncol(mammary_g1)))

save(gl_bb_g1,gl_bk_g1,gl_blv_g1,gl_blg_g1,gl_bJIMT_g1, gl_bHCC_g1, gl_bm_g1, file = "group_g1.Rdata")

# Matrix preparation
options(stringsAsFactors = F)

exp_bb_g1 <- cbind(brain_g1,bone_g1)
log_bb_g1 <- log2(exp_bb_g1 + 1)
designBB_g1 <- model.matrix(~0+factor(gl_bb_g1))
colnames(designBB_g1) = levels(factor(gl_bb_g1))
rownames(designBB_g1) = colnames(log_bb_g1)
designBB_g1

exp_bk_g1 <- cbind(brain_g1,kidney_g1)
log_bk_g1 <- log2(exp_bk_g1 + 1)
designBK_g1 <- model.matrix(~0+factor(gl_bk_g1))
colnames(designBK_g1) = levels(factor(gl_bk_g1))
rownames(designBK_g1) = colnames(log_bk_g1)
designBK_g1

exp_blv_g1 <- cbind(brain_g1,liver_g1)
log_blv_g1 <- log2(exp_blv_g1 + 1)
designBLV_g1 <- model.matrix(~0+factor(gl_blv_g1))
colnames(designBLV_g1) = levels(factor(gl_blv_g1))
rownames(designBLV_g1) = colnames(log_blv_g1)
designBLV_g1

exp_blg_g1 <- cbind(brain_g1,lung_g1)
log_blg_g1 <- log2(exp_blg_g1 + 1)
designBLG_g1 <- model.matrix(~0+factor(gl_blg_g1))
colnames(designBLG_g1) = levels(factor(gl_blg_g1))
rownames(designBLG_g1) = colnames(log_blg_g1)
designBLG_g1

exp_bJIMT_g1 <- cbind(brainJIMT_g1,JIMT_g1)
log_bJIMT_g1 <- log2(exp_bJIMT_g1 + 1)
designBJIMT_g1 <- model.matrix(~0+factor(gl_bJIMT_g1))
colnames(designBJIMT_g1) = levels(factor(gl_bJIMT_g1))
rownames(designBJIMT_g1) = colnames(log_bJIMT_g1)
designBJIMT_g1

exp_bHCC_g1 <- cbind(brainHCC_g1,HCC_g1)
log_bHCC_g1 <- log2(exp_bHCC_g1 + 1)
designBHCC_g1 <- model.matrix(~0+factor(gl_bHCC_g1))
colnames(designBHCC_g1) = levels(factor(gl_bHCC_g1))
rownames(designBHCC_g1) = colnames(log_bHCC_g1)
designBHCC_g1

exp_bm_g1 <- cbind(brainJIMT_g1,mammary_g1)
log_bm_g1 <- log2(exp_bm_g1 + 1)
designBM_g1 <- model.matrix(~0+factor(gl_bm_g1))
colnames(designBM_g1) = levels(factor(gl_bm_g1))
rownames(designBM_g1) = colnames(log_bm_g1)
designBM_g1


contrast.matrix <- makeContrasts(paste0(c("treat","control"),collapse = "-"),levels = designBM_g1)
contrast.matrix


# Differentially expressed gene analysis with limma
fitBB_g1 <- lmFit(log_bb_g1,designBB_g1)
fitBB_g1.2 <- contrasts.fit(fitBB_g1,contrast.matrix)
fitBB_g1.2 <- eBayes(fitBB_g1.2)
tempOutput = topTable(fitBB_g1.2, coef = 1, n=Inf)
nrDEGBB_g1 = na.omit(tempOutput)
head(nrDEGBB_g1)

fitBK_g1 <- lmFit(log_bk_g1,designBK_g1)
fitBK_g1.2 <- contrasts.fit(fitBK_g1,contrast.matrix)
fitBK_g1.2 <- eBayes(fitBK_g1.2)
tempOutput = topTable(fitBK_g1.2, coef = 1, n=Inf)
nrDEGBK_g1 = na.omit(tempOutput)
head(nrDEGBK_g1)

fitBLV_g1 <- lmFit(log_blv_g1,designBLV_g1)
fitBLV_g1.2 <- contrasts.fit(fitBLV_g1,contrast.matrix)
fitBLV_g1.2 <- eBayes(fitBLV_g1.2)
tempOutput = topTable(fitBLV_g1.2, coef = 1, n=Inf)
nrDEGBLV_g1 = na.omit(tempOutput)
head(nrDEGBLV_g1)

fitBLG_g1 <- lmFit(log_blg_g1,designBLG_g1)
fitBLG_g1.2 <- contrasts.fit(fitBLG_g1,contrast.matrix)
fitBLG_g1.2 <- eBayes(fitBLG_g1.2)
tempOutput = topTable(fitBLG_g1.2, coef = 1, n=Inf)
nrDEGBLG_g1 = na.omit(tempOutput)
head(nrDEGBLG_g1)

fitBJIMT_g1 <- lmFit(log_bJIMT_g1,designBJIMT_g1)
fitBJIMT_g1.2 <- contrasts.fit(fitBJIMT_g1,contrast.matrix)
fitBJIMT_g1.2 <- eBayes(fitBJIMT_g1.2)
tempOutput = topTable(fitBJIMT_g1.2, coef = 1, n=Inf)
nrDEGBJIMT_g1 = na.omit(tempOutput)
head(nrDEGBJIMT_g1)

fitBHCC_g1 <- lmFit(log_bHCC_g1,designBHCC_g1)
fitBHCC_g1.2 <- contrasts.fit(fitBHCC_g1,contrast.matrix)
fitBHCC_g1.2 <- eBayes(fitBHCC_g1.2)
tempOutput = topTable(fitBHCC_g1.2, coef = 1, n=Inf)
nrDEGBHCC_g1 = na.omit(tempOutput)
head(nrDEGBHCC_g1)

fitBM_g1 <- lmFit(log_bm_g1,designBM_g1)
fitBM_g1.2 <- contrasts.fit(fitBM_g1,contrast.matrix)
fitBM_g1.2 <- eBayes(fitBM_g1.2)
tempOutput = topTable(fitBM_g1.2, coef = 1, n=Inf)
nrDEGBM_g1 = na.omit(tempOutput)
head(nrDEGBM_g1)

save(nrDEGBB_g1,nrDEGBK_g1,nrDEGBLV_g1,nrDEGBLG_g1, nrDEGBJIMT_g1, nrDEGBHCC_g1, file = "DEG_g1.Rdata")
save(exp_bb_g1,exp_bk_g1,exp_blg_g1,exp_blv_g1, exp_bJIMT_g1, exp_bHCC_g1,file = "exp_g1.Rdata")
save(log_bb_g1,log_bk_g1,log_blg_g1,log_blv_g1,log_bJIMT_g1, log_bHCC_g1,file = "log_g1.Rdata")

rm(list = ls())

load("DEG_g1.Rdata")
load("exp_g1.Rdata")
load("log_g1.Rdata")
load("group_g1.Rdata")

gl_bb_g1[gl_bb_g1 == "control"] <- "brain"
gl_bb_g1[gl_bb_g1 == "treat"] <- "bone"

gl_bk_g1[gl_bk_g1 == "control"] <- "brain"
gl_bk_g1[gl_bk_g1 == "treat"] <- "kidney"

gl_blg_g1[gl_blg_g1 == "control"] <- "brain"
gl_blg_g1[gl_blg_g1 == "treat"] <- "lung"

gl_blv_g1[gl_blv_g1 == "control"] <- "brain"
gl_blv_g1[gl_blv_g1 == "treat"] <- "liver"

gl_bJIMT_g1[gl_bJIMT_g1 == "control"] <- "brain"
gl_bJIMT_g1[gl_bJIMT_g1 == "treat"] <- "JIMT1"

gl_bHCC_g1[gl_bHCC_g1 == "control"] <- "brain"
gl_bHCC_g1[gl_bHCC_g1 == "treat"] <- "HCC1806"

gl_bm_g1[gl_bm_g1 == "control"] <- "brain"
gl_bm_g1[gl_bm_g1 == "treat"] <- "mammary_tumor"


new_colnames <- lapply(seq_along(exp_bb_g1), function(i) {
  paste0(gl_bb_g1[i], "_", colnames(exp_bb_g1)[i])
})
colnames(exp_bb_g1) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_bk_g1), function(i) {
  paste0(gl_bk_g1[i], "_", colnames(exp_bk_g1)[i])
})
colnames(exp_bk_g1) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_blv_g1), function(i) {
  paste0(gl_blv_g1[i], "_", colnames(exp_blv_g1)[i])
})
colnames(exp_blv_g1) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_blg_g1), function(i) {
  paste0(gl_blg_g1[i], "_", colnames(exp_blg_g1)[i])
})
colnames(exp_blg_g1) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_bJIMT_g1), function(i) {
  paste0(gl_bJIMT_g1[i], "_", colnames(exp_bJIMT_g1)[i])
})
colnames(exp_bJIMT_g1) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_bHCC_g1), function(i) {
  paste0(gl_bHCC_g1[i], "_", colnames(exp_bHCC_g1)[i])
})
colnames(exp_bHCC_g1) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_bm_g1), function(i) {
  paste0(gl_bm_g1[i], "_", colnames(exp_bm_g1)[i])
})
colnames(exp_bm_g1) <- unlist(new_colnames)

exp_bb_g1$logAve_brain <- rowMeans(log_bb_g1[,1:6], na.rm = T)
exp_bb_g1$logAve_bone <-rowMeans(log_bb_g1[,7:14], na.rm = T)
sorted_row <- order(rownames(exp_bb_g1))
exp_bb_g1 <- exp_bb_g1[sorted_row,]
sorted_row <- order(rownames(nrDEGBB_g1))
nrDEGBB_g1 <- nrDEGBB_g1[sorted_row,]
bb_g1 <- cbind(exp_bb_g1,nrDEGBB_g1)

exp_bk_g1$logAve_brain <- rowMeans(log_bk_g1[,1:6], na.rm = T)
exp_bk_g1$logAve_kidney <-rowMeans(log_bk_g1[,7:10], na.rm = T)
sorted_row <- order(rownames(exp_bk_g1))
exp_bk_g1 <- exp_bk_g1[sorted_row,]
sorted_row <- order(rownames(nrDEGBK_g1))
nrDEGBK_g1 <- nrDEGBK_g1[sorted_row,]
bk_g1 <- cbind(exp_bk_g1,nrDEGBK_g1)

exp_blv_g1$logAve_brain <- rowMeans(log_blv_g1[,1:6], na.rm = T)
exp_blv_g1$logAve_liver <-rowMeans(log_blv_g1[,7:13], na.rm = T)
sorted_row <- order(rownames(exp_blv_g1))
exp_blv_g1 <- exp_blv_g1[sorted_row,]
sorted_row <- order(rownames(nrDEGBLV_g1))
nrDEGBLV_g1 <- nrDEGBLV_g1[sorted_row,]
blv_g1 <- cbind(exp_blv_g1,nrDEGBLV_g1)

exp_blg_g1$logAve_brain <- rowMeans(log_blg_g1[,1:6], na.rm = T)
exp_blg_g1$logAve_lung <-rowMeans(log_blg_g1[,7:11], na.rm = T)
sorted_row <- order(rownames(exp_blg_g1))
exp_blg_g1 <- exp_blg_g1[sorted_row,]
sorted_row <- order(rownames(nrDEGBLG_g1))
nrDEGBLG_g1 <- nrDEGBLG_g1[sorted_row,]
blg_g1 <- cbind(exp_blg_g1,nrDEGBLG_g1)

exp_bJIMT_g1$logAve_brain <- rowMeans(log_bJIMT_g1[,1:4], na.rm = T)
exp_bJIMT_g1$logAve_JIMT1 <-rowMeans(log_bJIMT_g1[,5:6], na.rm = T)
sorted_row <- order(rownames(exp_bJIMT_g1))
exp_bJIMT_g1 <- exp_bJIMT_g1[sorted_row,]
sorted_row <- order(rownames(nrDEGBJIMT_g1))
nrDEGBJIMT_g1 <- nrDEGBJIMT_g1[sorted_row,]
bJIMT_g1 <- cbind(exp_bJIMT_g1,nrDEGBJIMT_g1)

exp_bHCC_g1$logAve_brain <- rowMeans(log_bHCC_g1[,1:2], na.rm = T)
exp_bHCC_g1$logAve_HCC1806 <-rowMeans(log_bHCC_g1[,3:4], na.rm = T)
sorted_row <- order(rownames(exp_bHCC_g1))
exp_bHCC_g1 <- exp_bHCC_g1[sorted_row,]
sorted_row <- order(rownames(nrDEGBHCC_g1))
nrDEGBHCC_g1 <- nrDEGBHCC_g1[sorted_row,]
bHCC_g1 <- cbind(exp_bHCC_g1,nrDEGBHCC_g1)

exp_bm_g1$logAve_brain <- rowMeans(log_bm_g1[,1:4], na.rm = T)
exp_bm_g1$logAve_mammary <-rowMeans(log_bm_g1[,5:8], na.rm = T)
sorted_row <- order(rownames(exp_bm_g1))
exp_bm_g1 <- exp_bm_g1[sorted_row,]
sorted_row <- order(rownames(nrDEGBM_g1))
nrDEGBM_g1 <- nrDEGBM_g1[sorted_row,]
bm_g1 <- cbind(exp_bm_g1,nrDEGBM_g1)

library(openxlsx)
sheets <- list('brain_bone' = bb_g1, "brain_kidney" = bk_g1, "brain_liver" = blv_g1,
               "brain_lung" = blg_g1, "brain_JIMT1" = bJIMT_g1, "brain_HCC1806" = bHCC_g1)


write.xlsx(sheets, "/Users/xujingyi/Documents/szbl/metastatic cell genomic/results/met map/g1.xlsx",rowNames = T)
write.xlsx(bm_g1, "/Users/xujingyi/Documents/szbl/metastatic cell genomic/results/met map/g1_mammary.xlsx",rowNames = T)







