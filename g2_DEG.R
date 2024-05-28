setwd("/Users/xujingyi/Documents/szbl/metastatic cell genomic/code/met map")
rm(list = ls())

library(readxl)
library(limma)
library(edgeR)
data148283 <- read.csv("/Users/xujingyi/Documents/szbl/metastatic cell genomic/dataset/a metastasis map/GSE148283_all.count.csv",header=1)
rownames(data148283) <- data148283[,1]
data148283 <- data148283[,-1]


organ <- c("brain","bone","kidney","liver","lung")
information <- read_excel("/Users/xujingyi/Documents/szbl/metastatic cell genomic/dataset/a metastasis map/met map gene set.xlsx")
information <- information[information$Dataset_id == "GSE148283" & information$tissue_or_cell_type %in% organ,]
data148283 <- data148283[colnames(data148283) %in% information$Sample_name]
sorted_colnames <- sort(colnames(data148283))
data148283 <- data148283[,order(colnames(data148283))]
colnames(data148283) <- information$Sample_id
save(data148283,file = "data148283.Rdata")
save(information, file = "information.Rdata")

# Separate gene sets
rm(list = ls())
load("data148283.Rdata")
load("information.Rdata")
brain_g2 <- data148283[,information[information$Sample_characteristics == "group 2 (g2)" & information$tissue_or_cell_type == "brain",]$Sample_id]
bone_g2 <- data148283[,information[information$Sample_characteristics == "group 2 (g2)" & information$tissue_or_cell_type == "bone",]$Sample_id]
kidney_g2 <- data148283[,information[information$Sample_characteristics == "group 2 (g2)" & information$tissue_or_cell_type == "kidney",]$Sample_id]
liver_g2 <- data148283[,information[information$Sample_characteristics == "group 2 (g2)" & information$tissue_or_cell_type == "liver",]$Sample_id]
lung_g2 <- data148283[,information[information$Sample_characteristics == "group 2 (g2)" & information$tissue_or_cell_type == "lung",]$Sample_id]

# Get group information
gl_bb_g2 <- c(rep("control",ncol(brain_g2)),rep("treat",ncol(bone_g2)))
gl_bk_g2 <- c(rep("control",ncol(brain_g2)),rep("treat",ncol(kidney_g2)))
gl_blv_g2 <- c(rep("control",ncol(brain_g2)),rep("treat",ncol(liver_g2)))
gl_blg_g2 <- c(rep("control",ncol(brain_g2)),rep("treat",ncol(lung_g2)))

save(gl_bb_g2,gl_bk_g2,gl_blv_g2,gl_blg_g2, file = "group_g2.Rdata")

# Matrix preparation
options(stringsAsFactors = F)

exp_bb_g2 <- cbind(brain_g2,bone_g2)
log_bb_g2 <- log2(exp_bb_g2+1)
designBB_g2 <- model.matrix(~0+factor(gl_bb_g2))
colnames(designBB_g2) = levels(factor(gl_bb_g2))
rownames(designBB_g2) = colnames(log_bb_g2)
designBB_g2

exp_bk_g2 <- cbind(brain_g2,kidney_g2)
log_bk_g2 <- log2(exp_bk_g2+1)
designBK_g2 <- model.matrix(~0+factor(gl_bk_g2))
colnames(designBK_g2) = levels(factor(gl_bk_g2))
rownames(designBK_g2) = colnames(log_bk_g2)
designBK_g2

exp_blv_g2 <- cbind(brain_g2,liver_g2)
log_blv_g2 <- log2(exp_blv_g2+1)
designBLV_g2 <- model.matrix(~0+factor(gl_blv_g2))
colnames(designBLV_g2) = levels(factor(gl_blv_g2))
rownames(designBLV_g2) = colnames(log_blv_g2)
designBLV_g2

exp_blg_g2 <- cbind(brain_g2,lung_g2)
log_blg_g2 <- log2(exp_blg_g2+1)
designBLG_g2 <- model.matrix(~0+factor(gl_blg_g2))
colnames(designBLG_g2) = levels(factor(gl_blg_g2))
rownames(designBLG_g2) = colnames(log_blg_g2)
designBLG_g2

contrast.matrix <- makeContrasts(paste0(c("treat","control"),collapse = "-"),levels = designBB_g2)
contrast.matrix


# Differentially expressed gene analysis with limma
fitBB_g2 <- lmFit(log_bb_g2,designBB_g2)
fitBB_g2.2 <- contrasts.fit(fitBB_g2,contrast.matrix)
fitBB_g2.2 <- eBayes(fitBB_g2.2)
tempOutput = topTable(fitBB_g2.2, coef = 1, n=Inf)
nrDEGBB_g2 = na.omit(tempOutput)
head(nrDEGBB_g2)

fitBK_g2 <- lmFit(log_bk_g2,designBK_g2)
fitBK_g2.2 <- contrasts.fit(fitBK_g2,contrast.matrix)
fitBK_g2.2 <- eBayes(fitBK_g2.2)
tempOutput = topTable(fitBK_g2.2, coef = 1, n=Inf)
nrDEGBK_g2 = na.omit(tempOutput)
head(nrDEGBK_g2)

fitBLV_g2 <- lmFit(log_blv_g2,designBLV_g2)
fitBLV_g2.2 <- contrasts.fit(fitBLV_g2,contrast.matrix)
fitBLV_g2.2 <- eBayes(fitBLV_g2.2)
tempOutput = topTable(fitBLV_g2.2, coef = 1, n=Inf)
nrDEGBLV_g2 = na.omit(tempOutput)
head(nrDEGBLV_g2)

fitBLG_g2 <- lmFit(log_blg_g2,designBLG_g2)
fitBLG_g2.2 <- contrasts.fit(fitBLG_g2,contrast.matrix)
fitBLG_g2.2 <- eBayes(fitBLG_g2.2)
tempOutput = topTable(fitBLG_g2.2, coef = 1, n=Inf)
nrDEGBLG_g2 = na.omit(tempOutput)
head(nrDEGBLG_g2)


save(nrDEGBB_g2,nrDEGBK_g2,nrDEGBLV_g2,nrDEGBLG_g2, file = "DEG_g2.Rdata")
save(exp_bb_g2,exp_bk_g2,exp_blg_g2,exp_blv_g2, file = "exp_g2.Rdata")
save(log_bb_g2,log_bk_g2,log_blg_g2,log_blv_g2, file = "log_g2.Rdata")

rm(list = ls())

load("DEG_g2.Rdata")
load("exp_g2.Rdata")
load("log_g2.Rdata")
load("group_g2.Rdata")

gl_bb_g2[gl_bb_g2 == "control"] <- "brain"
gl_bb_g2[gl_bb_g2 == "treat"] <- "bone"

gl_bk_g2[gl_bk_g2 == "control"] <- "brain"
gl_bk_g2[gl_bk_g2 == "treat"] <- "kidney"

gl_blg_g2[gl_blg_g2 == "control"] <- "brain"
gl_blg_g2[gl_blg_g2 == "treat"] <- "lung"

gl_blv_g2[gl_blv_g2 == "control"] <- "brain"
gl_blv_g2[gl_blv_g2 == "treat"] <- "liver"


new_colnames <- lapply(seq_along(exp_bb_g2), function(i) {
  paste0(gl_bb_g2[i], "_", colnames(exp_bb_g2)[i])
})
colnames(exp_bb_g2) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_bk_g2), function(i) {
  paste0(gl_bk_g2[i], "_", colnames(exp_bk_g2)[i])
})
colnames(exp_bk_g2) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_blv_g2), function(i) {
  paste0(gl_blv_g2[i], "_", colnames(exp_blv_g2)[i])
})
colnames(exp_blv_g2) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_blg_g2), function(i) {
  paste0(gl_blg_g2[i], "_", colnames(exp_blg_g2)[i])
})
colnames(exp_blg_g2) <- unlist(new_colnames)

exp_bb_g2$logAve_brain <- rowMeans(log_bb_g2[,1:5], na.rm = T)
exp_bb_g2$logAve_bone <-rowMeans(log_bb_g2[,6:10], na.rm = T)
sorted_row <- order(rownames(exp_bb_g2))
exp_bb_g2 <- exp_bb_g2[sorted_row,]
sorted_row <- order(rownames(nrDEGBB_g2))
nrDEGBB_g2 <- nrDEGBB_g2[sorted_row,]
bb_g2 <- cbind(exp_bb_g2,nrDEGBB_g2)

exp_bk_g2$logAve_brain <- rowMeans(log_bk_g2[,1:5], na.rm = T)
exp_bk_g2$logAve_kidney <-rowMeans(log_bk_g2[,6:8], na.rm = T)
sorted_row <- order(rownames(exp_bk_g2))
exp_bk_g2 <- exp_bk_g2[sorted_row,]
sorted_row <- order(rownames(nrDEGBK_g2))
nrDEGBK_g2 <- nrDEGBK_g2[sorted_row,]
bk_g2 <- cbind(exp_bk_g2,nrDEGBK_g2)

exp_blv_g2$logAve_brain <- rowMeans(log_blv_g2[,1:5], na.rm = T)
exp_blv_g2$logAve_liver <-rowMeans(log_blv_g2[,6:12], na.rm = T)
sorted_row <- order(rownames(exp_blv_g2))
exp_blv_g2 <- exp_blv_g2[sorted_row,]
sorted_row <- order(rownames(nrDEGBLV_g2))
nrDEGBLV_g2 <- nrDEGBLV_g2[sorted_row,]
blv_g2 <- cbind(exp_blv_g2,nrDEGBLV_g2)

exp_blg_g2$logAve_brain <- rowMeans(log_blg_g2[,1:5], na.rm = T)
exp_blg_g2$logAve_lung <-rowMeans(log_blg_g2[,6:12], na.rm = T)
sorted_row <- order(rownames(exp_blg_g2))
exp_blg_g2 <- exp_blg_g2[sorted_row,]
sorted_row <- order(rownames(nrDEGBLG_g2))
nrDEGBLG_g2 <- nrDEGBLG_g2[sorted_row,]
blg_g2 <- cbind(exp_blg_g2,nrDEGBLG_g2)

library(openxlsx)
sheets <- list('brain_bone' = bb_g2, "brain_kidney" = bk_g2, "brain_liver" = blv_g2,
               "brain_lung" = blg_g2)


write.xlsx(sheets, "/Users/xujingyi/Documents/szbl/metastatic cell genomic/results/met map/g2.xlsx",rowNames = T)







 