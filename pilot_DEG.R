setwd("/Users/xujingyi/Documents/szbl/metastatic cell genomic/code/met map")
rm(list = ls())

library(readxl)
library(limma)
library(edgeR)
data148283 <- read.csv("/Users/xujingyi/Documents/szbl/metastatic cell genomic/dataset/a metastasis map/GSE148283_all.count.csv",header=1)
rownames(data148283) <- data148283[,1]
data148283 <- data148283[,-1]
data148283[,"pilot.S57"] <- data148283[,"pilot.S5758"]
colnames(data148283)[colnames(data148283) == "pilot.S5758"] <- "pilot.S58"



organ <- c("brain","bone","kidney","liver","lung","HCC1954_BREAST","JIMT1_BREAST","HCC1806_BREAST","mammary_tumor")
information <- read_excel("/Users/xujingyi/Documents/szbl/metastatic cell genomic/dataset/a metastasis map/met map gene set.xlsx")
information <- information[information$Dataset_id == "GSE148283" & information$tissue_or_cell_type %in% organ,]
data148283 <- data148283[colnames(data148283) %in% information$Sample_name]
sorted_rownames <- sort(colnames(data148283))
data148283 <- data148283[,order(colnames(data148283))]
colnames(data148283) <- information$Sample_id
save(data148283,file = "data148283.Rdata")
save(information, file = "information.Rdata")

# Separate gene sets
rm(list = ls())
load("data148283.Rdata")
load("information.Rdata")
brain_pilot <- data148283[,information[information$Sample_characteristics == "pilot" & information$tissue_or_cell_type == "brain",]$Sample_id]
bone_pilot <- data148283[,information[information$Sample_characteristics == "pilot" & information$tissue_or_cell_type == "bone",]$Sample_id]
kidney_pilot <- data148283[,information[information$Sample_characteristics == "pilot" & information$tissue_or_cell_type == "kidney",]$Sample_id]
liver_pilot <- data148283[,information[information$Sample_characteristics == "pilot" & information$tissue_or_cell_type == "liver",]$Sample_id]
lung_pilot <- data148283[,information[information$Sample_characteristics == "pilot" & information$tissue_or_cell_type == "lung",]$Sample_id]
HCC_pilot <- data148283[,information[information$Sample_characteristics == "pilot" & information$tissue_or_cell_type == "HCC1954_BREAST",]$Sample_id]

# Get group information
gl_bb_pilot <- c(rep("control",ncol(brain_pilot)),rep("treat",ncol(bone_pilot)))
gl_bk_pilot <- c(rep("control",ncol(brain_pilot)),rep("treat",ncol(kidney_pilot)))
gl_blv_pilot <- c(rep("control",ncol(brain_pilot)),rep("treat",ncol(liver_pilot)))
gl_blg_pilot <- c(rep("control",ncol(brain_pilot)),rep("treat",ncol(lung_pilot)))
gl_b1954_pilot <- c(rep("control",ncol(brain_pilot)),rep("treat",ncol(HCC_pilot)))

save(gl_bb_pilot,gl_bk_pilot,gl_blv_pilot,gl_blg_pilot,gl_b1954_pilot, file = "group_pilot.Rdata")

# Matrix preparation
options(stringsAsFactors = F)

exp_bb_pilot <- cbind(brain_pilot,bone_pilot)
log_bb_pilot <- log2(exp_bb_pilot + 1)
designBB_pilot <- model.matrix(~0+factor(gl_bb_pilot))
colnames(designBB_pilot) = levels(factor(gl_bb_pilot))
rownames(designBB_pilot) = colnames(log_bb_pilot)
designBB_pilot

exp_bk_pilot <- cbind(brain_pilot,kidney_pilot)
log_bk_pilot <- log2(exp_bk_pilot + 1)
designBK_pilot <- model.matrix(~0+factor(gl_bk_pilot))
colnames(designBK_pilot) = levels(factor(gl_bk_pilot))
rownames(designBK_pilot) = colnames(log_bk_pilot)
designBK_pilot

exp_blv_pilot <- cbind(brain_pilot,liver_pilot)
log_blv_pilot <- log2(exp_blv_pilot + 1)
designBLV_pilot <- model.matrix(~0+factor(gl_blv_pilot))
colnames(designBLV_pilot) = levels(factor(gl_blv_pilot))
rownames(designBLV_pilot) = colnames(log_blv_pilot)
designBLV_pilot

exp_blg_pilot <- cbind(brain_pilot,lung_pilot)
log_blg_pilot <- log2(exp_blg_pilot + 1)
designBLG_pilot <- model.matrix(~0+factor(gl_blg_pilot))
colnames(designBLG_pilot) = levels(factor(gl_blg_pilot))
rownames(designBLG_pilot) = colnames(log_blg_pilot)
designBLG_pilot

exp_b1954_pilot <- cbind(brain_pilot,HCC_pilot)
log_b1954_pilot <- log2(exp_b1954_pilot + 1)
designB1954_pilot <- model.matrix(~0+factor(gl_b1954_pilot))
colnames(designB1954_pilot) = levels(factor(gl_b1954_pilot))
rownames(designB1954_pilot) = colnames(log_b1954_pilot)
designB1954_pilot

contrast.matrix <- makeContrasts(paste0(c("treat","control"),collapse = "-"),levels = designBB_pilot)
contrast.matrix


# Differentially expressed gene analysis with limma
fitBB_pilot <- lmFit(log_bb_pilot,designBB_pilot)
fitBB_pilot.2 <- contrasts.fit(fitBB_pilot,contrast.matrix)
fitBB_pilot.2 <- eBayes(fitBB_pilot.2)
tempOutput = topTable(fitBB_pilot.2, coef = 1, n=Inf)
nrDEGBB_pilot = na.omit(tempOutput)
head(nrDEGBB_pilot)

fitBK_pilot <- lmFit(log_bk_pilot,designBK_pilot)
fitBK_pilot.2 <- contrasts.fit(fitBK_pilot,contrast.matrix)
fitBK_pilot.2 <- eBayes(fitBK_pilot.2)
tempOutput = topTable(fitBK_pilot.2, coef = 1, n=Inf)
nrDEGBK_pilot = na.omit(tempOutput)
head(nrDEGBK_pilot)

fitBLV_pilot <- lmFit(log_blv_pilot,designBLV_pilot)
fitBLV_pilot.2 <- contrasts.fit(fitBLV_pilot,contrast.matrix)
fitBLV_pilot.2 <- eBayes(fitBLV_pilot.2)
tempOutput = topTable(fitBLV_pilot.2, coef = 1, n=Inf)
nrDEGBLV_pilot = na.omit(tempOutput)
head(nrDEGBLV_pilot)

fitBLG_pilot <- lmFit(log_blg_pilot,designBLG_pilot)
fitBLG_pilot.2 <- contrasts.fit(fitBLG_pilot,contrast.matrix)
fitBLG_pilot.2 <- eBayes(fitBLG_pilot.2)
tempOutput = topTable(fitBLG_pilot.2, coef = 1, n=Inf)
nrDEGBLG_pilot = na.omit(tempOutput)
head(nrDEGBLG_pilot)

fitB1954_pilot <- lmFit(log_b1954_pilot,designB1954_pilot)
fitB1954_pilot.2 <- contrasts.fit(fitB1954_pilot,contrast.matrix)
fitB1954_pilot.2 <- eBayes(fitB1954_pilot.2)
tempOutput = topTable(fitB1954_pilot.2, coef = 1, n=Inf)
nrDEGB1954_pilot = na.omit(tempOutput)
head(nrDEGB1954_pilot)

save(nrDEGBB_pilot,nrDEGBK_pilot,nrDEGBLV_pilot,nrDEGBLG_pilot, nrDEGB1954_pilot, file = "DEG_pilot.Rdata")
save(exp_bb_pilot,exp_bk_pilot,exp_blg_pilot,exp_blv_pilot, exp_b1954_pilot, file = "exp_pilot.Rdata")
save(log_bb_pilot,log_bk_pilot,log_blg_pilot,log_blv_pilot,log_b1954_pilot, file = "log_pilot.Rdata")

rm(list = ls())

load("DEG_pilot.Rdata")
load("exp_pilot.Rdata")
load("log_pilot.Rdata")
load("group_pilot.Rdata")

gl_bb_pilot[gl_bb_pilot == "control"] <- "brain"
gl_bb_pilot[gl_bb_pilot == "treat"] <- "bone"

gl_bk_pilot[gl_bk_pilot == "control"] <- "brain"
gl_bk_pilot[gl_bk_pilot == "treat"] <- "kidney"

gl_blg_pilot[gl_blg_pilot == "control"] <- "brain"
gl_blg_pilot[gl_blg_pilot == "treat"] <- "lung"

gl_blv_pilot[gl_blv_pilot == "control"] <- "brain"
gl_blv_pilot[gl_blv_pilot == "treat"] <- "liver"

gl_b1954_pilot[gl_b1954_pilot == "control"] <- "brain"
gl_b1954_pilot[gl_b1954_pilot == "treat"] <- "HCC1954"

new_colnames <- lapply(seq_along(exp_bb_pilot), function(i) {
  paste0(gl_bb_pilot[i], "_", colnames(exp_bb_pilot)[i])
})
colnames(exp_bb_pilot) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_bk_pilot), function(i) {
  paste0(gl_bk_pilot[i], "_", colnames(exp_bk_pilot)[i])
})
colnames(exp_bk_pilot) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_blv_pilot), function(i) {
  paste0(gl_blv_pilot[i], "_", colnames(exp_blv_pilot)[i])
})
colnames(exp_blv_pilot) <- unlist(new_colnames)

new_colnames <- lapply(seq_along(exp_blg_pilot), function(i) {
  paste0(gl_blg_pilot[i], "_", colnames(exp_blg_pilot)[i])
})
colnames(exp_blg_pilot) <- unlist(new_colnames)


new_colnames <- lapply(seq_along(exp_b1954_pilot), function(i) {
  paste0(gl_b1954_pilot[i], "_", colnames(exp_b1954_pilot)[i])
})
colnames(exp_b1954_pilot) <- unlist(new_colnames)

exp_bb_pilot$logAve_brain <- rowMeans(log_bb_pilot[,1:8], na.rm = T)
exp_bb_pilot$logAve_bone <-rowMeans(log_bb_pilot[,9:13], na.rm = T)
sorted_row <- order(rownames(exp_bb_pilot))
exp_bb_pilot <- exp_bb_pilot[sorted_row,]
sorted_row <- order(rownames(nrDEGBB_pilot))
nrDEGBB_pilot <- nrDEGBB_pilot[sorted_row,]
bb_pilot <- cbind(exp_bb_pilot,nrDEGBB_pilot)

exp_bk_pilot$logAve_brain <- rowMeans(log_bk_pilot[,1:8], na.rm = T)
exp_bk_pilot$logAve_kidney <-rowMeans(log_bk_pilot[,9:12], na.rm = T)
sorted_row <- order(rownames(exp_bk_pilot))
exp_bk_pilot <- exp_bk_pilot[sorted_row,]
sorted_row <- order(rownames(nrDEGBK_pilot))
nrDEGBK_pilot <- nrDEGBK_pilot[sorted_row,]
bk_pilot <- cbind(exp_bk_pilot,nrDEGBK_pilot)

exp_blv_pilot$logAve_brain <- rowMeans(log_blv_pilot[,1:8], na.rm = T)
exp_blv_pilot$logAve_liver <-rowMeans(log_blv_pilot[,9:14], na.rm = T)
sorted_row <- order(rownames(exp_blv_pilot))
exp_blv_pilot <- exp_blv_pilot[sorted_row,]
sorted_row <- order(rownames(nrDEGBLV_pilot))
nrDEGBLV_pilot <- nrDEGBLV_pilot[sorted_row,]
blv_pilot <- cbind(exp_blv_pilot,nrDEGBLV_pilot)

exp_blg_pilot$logAve_brain <- rowMeans(log_blg_pilot[,1:8], na.rm = T)
exp_blg_pilot$logAve_lung <-rowMeans(log_blg_pilot[,9:14], na.rm = T)
sorted_row <- order(rownames(exp_blg_pilot))
exp_blg_pilot <- exp_blg_pilot[sorted_row,]
sorted_row <- order(rownames(nrDEGBLG_pilot))
nrDEGBLG_pilot <- nrDEGBLG_pilot[sorted_row,]
blg_pilot <- cbind(exp_blg_pilot,nrDEGBLG_pilot)

exp_b1954_pilot$logAve_brain <- rowMeans(log_b1954_pilot[,1:8], na.rm = T)
exp_b1954_pilot$logAve_HCC1954 <-rowMeans(log_b1954_pilot[,9:11], na.rm = T)
sorted_row <- order(rownames(exp_b1954_pilot))
exp_b1954_pilot <- exp_b1954_pilot[sorted_row,]
sorted_row <- order(rownames(nrDEGB1954_pilot))
nrDEGB1954_pilot <- nrDEGB1954_pilot[sorted_row,]
b1954_pilot <- cbind(exp_b1954_pilot,nrDEGB1954_pilot)

library(openxlsx)
sheets <- list('brain_bone' = bb_pilot, "brain_kidney" = bk_pilot, "brain_liver" = blv_pilot,
               "brain_lung" = blg_pilot, "brain_HCC1954" = b1954_pilot)


write.xlsx(sheets, "/Users/xujingyi/Documents/szbl/metastatic cell genomic/results/met map/pilot.xlsx",rowNames = T)







