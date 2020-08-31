library(data.table)
library(magrittr)
library(edgeR)
library(stringr)
library(clusterProfiler)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(pheatmap)
library("RColorBrewer")
library("PoiClaClu")
library(org.Mm.eg.db)
library(tidyverse)
library(SummarizedExperiment)


setwd("C:/Users/tvb217/Documents/R/Comparison of liver models/Bass 2020")

redirectToFile <- function(y, file, fun, last = TRUE, ...)
{
  pdf(file, width = 4, height = 4) 
  fun(y, ...) 
  if(last) dev.off()
}

counts <- read.xlsx("Data_ZT8_Nampt_ko.xlsx")
group <- read.xlsx("Metadata_ZT8_Nampt_KO.xlsx", colNames = T)
group <- group %>%
  separate(Metadata, c("Sample_ID","ZT","Genotype"), sep = " ")
rownames(group) <- group$Sample_ID
group <- group %>%
  dplyr::select(Genotype)
group_factor <- factor(c(1,1,1,2,2,2,2))

countsMat <- counts[, -1] %>% as.matrix
rownames(countsMat) <- counts$Geneid
View(countsMat)
View(countsMat)
y = DGEList(counts = countsMat, group = group_factor)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = F]
y <- calcNormFactors(y)
design <- model.matrix(~group_factor)
y <- estimateDisp(y, design)%T>%
  redirectToFile(file.path("qcFigures", "BCV.pdf"), plotBCV)

efit <- glmQLFit(y, design) %T>%
  redirectToFile(file.path("qcFigures", "QLDispersion.pdf"), plotQLDisp)
colnames(design)[2] <- "Genotype"

dgeResults <- glmQLFTest(glmfit = efit, coef = 2)
topTags(dgeResults, n = Inf, p.value = 1) %>% 
  extract2("table") %>% 
  as.data.table(keep.rownames = TRUE)

topTags(dgeResults)
View(dgeResults)
dgeResults_export <- topTags(dgeResults, n = Inf, p.value = 1) %>% 
  extract2("table") %>% 
  as.data.table(keep.rownames = TRUE)

dir.create("edgeR_results", showWarnings = FALSE)
write.xlsx(dgeResults_export, file = "edgeR_results/Bass_WT_vs_KO.xlsx", asTable = TRUE)
View(dgeResults_export)
####GO analysis####


KO_control_list <- dgeResults_export
KO_control_list_significant <- KO_control_list %>%
  filter(FDR<0.05) 

KO_list <- KO_control_list_significant[,1]

KO_background <- KO_control_list[,1]

eg= bitr(KO_list$rn, 
         fromType = "SYMBOL", 
         toType = "ENTREZID", 
         OrgDb = "org.Mm.eg.db",
         drop = T)


bg = bitr(KO_background$rn, 
          fromType = "SYMBOL", 
          toType = "ENTREZID", 
          OrgDb = "org.Mm.eg.db",
          drop = T)
View(eg)
goResults <- enrichGO(gene = eg$ENTREZID,
                      universe = bg$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")
head(goResults@result)
dotplot(goResults)
goResults <- setReadable(goResults, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(goResults)



?redirectToFile
####Import of datafile####
edgeR_results <- read.xlsx("C:/Users/tvb217/Documents/R/Comparison of liver models/Bass 2020/edgeR_results/Bass_WT_vs_KO.xlsx", i)

####Test of NR effect####
setwd("C:/Users/tvb217/Documents/R/Comparison of liver models/Bass 2020")

redirectToFile <- function(y, file, fun, last = TRUE, ...)
{
  pdf(file, width = 4, height = 4) 
  fun(y, ...) 
  if(last) dev.off()
}

counts_NR <- read.xlsx("Data_ZT8.xlsx")
group_NR <- read.xlsx("Metadata_ZT8.xlsx", colNames = T)
group_NR <- group_NR %>%
  separate(Metadata, c("Sample_ID","ZT","Genotype", "Age", "Treatment"), sep = " ")
View(group_NR)
rownames(group_NR) <- group_NR$Sample_ID
group_NR <- group_NR %>%
  dplyr::select(Treatment)
group_factor <- factor(c(1,1,1,2,2,2))

countsMat_NR <- counts_NR[, -1] %>% as.matrix
rownames(countsMat_NR) <- counts_NR$Geneid
View(countsMat_NR)
y = DGEList(counts = countsMat_NR, group = group_factor)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = F]
y <- calcNormFactors(y)
design <- model.matrix(~group_factor)
y <- estimateDisp(y, design)%T>%
  redirectToFile(file.path("qcFigures", "BCV.pdf"), plotBCV)

efit <- glmQLFit(y, design) %T>%
  redirectToFile(file.path("qcFigures", "QLDispersion.pdf"), plotQLDisp)
colnames(design)[2] <- "Genotype"

dgeResults <- glmQLFTest(glmfit = efit, coef = 2)
topTags(dgeResults, n = Inf, p.value = 1) %>% 
  extract2("table") %>% 
  as.data.table(keep.rownames = TRUE)

topTags(dgeResults)

dgeResults_export <- topTags(dgeResults, n = Inf, p.value = 1) %>% 
  extract2("table") %>% 
  as.data.table(keep.rownames = TRUE)

dir.create("edgeR_results_NR", showWarnings = FALSE)
write.xlsx(dgeResults_export, file = "edgeR_results_NR/Bass_NR_Effect.xlsx", asTable = TRUE)
View(dgeResults_export)

#No terms are significant following FDR correction

####Testing the code for another setup at ZT=4####
counts <- read.xlsx("Data_ZT4_Nampt_ko.xlsx")
group <- read.xlsx("Metadata_ZT4_Nampt_KO.xlsx", colNames = T)
group <- group %>%
  separate(Metadata, c("Sample_ID","ZT","Genotype"), sep = " ")
View(group)
rownames(group) <- group$Sample_ID
group <- group %>%
  dplyr::select(Genotype)
group_factor <- factor(c(1,1,1,2,2,2))

countsMat <- counts[, -1] %>% as.matrix
rownames(countsMat) <- counts$Geneid
View(countsMat)
View(countsMat)
y = DGEList(counts = countsMat, group = group_factor)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = F]
y <- calcNormFactors(y)
design <- model.matrix(~group_factor)
y <- estimateDisp(y, design)%T>%
  redirectToFile(file.path("qcFigures", "BCV.pdf"), plotBCV)

efit <- glmQLFit(y, design) %T>%
  redirectToFile(file.path("qcFigures", "QLDispersion.pdf"), plotQLDisp)
colnames(design)[2] <- "Genotype"

dgeResults <- glmQLFTest(glmfit = efit, coef = 2)
topTags(dgeResults, n = Inf, p.value = 1) %>% 
  extract2("table") %>% 
  as.data.table(keep.rownames = TRUE)

topTags(dgeResults)
View(dgeResults)
dgeResults_export <- topTags(dgeResults, n = Inf, p.value = 1) %>% 
  extract2("table") %>% 
  as.data.table(keep.rownames = TRUE)


write.xlsx(dgeResults_export, file = "edgeR_results/Bass_WT_vs_KO_ZT4.xlsx", asTable = TRUE)
View(dgeResults_export)

####analyzing data for main effect of NR treatment####
setwd("C:/Users/tvb217/Documents/R/Comparison of liver models/Bass 2020")
counts <- read.xlsx("Data_NR_main.xlsx")
group <- read.xlsx("metadata_main_NR.xlsx", colNames = T)
group <- group %>%
  separate(Metadata, c("Sample_ID","ZT","Genotype", "age", "months", "Treatment"), sep = " ")

View(group)
rownames(group) <- group$Sample_ID
group <- group %>%
  dplyr::select(Treatment)
group_factor <- if_else(group$Treatment == "\"H2O\"", "1", "2")


countsMat <- counts[, -1] %>% as.matrix
rownames(countsMat) <- counts$Geneid

y = DGEList(counts = countsMat, group = group_factor)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = F]
y <- calcNormFactors(y)
design <- model.matrix(~group_factor)
y <- estimateDisp(y, design)%T>%
  redirectToFile(file.path("qcFigures", "BCV.pdf"), plotBCV)

efit <- glmQLFit(y, design) %T>%
  redirectToFile(file.path("qcFigures", "QLDispersion.pdf"), plotQLDisp)
colnames(design)[2] <- "Genotype"

dgeResults <- glmQLFTest(glmfit = efit, coef = 2)
topTags(dgeResults, n = Inf, p.value = 1) %>% 
  extract2("table") %>% 
  as.data.table(keep.rownames = TRUE)

topTags(dgeResults)
View(dgeResults)
dgeResults_export <- topTags(dgeResults, n = Inf, p.value = 1) %>% 
  extract2("table") %>% 
  as.data.table(keep.rownames = TRUE)


write.xlsx(dgeResults_export, file = "edgeR_results/Bass_NR_main.xlsx", asTable = TRUE)
View(dgeResults_export)

NR_list_data <- dgeResults_export
NR_list_significant <- NR_list %>%
  filter(FDR<0.05) 

NR_list <- NR_list_significant[,1]

NR_background <- NR_list_data[,1]

eg= bitr(NR_list$rn, 
         fromType = "SYMBOL", 
         toType = "ENTREZID", 
         OrgDb = "org.Mm.eg.db",
         drop = T)


bg = bitr(NR_background$rn, 
          fromType = "SYMBOL", 
          toType = "ENTREZID", 
          OrgDb = "org.Mm.eg.db",
          drop = T)
View(bg)
goResults <- enrichGO(gene = eg$ENTREZID,
                      universe = bg$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")
head(goResults@result)
dotplot(goResults)
goResults <- setReadable(goResults, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(goResults)
View(goResults@result)
