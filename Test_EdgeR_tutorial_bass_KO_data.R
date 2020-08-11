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




