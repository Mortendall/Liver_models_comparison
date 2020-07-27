library(openxlsx)
library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(clipr)
library(limma)
library(tidyverse)
library(clusterProfiler)
library(biomaRt)
library(org.Mm.eg.db)
library(readxl)
#import data from 12W HFHS cohort

Western_diet_data <- read_excel("12W_western_silva_2018.xlsx")

Western_diet_data_sig <- Western_diet_data %>%
  filter(FDR < 0.05)
View(Western_diet_data_sig)

western_enztrez= bitr(Western_diet_data_sig$Gene_Symbol, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = "org.Mm.eg.db",
                    drop = T)

western_enztrez_bg= bitr(Western_diet_data$Gene_Symbol, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Mm.eg.db",
                      drop = T)

goResults_western <- enrichGO(gene = western_enztrez$ENTREZID,
                            universe = western_enztrez_bg$ENTREZID,
                            OrgDb = org.Mm.eg.db,
                            ont = "BP")
dotplot(goResults_western)
goResults_western <- setReadable(goResults_western, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
cnetplot(goResults_western)



