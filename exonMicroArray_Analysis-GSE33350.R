# The following code is for Exon microarray analysis (normalisation, QC, differential gene expression and gene ontology)...
# ... for the dataset GSE33350: https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/query/acc.cgi?acc=GSE33350

# Installs/library calls ----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("oligo")

library(oligo) # package designed specifically for exon microarray analysis
library(GEOquery)
library(affycoretools)
library(dplyr)
library(tibble)
library(preprocessCore)
library(affy)
library(vsn)

# Background subtraction and normalisation to probeset level ----
getGEOSuppFiles("GSE33350") # downloads .cel files for downstream analysis
untar("M:/GSE33350/GSE33350/GSE33350_RAW.tar", list = F)
celFiles <- oligoClasses::list.celfiles('M:/GSE33350/', listGzipped = T) # enter in full file directory here
raw.data.oligo <- read.celfiles(celFiles)
gse <- getGEO("GSE33350") # we'll mine this for phenotypic data later on.

oligo::boxplot(raw.data.oligo, target = "probeset", transfo=log2) #pre-processed Log2 transformed probset data

sampleNames <- sub(".CEL.gz$", "", celFiles)
sampleNames(raw.data.oligo) <- sampleNames

normDataprobeset <- oligo::rma(raw.data.oligo, background = T, target = "probeset", normalize = T) # Uses quantile normalisation
normDataprobeset <- affycoretools::annotateEset(normDataprobeset, pd.hugene.1.0.st.v1, type = "probeset")

gene_data_frame = fData(normDataprobeset)
gene_data_frame <- gene_data_frame[!is.na(gene_data_frame$SYMBOL),]

normDataAnno <- as.data.frame(normDataprobeset@assayData$exprs) %>%
  rownames_to_column(var = "PROBEID")
normDataAnno <- dplyr::inner_join(normDataAnno, gene_data_frame, by = "PROBEID")

# Subset top 50% of probesets with highest interquartile range across all samples, as per the methodology of the original paper 
pb = txtProgressBar(min = 1, max = length(unique(normDataAnno$SYMBOL)), initial = 0, char = "=") #creating a progress bar
normDataGeneExprs <- data.frame(matrix(nrow= length(unique(normDataAnno$SYMBOL)),
                                       ncol = length(c(sampleNames, "SYMBOL"))))
runTime = 0

# here we create a blank dataframe to make the for loop faster
colnames(normDataGeneExprs) <- c(sampleNames, "SYMBOL")
normDataGeneExprs$SYMBOL <- BiocGenerics::unique(normDataAnno$SYMBOL)
normDataGeneExprs <- normDataGeneExprs %>%
  column_to_rownames("SYMBOL")

# For loop takes the upper interquartile range of probesets and performed median polish summarisation
for(i in unique(normDataAnno$SYMBOL)){
  
  tmp <- normDataAnno[normDataAnno$SYMBOL %in% i,]
  tmp$aveExprs <- apply(as.matrix(tmp[ , sampleNames]), 1, mean)
  
  if(as.numeric(nrow(tmp)) >= 3){
    tmp <- tmp %>%
      filter(aveExprs >= quantile(aveExprs, .5) & aveExprs <= quantile(aveExprs, .75))
  } else {
    tmp <- tmp %>%
    filter(aveExprs >= quantile(aveExprs, .5))
  }
  
  tmp <- as.data.frame(tmp)
  
  normDataGeneExprs[i, sampleNames] <- colSummarizeMedianpolish(as.matrix(tmp[, sampleNames]))$Estimates
  
  runTime = runTime + 1
  setTxtProgressBar(pb, runTime)
}  # This loop takes a long time. Go get a cuppa.
close(pb)
rm(i, runTime, tmp, pb, runTime)

normDataGeneExprs["ASS1" ,] # shows data for ASS1

# QC QA plots to check normalisation ----
# meanSD plot should show horizontal line since we assume that variance among expression of all genes...
# is equal since most are not differentially expressed
# see for explaination -> https://www.bioconductor.org/packages/devel/bioc/vignettes/vsn/inst/doc/A-vsn.html

fitQC_raw <- vsnMatrix(as.matrix(raw.data.oligo@assayData[["exprs"]]))
meanSdPlot(fitQC_raw) # meanSDplot of raw data
ggsave("meanSD_raw.png")

fitQC_norm <- vsnMatrix(as.matrix(normDataGeneExprs))
meanSdPlot(fitQC_norm)
ggsave("meanSD_normalisedGeneExprs.png")

boxplot(raw.data.oligo, target = "core")
ggsave("boxplot_raw.png")
boxplot(normDataGeneExprs, target = "core")
ggsave("boxplot_norm.png")

hist(raw.data.oligo, target = "core")
ggsave("hist_raw.png")
hist(normDataprobeset, target = "core")
ggsave("hist_norm.png")

# LIMMA for differential gene expression ----
library(pheatmap)
library(limma)

eset4limma <- normDataGeneExprs
groups <- pData(gse[[1]])
groups <- data.frame(sampleNames, groups[ "source_name_ch1"])
groups$source_name_ch1 <- gsub(" ", "", groups$source_name_ch1)# creating group name vector. Check sample IDs vs groups.
groupsF <- factor(groups[ , "source_name_ch1"])


mm <- model.matrix(~0 + groupsF)
pheatmap(mm, cluster_rows = F, cluster_cols = F) # double checks which group each sample belongs to.

fit <- lmFit(eset4limma, mm)
coef.fit <- fit$coefficients # fits to linear model to generate a summary average expression value for the biological replicates
colnames(coef.fit) <- gsub("groupsF", "", colnames(coef.fit))

contr <- makeContrasts(Lungmetastasisderivative - Parental, levels = groupsF)
colnames(mm) <- colnames(coef.fit) # tries to match column names to avoid warning message later on but is unsuccessful. To fix later.
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

x = 1
for(i in top.table$logFC){
  
  top.table[x, "FC"] <- 2^i
  x <- x+1
} # unlogs fold change values
rm(x, i)

head(top.table)

# GSEA analysis ----
library(pathfindR) # alternate packages include fgsea
top.table <- top.table %>% 
  rownames_to_column("SYMBOL")

gseaTopTablefilter <- top.table[, c("SYMBOL", "logFC", "adj.P.Val")]
gseaTopTablefilter <- gseaTopTablefilter[1:1000, ] # Picks up to 1000 differentially expressed genes.

# The run pathfindR function takes a long time.
output_KEGG <- run_pathfindR(gseaTopTablefilter, p_val_threshold = 0.05, gene_sets = "KEGG")
output_Reactome <- run_pathfindR(gseaTopTablefilter, p_val_threshold = 0.05, gene_sets = "Reactome")
output_GO_BioCarta <- run_pathfindR(gseaTopTablefilter, p_val_threshold = 0.05, gene_sets = "BioCarta")
output_GO_All <- run_pathfindR(gseaTopTablefilter, p_val_threshold = 0.05, gene_sets = "GO-All")
output_Biological_Process <- run_pathfindR(gseaTopTablefilter, p_val_threshold = 0.05, gene_sets = "GO-BP")
output_Molecular_Function <- run_pathfindR(gseaTopTablefilter, p_val_threshold = 0.05, gene_sets = "GO-MF")

enrichment_chart(output_KEGG[output_KEGG$ID %in%
                                    c("hsa00910", "hsa00320", "hsa01040", "hsa00330", "hsa04927", "hsa04072", "hsa04010"), ], top_terms = Inf)

library(openxlsx)
df_names <- list("DEG" = top.table, "output_KEGG" = output_KEGG, "output_Reactome" = output_Reactome, "output_GO_BioCarta" = output_GO_BioCarta , 
                 "output_GO_All" = output_GO_All, "output_GO_Biological_Process" = output_Biological_Process,
                 "output_Molecular_Function" = output_Molecular_Function)
openxlsx::write.xlsx(df_names, file = "Parental-Lung_GSEA_Analysis.xlsx") # writes the pathfindR outputs to an excel sheet with multiple tabs


