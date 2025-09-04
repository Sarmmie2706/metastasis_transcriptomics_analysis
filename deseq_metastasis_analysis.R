# Load necessary libraries
library(tidyverse)
library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)

# Building a query to retrieve primary and metastatic data
query <- GDCquery(
  project = 'TCGA-SKCM', # Skin cutaneous melanoma
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  access = 'open',
  sample.type = c('Primary Tumor', 'Metastatic') 
)

# View sample metadata before downloading (Optional)
output_query <- getResults(query)
View(output_query)

# download data
GDCdownload(query, files.per.chunk = 165)

# prepare data
## run to save prepared data as TCGA_BRCA.rda
skcm_data <- GDCprepare(query, 
                        summarizedExperiment = TRUE, #Default
                        save=TRUE, 
                        save.filename = 'TCGA_SKCM2.rda')
View(skcm_data)


# Retrieving Count matrix
expression_count <- assay(skcm_data)
View(expression_count)


# Retrieving Sample info
sample_info <- colData(skcm_data) |>
  as.data.frame() |>
  apply(2, as.character) |>
  as.data.frame()
View(sample_info)

#Subsetting sample_info
sample_info <- sample_info %>% 
  select(barcode,
         patient,
         definition,
         sample_type,
         tumor_descriptor,
         ajcc_pathologic_stage,
         vital_status,
         gender,
         race,
         ethnicity,
         age_at_index,
         days_to_death) %>% 
  mutate( #Change pathological stage to NA if tissue type is metastatic
    ajcc_pathologic_stage = ifelse(
      tumor_descriptor == "Metastatic",
      "NA",
      ajcc_pathologic_stage
    )
)

#Checking number of tumor and normal samples
table(sample_info$tumor_descriptor)

#Making the samples (columns) of expression data align with samples (rows) of sample_info
sample_info <- sample_info %>%
  rownames_to_column(var = "number_id") %>% #Converts numbered rownames to column
  mutate(barcode_copy = barcode) %>%       # Create a new column to keep barcodes
  column_to_rownames(var = "barcode") %>%  # Move original 'barcode' to rownames
  select(barcode_copy, everything(), -number_id) #Remove the newly formed numbered column and moves the barcode_copy to the start

# Selecting a third of both types of samples to reduce computational strain, you can skip all 
# these steps and save continue working with the sample_info variable

# Selecting and subsetting the samples
primary_tumor <- sample_info[sample_info$tumor_descriptor == "Primary", ]
metastatic_tumor <- sample_info[sample_info$tumor_descriptor == "Metastatic", ]

# Sampling and picking a third of each
primary_tumor <- primary_tumor[sample(nrow(primary_tumor), floor(nrow(primary_tumor) / 3)), ]
metastatic_tumor <- metastatic_tumor[sample(nrow(metastatic_tumor), floor(nrow(metastatic_tumor) / 3)), ]

# Joining both dataframes
working_sample_info <- rbind(primary_tumor, metastatic_tumor)

#Select the desired samples based on the selected samples for working_sample_info
expression_count <- expression_count[, rownames(working_sample_info)]
dim(expression_count)

#Checking if matching was done correctly and in the right order
table(rownames(working_sample_info) %in% colnames(expression_count))
table(rownames(working_sample_info) == colnames(expression_count))

#Save your datasets into rds files. This preserves structure and object types when loading back into R
saveRDS(working_sample_info, "working_sample_info.rds")
saveRDS(expression_count, "expression_count.rds")

#Read back your data into R
working_sample_info <- readRDS("working_sample_info.rds")
expression_count <- readRDS("expression_count.rds")

#Checking if matching was done correctly and in the right order
table(rownames(working_sample_info) %in% colnames(expression_count))
table(rownames(working_sample_info) == colnames(expression_count))

#Filtering out lowly expressed genes
expression_count <- expression_count[rowMeans(expression_count) >= 10, ]
dim(expression_count)

####################################################################################
######################## DIFFERENTIAL EXPRESSION ANALYSIS ##########################
####################################################################################

#Creating a DESeq2 Object
dds <- DESeqDataSetFromMatrix(
  countData = expression_count,
  colData = working_sample_info,
  design = ~ tumor_descriptor
)

#Setting reference levels
dds$tumor_descriptor <- relevel(dds$tumor_descriptor, ref = "Primary")

#Collapse technical replicates if you have them before running DESeq
#Running DESeq proper
dds <- DESeq(dds)

#Stabilizing your dds result. This will be used for PCA and Clustering later on
dds_norm <- varianceStabilizingTransformation(dds)

#Save and view result statistics
res <- results(dds)
res
summary(res)

saveRDS(res, file = "res.rds")

#Setting FDR threshold to 0.01
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

resultsNames(dds) #To check the different checked conditions

#In cases of more than two compared conditions and you want to check between two specific conditions
#e.g.: Primary, Metastatic, Recurrent
#Essentially saying log2foldchange = Metastatic - Recurrent
#results(dds, contrast = c("tumor_descriptor", "Metastatic", "Recurrent"))

#Filtering the most significantly differentiated genes
res <- res[res$padj < 0.01 & abs(res$log2FoldChange) > 2, ]

res_up <- res[order(res$log2FoldChange, decreasing = T), ]
top_upregulated_genes <- head(res_up, 20) %>% as.data.frame()
top_upregulated_genes <- rownames(top_upregulated_genes)

res_down <- res[order(res$log2FoldChange), ]
top_downregulated_genes <- head(res_down, 20) %>% as.data.frame()
top_downregulated_genes <- rownames(top_downregulated_genes)

#Saving your genes to csv files
write.csv(top_upregulated_genes, file = "top_upregulates_genelist.csv")
write.csv(top_downregulated_genes, file = "top_downregulates_genelist.csv")

####################################################################################
######################### FUNCTIONAL ENRICHMENT ANALYSIS ###########################
####################################################################################
library(clusterProfiler)
library(org.Hs.eg.db)  

# Extract gene IDs
genes <- rownames(res)

# Remove the Ensembl version numbers after the dot as bitr can't recognise them
genes <- gsub("\\..*", "", genes)

# Convert Ensembl IDs to Entrez IDs as enrichKEGG only accepts Entrez gene IDs
gene_conversion <- bitr(genes, 
                        fromType = "ENSEMBL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db
                        )
View(gene_conversion)

# Create gene list for enrichment
entrez_gene_list <- gene_conversion$ENTREZID

enrich_go <- enrichGO(gene          = entrez_gene_list,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",         # Biological Process, could also be Cellular Component (CC), or Molecular Function (MF)
                pAdjustMethod = "BH",         
                pvalueCutoff  = 0.05,         
                qvalueCutoff  = 0.05,         
                readable      = TRUE)  # This converts Entrez IDs back to ENSEMBL gene symbols

View(enrich_go)

saveRDS(enrich_go, "enrich_go_results.rds")




























