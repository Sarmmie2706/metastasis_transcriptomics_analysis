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
                        save.filename = 'TCGA_SKCM.rda')
View(skcm_data)


# Count matrix
expression_count <- assay(skcm_data)
View(expression_count)


#Sample_info
sample_info <- colData(skcm_data) |>
  as.data.frame() |>
  apply(2, as.character) |>
  as.data.frame()
View(sample_info)

#Subsetting sample_info and removing rows with NA
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

#Save your datasets into scv files
write.csv(working_sample_info, "working_sample_info.csv")
write.csv(expression_count, "expression_count.csv")






