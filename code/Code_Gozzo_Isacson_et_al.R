
# This is the R code used for Gozzo & Isacson et al. 2026
# Here is included omics and experimental analysis contained in the manuscript. 
# If you have any questions or concerns, please reach out to me at
# alessandro.gozzo@liu.se or here on GitHub
#
#Across the code Piwi/Aub/GFP refer respectively to each RNAi or paternal depletion line
#D4 and D6 refer to CD and HSD


# progress so far;
# done:
#   creation of PAC for sRNA from suptab
# fig 1A, B, C, fig 2 B, 
# 
# Where I need additional info;
# Fig 1E (need object Gene_target, could we fetch this info from elsewhere?)
# 
# half-done;
# Fig S1A - PAC object is done, needs to be reannotated

#loading of packages ###########################################################

library(seqpac)
library(stringr)
library(biomaRt)
library(circlize)
library(tidyr)
library(ggplot2)
library(tibble)
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(limma)
library(ggpubr)
library(Seurat)
library(SingleR)
library(KEGGREST)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(purrr)
library(viridis)
library(tximport)
library(FactoMineR)
library(factoextra)
library(pheatmap)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
library(GO.db)
library(readr)
library(rstatix)
library(readxl)
library(scmap)
library(SingleCellExperiment)
library(openxlsx)
library(SingleR)
library(scrapper)
library(RColorBrewer)
organism = "org.Dm.eg.db"
library(organism, character.only = TRUE)
library(patchwork)
library(VennDiagram)
library(igraph)
library(ggraph)
library(ggrepel)
library(scales)


#PREPROCESSING

############################## ##############################    sRNA seq    ############################## ##############################


#ADPATOR TRIMMING 
input <- "~/Merged"
output <- "~/Trimmed"

prog_report  <-  make_trim(
  input=input, output=output, threads=4, 
  adapt_3_set=c(type="hard_rm", min=10, mismatch=0.1), 
  adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCA", 
  polyG=c(type="hard_trim", min=20, mismatch=0.1),
  seq_range=c(min=14, max=70),
  quality=c(threshold=20, percent=0.8))


### Count sequences ##############################################################

path <- "~/Trimmed"
input <- list.files(path = path, pattern = "fastq.gz", all.files = FALSE,
                    full.names = TRUE)

counts  <- make_counts(input, threads=6, parse="default_neb", plot=TRUE, trimming = NULL,
                       evidence=c(experiment=2, sample=1))



### Pheno table ##################################################################
## Generate a Pheno table using the file names of test fastq 


Sample_ID <- colnames(counts$counts)
pheno_fl <- data.frame(row.names = Sample_ID,
                       Geno=rep(c("Aub", "GFP", "Piwi"), times=c(8,8,6)))

pheno_fl$diet <- ifelse(grepl("D4", Sample_ID), "D4", "D6")
table(pheno_fl$diet)


### Generate PAC object ##########################################################
## Generate PAC object

pac_raw <- make_PAC(pheno=pheno_fl, counts=counts$counts)

## If TRUE the PAC object is ok 
PAC_check(pac_raw)


### Filtering ####################################################################

pac_filt <- PAC_filter(pac_raw,nucleotide_range = c(15,55), norm="counts", 
                       threshold=1, coverage=90, stat = TRUE)


### Mapping ######################################################################
#Genome mapping 

outpath_genome <- "~/Genome"

ref_paths <- list(nuclear_genome = "~/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa",
                  Gene_map = "~/dmel-all-gene-r6.60.fa")

map_reanno(PAC = pac_filt, input = ref_paths, output = outpath_genome, 
           type = "internal", mismatches = 3, import = "genome", 
           threads = 6, keep_temp = TRUE,override = TRUE)


# Biotype mapping 
outpath_biotype <- "~/Reanno"

ref_paths <- list(piRNA = "~/piR_dme.fa",
                  miRNA ="~/miRBase_21-dme.fa",
                  tRNA = "~/dmel-all-tRNA-r6.60.fa",
                  rRNA ="~/drosophila_rRNA_all.fa",
                  sno_RNA = "~/Drosophila_Ensembl_snoRNA.fa",
                  snRNA = "~/snrna.fa",
                  CDS = "~/dmel-all-CDS-r6.60.fa",
                  five_prime = "~/dmel-all-five_prime_UTR-r6.60.fa",
                  three_prime = "~/dmel-all-three_prime_UTR-r6.60.fa",
                  all_genes = "~/dmel-all-gene-r6.60.fa",
                  all_introns = "~/dmel-all-intron-r6.60.fa",
                  all_pseudo = "~/dmel-all-pseudogene-r6.60.fa",
                  all_transcripts = "~/dmel-all-transcript-r6.60.fa",
                  all_TEs = "~/dmel-all-transposon-r6.60.fa",
                  lncRNA = "~/lncRNA_Central.fa")

map_reanno(pac_filt, input= ref_paths, output=outpath_biotype, 
           type="internal", mismatches=3,  import="biotype", 
           threads=6, keep_temp=TRUE,override = TRUE)



### Generate reanno objects ######################################################

### Generate a reanno object with make_reanno 
reanno_genome <- make_reanno(outpath_genome, PAC=pac_filt, mis_fasta_check = TRUE)
reanno_biotype <- make_reanno(outpath_biotype, PAC=pac_filt, mis_fasta_check = TRUE)

### Output merged with provided PAC object
pac_master <- add_reanno(reanno_genome, type="genome", 
                         genome_max=10, mismatches=3, merge_pac=pac_filt)

bio_search <- list(
  rRNA=c("5S", "5.8S", "12S", "16S", "18S", "28S", "pre_S45", "rRNA_other"),
  tRNA =c("_tRNA"),
  miRNA=c("dme-mir", "dme-mir-10" ),
  piRNA =c("piR"),
  CDS = c("type=CDS"),
  five_prime = c("type=five_prime_UTR"),
  all_genes = c("type=gene"),
  all_introns = c("type=intron"),
  all_pseudo = c("type=pseudogene"),
  three_prime = c("type=three_prime_UTR"),
  all_transcripts = c("type=mRNA"),
  sno_RNA = c("snoRNA"),
  snRNA = c("snRNA")
)


pac_master <- add_reanno(reanno_biotype, bio_search=bio_search, 
                         type="biotype", bio_perfect=FALSE, mismatches = 3,
                         merge_pac=pac_master)


hierarchy <- list(miRNA= "miRNA",
                  rRNA="rRNA",
                  tRNA="tRNA",
                  piRNA="piRNA" ,
                  sno_RNA = "snoRNA",
                  snRNA = "snRNA",
                  from_gene= "CDS|five_prime|all_genes|all_introns|all_pseudo|three_prime|all_transcripts" 
)


pac_master <- simplify_reanno(input=pac_master, hierarchy=hierarchy, mismatches=0, 
                              bio_name="Biotypes_mis0", merge_pac=TRUE)
pac_master <- simplify_reanno(input=pac_master, hierarchy=hierarchy, mismatches=1, 
                              bio_name="Biotypes_mis1", merge_pac=TRUE)
pac_master <- simplify_reanno(input=pac_master, hierarchy=hierarchy, mismatches=3, 
                              bio_name="Biotypes_mis3", merge_pac=TRUE)


pac_master <- PAC_filter(pac_master, 
                         anno_target= list("Any_genome" , "Hit"))


#Bind Geno and Diet
pac_master@Pheno <- pac_master@Pheno %>% 
  mutate(Geno_diet = paste(Geno, diet, sep = "_"))

### Norm + filtering + summary on mismatch 0 #####################################

pac_master <- PAC_norm(pac_master, norm= "cpm" , merge_pac = TRUE)
pac_master <- PAC_filter(pac_master, threshold=1, coverage=90, norm = "cpm",stat = F )
pac_master <- PAC_filter(pac_master, threshold=10,coverage=15, norm = "cpm",stat = F )


pac_master <- PAC_summary(pac_master, norm = "cpm", type = "means", 
                          pheno_target=list("Geno_diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "se", 
                          pheno_target=list("Geno_diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FC", 
                          pheno_target=list("Geno_diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FCgrand", 
                          pheno_target=list("Geno_diet"), merge_pac=TRUE)


pac_master <- PAC_summary(pac_master, norm = "cpm", type = "means", 
                          pheno_target=list("Geno"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "se", 
                          pheno_target=list("Geno"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FC", 
                          pheno_target=list("Geno"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FCgrand", 
                          pheno_target=list("Geno"), merge_pac=TRUE)



pac_master <- PAC_summary(pac_master, norm = "cpm", type = "means", 
                          pheno_target=list("diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "se", 
                          pheno_target=list("diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FC", 
                          pheno_target=list("diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FCgrand", 
                          pheno_target=list("diet"), merge_pac=TRUE)


pac_master@Anno$genome <- ifelse(
  grepl("^mito", pac_master@Anno$mis0_nuclear_genome) |
    grepl("^mito", pac_master@Anno$mis1_nuclear_genome) |
    grepl("^mito", pac_master@Anno$mis2_nuclear_genome) |
    grepl("^mito", pac_master@Anno$mis3_nuclear_genome),
  "mt", 
  "nu"
)

pac_mitos <- PAC_filter(pac_master, anno_target = list("genome", "mt"))
pac_nu <- PAC_filter(pac_master, anno_target = list("genome", "nu"))

Anno_df <- pac_mitos@Anno %>%
  rownames_to_column("sequence")

cols <- c("mis0_Gene_map_genome",
          "mis1_Gene_map_genome",
          "mis2_Gene_map_genome",
          "mis3_Gene_map_genome")

Anno_df <- Anno_df %>%
  rowwise() %>%
  mutate(
    FBgn_id = {
      vals <- c_across(all_of(cols))
      fbgn <- str_extract(vals, "FBgn\\d+")
      fbgn <- fbgn[!is.na(fbgn)]
      if (length(fbgn) > 0) fbgn[1] else NA_character_
    }
  ) %>%
  ungroup()


mart <- useEnsembl(
  biomart = "ensembl",
  dataset = "dmelanogaster_gene_ensembl"
)

mapping <- getBM(
  attributes = c("flybase_gene_id", "external_gene_name"),
  filters = "flybase_gene_id",
  values = unique(Anno_df$FBgn_id),
  mart = mart
)

Anno_df <- Anno_df %>%
  left_join(mapping, by = c("FBgn_id" = "flybase_gene_id")) %>%
  rename(Gene_target = external_gene_name)

pac_mitos@Anno <- Anno_df %>%
  column_to_rownames("sequence")

PAC_check(pac_mitos)


pac_master@Anno$Bio_all <- paste(
  pac_master@Anno$Biotypes_mis3,
  pac_master@Anno$genome,sep = "_")

lfc <- pac_master@summary$Log2FC_Geno_diet

swap_cols <- grep("_D4_vs_.*_D6$", colnames(lfc))
for (i in swap_cols) {
  old_name <- colnames(lfc)[i]
  parts <- strsplit(old_name, "_vs_")[[1]]
  new_name <- paste(parts[2], "vs", parts[1], sep = "_")
  lfc[, i] <- -lfc[, i]
  colnames(lfc)[i] <- new_name
}

pac_master@summary$Log2FC_Geno_diet <- lfc


pac_mitos@Anno <- pac_mitos@Anno  %>%
  mutate(
    mt_class = case_when(
      grepl("^mt:tRNA", Gene_target) ~ "tRNA",
      Gene_target %in% c("mt:srRNA", "mt:lrRNA") ~ "rRNA",
      TRUE ~ "protein_coding"
    )
    
  )




### Note: Use counts from sperm pre-processing (Do not use counts from supp table)
counts <- as.data.frame(pac_master@Counts)
CPM <- as.data.frame(pac_master@norm$cpm)
anno <- as.data.frame(pac_master@Anno)
pheno <- as.data.frame(pac_master@Pheno)
pheno$Geno <- gsub(" RNAi", "", pheno$Geno)
pheno$Geno_diet <- paste(pheno$Geno, pheno$diet, sep = "_")

pheno$Geno <- factor(pheno$Geno, levels=c( "GFP", "Piwi" , "Aub" ))
pheno$diet <- factor(pheno$diet, levels = c("D4", "D6"))


dds <- DESeq2::DESeqDataSetFromMatrix(
  counts,
  colData = pheno,
  design = ~ Geno
)

dds <- DESeq(dds)
resultsNames(dds)

res_piwi_sperm <- results(dds,  name = "Geno_Piwi_vs_GFP")
res_aub_sperm <- results(dds,  name = "Geno_Aub_vs_GFP")


dds <- DESeq2::DESeqDataSetFromMatrix(
  counts,
  colData = pheno,
  design = ~ Geno_diet
)

dds <- DESeq(dds)
resultsNames(dds)

res_gfp_diet_sperm <- results(dds, contrast = c("Geno_diet", "GFP_D6", "GFP_D4"))
res_piwi_diet_sperm <- results(dds, contrast = c("Geno_diet", "Piwi_D6", "Piwi_D4"))
res_aub_diet_sperm <- results(dds, contrast = c("Geno_diet", "Aub_D6", "Aub_D4"))



############################## ##############################    testis RNA prerocessing     ############################## ##############################



#Load Data
base_dir <- "~/quant_files_directory"
sample_folders <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

files <- sapply(sample_folders, function(folder) {
  list.files(folder, pattern = "quant.sf$", full.names = TRUE)
})

sample_names <- basename(dirname(files))
print(sample_names)


metadata <- data.frame(
  genotype = c(rep("Aub", 16),     
               rep("GFP", 16),     
               rep("Piwi", 15)),  
  
  diet = c(rep("D4", 8), rep("D6", 8),   
           rep("D4", 8), rep("D6", 8),   
           rep("D4", 7), rep("D6", 8)),  
  
  row.names = c("A_D4_1", "A_D4_2", "A_D4_3", "A_D4_4", "A_D4_5", "A_D4_6", "A_D4_7", "A_D4_8",
                "A_D6_1", "A_D6_2", "A_D6_3", "A_D6_4", "A_D6_5", "A_D6_6", "A_D6_7", "A_D6_8",
                "G_D4_1", "G_D4_2", "G_D4_3", "G_D4_4", "G_D4_5", "G_D4_6", "G_D4_7", "G_D4_8",
                "G_D6_1", "G_D6_2", "G_D6_3", "G_D6_4", "G_D6_5", "G_D6_6", "G_D6_7", "G_D6_8",
                "P_D4_2", "P_D4_3", "P_D4_4", "P_D4_5", "P_D4_6", "P_D4_7", "P_D4_8",
                "P_D6_1", "P_D6_2", "P_D6_3", "P_D6_4", "P_D6_5", "P_D6_6", "P_D6_7", "P_D6_8"),
  
  batch = c(rep("batch1", 15), rep("batch2", 16),rep("batch3", 16))
)

sample_names <- c("A_D4_1", "A_D4_2", "A_D4_3", "A_D4_4", "A_D4_5", "A_D4_6", "A_D4_7", "A_D4_8",
                  "A_D6_1", "A_D6_2", "A_D6_3", "A_D6_4", "A_D6_5", "A_D6_6", "A_D6_7", "A_D6_8",
                  "G_D4_1", "G_D4_2", "G_D4_3", "G_D4_4", "G_D4_5", "G_D4_6", "G_D4_7", "G_D4_8",
                  "G_D6_1", "G_D6_2", "G_D6_3", "G_D6_4", "G_D6_5", "G_D6_6", "G_D6_7", "G_D6_8",
                  "P_D4_2", "P_D4_3", "P_D4_4", "P_D4_5", "P_D4_6", "P_D4_7", "P_D4_8",
                  "P_D6_1", "P_D6_2", "P_D6_3", "P_D6_4", "P_D6_5", "P_D6_6", "P_D6_7", "P_D6_8")

batch <- character(length(sample_names))
sample_end <- as.numeric(sub(".*_(\\d)$", "\\1", sample_names))

batch[sample_end %in% c(1,2,3)] <- "batch1"
batch[sample_end %in% c(4,5)] <- "batch2"
batch[sample_end %in% c(6,7,8)] <- "batch3"
batch[sample_names %in% c("P_D4_3", "P_D6_3")] <- "batch2"
batch[sample_names %in% c( "A_D6_6","A_D4_6")] <- "batch2"

metadata <- data.frame(
  genotype = c(rep("Aub", 16), rep("GFP", 16), rep("Piwi", 15)),
  diet     = c(rep("D4", 8), rep("D6", 8),
               rep("D4", 8), rep("D6", 8),
               rep("D4", 7), rep("D6", 8)),
  batch    = batch,
  row.names = sample_names
)

metadata <- metadata %>%
  mutate(genotype_diet = paste(genotype, diet, sep = "_"))
names(files) <- rownames(metadata)


tx2gene <- as.data.frame(read_delim("~/tx2gene.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE))
txi.salm2 <- tximport(files, type="salmon", tx2gene = tx2gene,countsFromAbundance = "scaledTPM" )


counts <- as.data.frame(txi.salm2$counts)
counts$FBgn <- rownames(counts)

fbgns <- rownames(txi.salm2$counts)
gene_info <- getBM(attributes = c('flybase_gene_id', 'external_gene_name'),
                   filters = 'flybase_gene_id', 
                   values = fbgns, 
                   mart = ensembl)


counts <- counts %>%
  left_join(gene_info, by = c("FBgn" = "flybase_gene_id")) %>%
  mutate(external_gene_name = coalesce(external_gene_name, FBgn))

rownames(counts) <- counts$external_gene_name
colnames(counts)
counts <- counts[,-c(48,49)]
colnames(counts) <- rownames(metadata)
counts[] <- sapply(counts, as.integer)
rownames(txi.salm2$counts) <- rownames(counts)

tpm <- as.data.frame(txi.salm2$abundance)
tpm$FBgn <- rownames(tpm)
fbgns <- rownames(txi.salm2$abundance)
gene_info <- getBM(attributes = c('flybase_gene_id', 'external_gene_name'),
                   filters = 'flybase_gene_id', 
                   values = fbgns, 
                   mart = ensembl)

tpm <- tpm %>%
  left_join(gene_info, by = c("FBgn" = "flybase_gene_id")) %>%
  mutate(external_gene_name = coalesce(external_gene_name, FBgn))

rownames(tpm) <- tpm$external_gene_name
colnames(tpm)
tpm <- tpm[,-c(48,49)]
colnames(tpm) <- rownames(metadata)
tpm[] <- sapply(tpm, as.integer)

txi.salm2$abundance  <- tpm
txi.salm2$counts  <- counts

metadata$genotype <- factor(metadata$genotype, levels = c("GFP", "Aub", "Piwi"))
metadata$diet <- factor(metadata$diet, levels = c("D4", "D6"))


#Make pac and filter
pac <- make_PAC(counts = counts,pheno =  metadata)
pac_norm <- PAC_norm(pac, norm = "cpm")
pac_norm@norm$tpm <- as.data.frame(tpm)

pac_norm <- PAC_filter(pac_norm,norm= "tpm", threshold = 5, coverage = 80,stat = F) 
pac_testes <- pac_norm

keep <- rownames(pac_norm@Anno)
txi.salm2$abundance  <-  txi.salm2$abundance[keep,]
txi.salm2$counts     <-  txi.salm2$counts[keep,]

# Make sure sample order matches metadata
metadata <- metadata[colnames(txi.salm2$counts), , drop = FALSE]
stopifnot(all(colnames(txi.salm2$counts) == rownames(metadata)))

txi.salm2$counts <- as.matrix(txi.salm2$counts)
txi.salm2$abundance <- as.matrix(txi.salm2$abundance)

storage.mode(txi.salm2$counts) <- "numeric"
storage.mode(txi.salm2$abundance) <- "numeric"

# check
if ("length" %in% names(txi.salm2)) {
  txi.salm2$length <- as.matrix(txi.salm2$length)
  storage.mode(txi.salm2$length) <- "numeric"
}

TPM <- pac_norm@norm$tpm
batch <- metadata$batch

#batch correction
log2TPM <- log2(TPM + 1)
mod <- model.matrix(~ genotype * diet, data = metadata)
log2TPM_limma <- removeBatchEffect(log2TPM, batch = metadata$batch,design = mod)

#TPM to supplementary material
TPM_corrected <- 2^log2TPM_limma - 1


#from testes data
dds <- DESeqDataSetFromTximport(txi.salm2,
                                colData = metadata,
                                design = ~ batch + genotype_diet)



dsq <- DESeq(dds)
resultsNames(dsq)

res_gfp_diet_testes <- results(dsq, contrast = c("genotype_diet", "GFP_D6", "GFP_D4"))
res_Piwi_diet_testes <- results(dsq, contrast = c("genotype_diet", "Piwi_D6", "Piwi_D4"))
res_Aub_diet_testes <- results(dsq, contrast = c("genotype_diet", "Aub_D6", "Aub_D4"))


res_gfp_diet_testes <- as.data.frame(res_gfp_diet_testes)
res_Piwi_diet_testes <- as.data.frame(res_Piwi_diet_testes)
res_Aub_diet_testes <- as.data.frame(res_Aub_diet_testes)


############################## ##############################    RNAi embryo RNA preprocessing     ############################## ##############################



base_dir <- "~/quant_files_directory"
sample_folders <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

files <- sapply(sample_folders, function(folder) {
  list.files(folder, pattern = "quant.sf$", full.names = TRUE)
}, simplify = FALSE)

files <- files[lengths(files) > 0]

sample_names <- basename(names(files))
print(sample_names)

files <- unlist(files)
names(files) <- basename(dirname(files))


metadata <- data.frame(
  genotype = c(rep("Aub", 20), rep("GFP", 22), rep("Piwi", 22)),
  
  row.names = c(
    "A1_S7_F2","A10_S21_F2","A11_S22_F2","A12_S23_F2",
    "A13_S13_F1","A13_S24_F2","A14_S14_F1","A14_S25_F2",
    "A15_S15_F1","A16_S16_F1","A16_S26_F2","A17_S39_F2",
    "A18_S40_F2","A20_S41_F2","A21_S42_F2","A5_S8_F2",
    "A6_S9_F2","A7_S10_F2","A8_S11_F2","A9_S12_F2",
    "G1_S1_F2","G10_S10_F1","G10_S17_F2","G11_S11_F1",
    "G11_S18_F2","G12_S12_F1","G12_S19_F2","G13_S33_F2",
    "G14_S34_F2","G15_S20_F2","G16_S35_F2","G19_S36_F2",
    "G2_S2_F2","G20_S37_F2","G21_S38_F2","G3_S3_F2",
    "G4_S4_F2","G7_S7_F1","G8_S5_F2","G8_S8_F1",
    "G9_S6_F2","G9_S9_F1",
    "P1_S1_F1","P10_S27_F2","P11_S28_F2","P12_S29_F2",
    "P13_S30_F2","P14_S31_F2","P15_S32_F2","P17_S44_F2",
    "P18_S45_F2","P19_S46_F2","P2_S2_F1","P20_S47_F2",
    "P21_S48_F2","P3_S13_F2","P3_S3_F1","P4_S14_F2",
    "P4_S4_F1","P5_S15_F2","P5_S5_F1","P6_S43_F2",
    "P6_S6_F1","P8_S16_F2"
  ),
  
  batch = c(
    1,2,2,2,
    4,2,4,2,
    4,4,2,3,
    3,3,3,1,
    1,1,1,1,
    1,4,2,4,
    2,4,2,3,
    3,2,3,3,
    1,3,3,1,
    1,4,1,4,
    1,4,
    4,2,2,2,
    2,2,2,3,
    3,3,4,3,
    3,1,4,1,
    4,1,4,3,
    4,1
  )
)


tx2gene <- as.data.frame(read_delim("~/tx2gene.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE))
txi.salm <- tximport(files, type="salmon", tx2gene = tx2gene,countsFromAbundance = "scaledTPM")   


TPM <- as.data.frame(txi.salm$abundance)
TPM$FBgn <- rownames(TPM)
fbgns <- rownames(TPM)
gene_info <- getBM(attributes = c('flybase_gene_id', 'external_gene_name'),
                   filters = 'flybase_gene_id', 
                   values = fbgns, 
                   mart = ensembl)

TPM <- TPM %>%
  left_join(gene_info, by = c("FBgn" = "flybase_gene_id")) %>%
  mutate(external_gene_name = coalesce(external_gene_name, FBgn))

rownames(TPM) <- TPM$external_gene_name
colnames(TPM)
TPM <- TPM[,-c(65,66)]
colnames(TPM) <- rownames(metadata)
TPM[] <- sapply(TPM, as.integer)

counts <- as.data.frame(txi.salm$counts)
counts$FBgn <- rownames(counts)
fbgns <- rownames(counts)
gene_info <- getBM(attributes = c('flybase_gene_id', 'external_gene_name'),
                   filters = 'flybase_gene_id', 
                   values = fbgns, 
                   mart = ensembl)


counts <- counts %>%
  left_join(gene_info, by = c("FBgn" = "flybase_gene_id")) %>%
  mutate(external_gene_name = coalesce(external_gene_name, FBgn))

rownames(counts) <- counts$external_gene_name
colnames(counts)
counts <- counts[,-c(65,66)]
colnames(counts) <- rownames(metadata)
counts[] <- sapply(counts, as.integer)


rownames(txi.salm$counts) <- rownames(counts)
rownames(txi.salm$abundance) <- rownames(TPM)

log2TPM <- log2(TPM + 1)
mod <- model.matrix(~ genotype, data = metadata)
TPM_limma <- removeBatchEffect(log2TPM, batch = metadata$batch,design = mod)

#TPM to supplementary material
TPM_corrected <- 2^TPM_limma-1


pac <- make_PAC(counts = counts,pheno =  metadata)
pac_norm <- PAC_norm(pac, norm = "cpm")
pac_norm@norm$tpm <- as.data.frame(TPM)

pac_norm <- PAC_filter(pac_norm,norm= "tpm", threshold = 1, coverage = 90, stat = F) 

keep <- rownames(pac_norm@Anno)

txi.salm$abundance  <-  txi.salm$abundance[keep,]
txi.salm$counts     <-  txi.salm$counts[keep,]

# Make sure sample order matches metadata
metadata <- metadata[colnames(txi.salm$counts), , drop = FALSE]
stopifnot(all(colnames(txi.salm$counts) == rownames(metadata)))

txi.salm$counts <- as.matrix(txi.salm$counts)
txi.salm$abundance <- as.matrix(txi.salm$abundance)

storage.mode(txi.salm$counts) <- "numeric"
storage.mode(txi.salm$abundance) <- "numeric"

if ("length" %in% names(txi.salm)) {
  txi.salm$length <- as.matrix(txi.salm$length)
  storage.mode(txi.salm$length) <- "numeric"
}


metadata$genotype <- factor(metadata$genotype, levels=c( "GFP", "Piwi" , "Aub"))
metadata$batch <- factor(metadata$batch)

dds <- DESeqDataSetFromTximport(txi.salm,
                                colData = metadata,
                                design = ~ batch + genotype)

n_samples <- ncol(DESeq2::counts(dds))
keep <- (
  rowSums(DESeq2::counts(dds) >= 5) >= 0.70 * n_samples)


dds <- dds[keep, ]

dsq <- DESeq(dds)
resultsNames(dsq)

#deseq tables to Supplementary material
res_Piwi <- results(dsq, contrast = c("genotype","Piwi","GFP"))
res_Aub <- results(dsq, contrast = c("genotype","Aub","GFP"))




############################## ##############################    Diet embryo RNA preprocessing     ############################## ##############################

#Load Data
base_dir <- "~/quant_files_directory"
sample_folders <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

files <- sapply(sample_folders, function(folder) {
  list.files(folder, pattern = "quant.sf$", full.names = TRUE)
})

files
metadata <- data.frame(
  row.names = c("D4_1", "D4_10", "D4_11", "D4_12", "D4_13", 
                "D4_2", "D4_4", "D4_5", "D4_6", "D4_7", "D4_8", "D4_9",
                "D6_1", "D6_10", "D6_11", "D6_12", 
                "D6_2", "D6_3", "D6_4", "D6_5", "D6_6", "D6_7", "D6_8", "D6_9"),
  diet = c(rep("D4", 12), rep("D6", 12))
)

batch1_samples <- c("D6_1","D4_1","D6_2","D4_2","D6_3","D4_13",
                    "D6_4","D4_4","D6_5","D4_5","D6_6","D4_6")

metadata$batch <- ifelse(rownames(metadata) %in% batch1_samples, "batch1", "batch2")

library(readxl)

names(files) <- rownames(metadata)
tx2gene <- as.data.frame(read_delim("~/tx2gene.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE))
tx2gene_diet <- tx2gene

#generate Txi
txi.salm2 <- tximport(files, type="salmon", tx2gene = tx2gene_diet,countsFromAbundance = "scaledTPM")   

counts <- as.data.frame(txi.salm2$counts)
counts$FBgn <- rownames(counts)

fbgns <- rownames(txi.salm2$counts)
gene_info <- getBM(attributes = c('flybase_gene_id', 'external_gene_name'),
                   filters = 'flybase_gene_id', 
                   values = fbgns, 
                   mart = ensembl)


counts <- counts %>%
  left_join(gene_info, by = c("FBgn" = "flybase_gene_id")) %>%
  mutate(external_gene_name = coalesce(external_gene_name, FBgn))

rownames(counts) <- counts$external_gene_name
colnames(counts)
counts <- counts[,-c(25,26)]
colnames(counts) <- rownames(metadata)
counts[] <- sapply(counts, as.integer)
rownames(txi.salm2$counts) <- rownames(counts)

TPM <- as.data.frame(txi.salm2$abundance)
TPM$FBgn <- rownames(TPM)
fbgns <- rownames(TPM)
gene_info <- getBM(attributes = c('flybase_gene_id', 'external_gene_name'),
                   filters = 'flybase_gene_id', 
                   values = fbgns, 
                   mart = ensembl)


TPM <- TPM %>%
  left_join(gene_info, by = c("FBgn" = "flybase_gene_id")) %>%
  mutate(external_gene_name = coalesce(external_gene_name, FBgn))

rownames(TPM) <- TPM$external_gene_name
colnames(TPM)
TPM <- TPM[,-c(25,26)]
colnames(TPM) <- rownames(metadata)
TPM[] <- sapply(TPM, as.integer)

txi.salm2$abundance  <- TPM
txi.salm2$counts  <- counts



pac <- make_PAC(counts = counts,pheno =  metadata)
pac_norm <- PAC_norm(pac, norm = "cpm")
pac_norm@norm$tpm <- as.data.frame(TPM)


pac_norm <- PAC_filter(pac_norm,norm= "tpm", threshold = 1, coverage = 90, stat = F) 


keep <- rownames(pac_norm@Anno)

txi.salm2$abundance  <-  txi.salm2$abundance[keep,]
txi.salm2$counts     <-  txi.salm2$counts[keep,]

# Make sure sample order matches metadata
metadata <- metadata[colnames(txi.salm2$counts), , drop = FALSE]

stopifnot(all(colnames(txi.salm2$counts) == rownames(metadata)))

# Convert back to matrices
txi.salm2$counts <- as.matrix(txi.salm2$counts)
txi.salm2$abundance <- as.matrix(txi.salm2$abundance)

storage.mode(txi.salm2$counts) <- "numeric"
storage.mode(txi.salm2$abundance) <- "numeric"

# If length exists, it must also match
if ("length" %in% names(txi.salm2)) {
  txi.salm2$length <- as.matrix(txi.salm2$length)
  storage.mode(txi.salm2$length) <- "numeric"
}

# From mmc1-2 pdf (of Perez-Mojica et al.)
df <- read.table(
  "~/Perez-Mojica et al..txt",
  header = F,
  sep = "",
  stringsAsFactors = FALSE
)

colnames(df) <- c(
  "Replicate",
  "Pseudotime",
  "ID",
  "Sex"
)

head(df)

# Duplicated Replicates
df[duplicated(df$Replicate) | duplicated(df$Replicate, fromLast = TRUE), ]

# Duplicated IDs
df[duplicated(df$ID) | duplicated(df$ID, fromLast = TRUE), ]

metadata_pseudotime <- df[
  !duplicated(df$Replicate) &
    !duplicated(df$Replicate, fromLast = TRUE) &
    !duplicated(df$ID) &
    !duplicated(df$ID, fromLast = TRUE),
]

nrow(metadata_pseudotime)

# From SRA (of Perez-Mojica et al.)
embryos <- read.csv("~/GSE214118_Pseudotime_normalized_reads.csv")
rownames(embryos) <- embryos$X
embryos$X <- NULL

#filter out replicate samples from data
embryos <- embryos[, colnames(embryos) %in% metadata_pseudotime$ID]
ncol(embryos)

# metadata of reads, ID and pseudotime clustering
embryos_id <- metadata_pseudotime

embryos <- embryos[rowSums(embryos)>0,]
ann1 <- data.frame(cell_type1=embryos_id$Pseudotime,
                   row.names = colnames(embryos))
embryos <- data.frame(embryos)

## Transform Pérez-Mojica data to singlecellexperiment format
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts = as.matrix(embryos)), colData = ann1)
logcounts(sce) <- log2(normcounts(sce) + 1)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]

## Perform deconvolution
pred.grun <- SingleR(test=tpm, ref=sce, labels=sce$cell_type1, de.method="wilcox")
table(pred.grun$labels)
plotScoreHeatmap(pred.grun, show_colnames = TRUE)
pred.grun$scores

# df to Supplementary material
scores <- as.data.frame(pred.grun)
scores <- cbind(
  Sample_ID = rownames(scores),
  scores
)

#################################################
metadata$labels <- as.numeric(scores$labels)
batch <- metadata$batch
metadata$labels <- scale(metadata$labels)
mod <- model.matrix(~ diet, data = metadata)

metadata$diet <- factor(metadata$diet, levels = c("D4", "D6"))

metadata$batch <- factor(metadata$batch)
metadata$diet <- factor(metadata$diet)
metadata$labels <- as.numeric(scale(metadata$labels))

class(metadata$batch)
class(metadata$diet)
class(metadata$labels)

#Batch correction
log2TPM <- log2(TPM + 1)

log2TPM_limma <- removeBatchEffect(
  log2TPM,
  batch = metadata$batch,
  batch2 = metadata$labels,
  design = mod
)
TPM_corrected <- 2^log2TPM_limma - 1


#Deseq2
dds <- DESeqDataSetFromTximport(txi.salm2,
                                colData = metadata,
                                design = ~ batch + labels + diet)


dsq <- DESeq(dds)
resultsNames(dsq)

#df to Supplementary material
res_diet <- results(dsq, name = "diet_D6_vs_D4")


############################## ############################## ############################## ############################## ##############################
############################## ############################## Load data from Supp Material   ############################## ##############################
############################## ############################## ############################## ############################## ##############################

############################## ##############################         sRNA from sup.tab     ############################## ##############################



st <- as.data.frame(read_excel(path="C:/Users/alego91/OneDrive - Linköpings universitet/Anita Öst's files - 4_piwi story/TOPSECRET_piwiproject/Figures/Supplementary Material.xlsx", sheet="TableS1", skip=3))
rownames(st)<-st$...1
st <- st[,-1]

ano <- st[,c(23:34)]
st <- st[,-c(23:34)]

pheno <- data.frame(row.names= colnames(st), Sample_ID = colnames(st))
pheno$diet <- ifelse(grepl("D4", pheno$Sample_ID), "D4", "D6")
pheno$Geno <- ifelse(grepl("GFP", pheno$Sample_ID), "GFP", "Piwi")
pheno$Geno <- ifelse(grepl("Aub", pheno$Sample_ID), "Aub", pheno$Geno)
pheno$Geno_diet <- paste0(pheno$Geno,  sep="_", pheno$diet)

pac_master <- make_PAC(counts=st, pheno=pheno, anno=ano)
pac_master <- PAC_norm(pac_master)
pac_master@norm$cpm <- pac_master@Counts

pac_master <- PAC_summary(pac_master, norm = "cpm", type = "means", 
                          pheno_target=list("Geno_diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "se", 
                          pheno_target=list("Geno_diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FC", 
                          pheno_target=list("Geno_diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FCgrand", 
                          pheno_target=list("Geno_diet"), merge_pac=TRUE)


pac_master <- PAC_summary(pac_master, norm = "cpm", type = "means", 
                          pheno_target=list("Geno"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "se", 
                          pheno_target=list("Geno"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FC", 
                          pheno_target=list("Geno"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FCgrand", 
                          pheno_target=list("Geno"), merge_pac=TRUE)



pac_master <- PAC_summary(pac_master, norm = "cpm", type = "means", 
                          pheno_target=list("diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "se", 
                          pheno_target=list("diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FC", 
                          pheno_target=list("diet"), merge_pac=TRUE)
pac_master <- PAC_summary(pac_master, norm = "cpm", type = "log2FCgrand", 
                          pheno_target=list("diet"), merge_pac=TRUE)

pac_mitos <- PAC_filter(pac_master, anno_target = list("genome", "mt"))
pac_nu <- PAC_filter(pac_master, anno_target = list("genome", "nu"))



sperm_df <- pac_mitos@norm$cpm
metadata_sperm <- pac_master@Pheno


############################## ##############################         testis seq from sup.tab    ############################## ##############################



st <- as.data.frame(read_excel(path = "C:/Users/alego91/OneDrive - Linköpings universitet/Anita Öst's files - 4_piwi story/TOPSECRET_piwiproject/Figures/Supplementary Material 2.xls", sheet = "TableS3", skip = 3))
rownames(st)<-st$GeneID
st <- st[,-1]

pheno <- data.frame(row.names= colnames(st), Sample_ID = colnames(st))
pheno$diet <- ifelse(grepl("D4", pheno$Sample_ID), "D4", "D6")
pheno$Geno <- ifelse(grepl("G_", pheno$Sample_ID), "GFP", "Piwi")
pheno$Geno <- ifelse(grepl("A_", pheno$Sample_ID), "Aub", pheno$Geno)

colnames(pheno)[colnames(pheno) == "Geno"] <- "genotype"

pheno$genotype_diet <- paste(
  pheno$genotype,
  pheno$diet,
  sep = "_"
)

pac_master <- make_PAC(counts=st, pheno=pheno)
pac_master@norm$cpm <- pac_master@Counts

testes_df <- pac_master@norm$cpm
metadata_testes <- pac_master@Pheno



############################## ##############################         RNAi embryo seq from sup.tab    ############################## ##############################



st <- as.data.frame(read_excel(path = "C:/Users/alego91/OneDrive - Linköpings universitet/Anita Öst's files - 4_piwi story/TOPSECRET_piwiproject/Figures/Supplementary Material 2.xls", sheet = "TableS4", skip = 3))
rownames(st)<-st$GeneID
st <- st[,-1]

pheno <- data.frame(row.names= colnames(st), Sample_ID = colnames(st))
pheno$Geno <- ifelse(grepl("G", pheno$Sample_ID), "GFP", "Piwi")
pheno$Geno <- ifelse(grepl("A", pheno$Sample_ID), "Aub", pheno$Geno)

pac_master <- make_PAC(counts=st, pheno=pheno)
pac_master@norm$cpm <- pac_master@Counts

RNAi_embryo_df <- pac_master@norm$cpm
metadata_embryo_RNAi <- pac_master@Pheno
colnames(metadata_embryo_RNAi)[colnames(metadata_embryo_RNAi) == "Geno"] <- "genotype"


res_Piwi <- as.data.frame(read_excel(path = "C:/Users/alego91/OneDrive - Linköpings universitet/Anita Öst's files - 4_piwi story/TOPSECRET_piwiproject/Figures/Supplementary Material 2.xls", sheet = "TableS5", skip = 3))
res_Aub <- as.data.frame(read_excel(path = "C:/Users/alego91/OneDrive - Linköpings universitet/Anita Öst's files - 4_piwi story/TOPSECRET_piwiproject/Figures/Supplementary Material 2.xls", sheet = "TableS6", skip = 3))

rownames(res_Piwi)<-res_Piwi$GeneID
res_Piwi <- res_Piwi[,-1]

rownames(res_Aub)<-res_Aub$GeneID
res_Aub <- res_Aub[,-1]


############################## ##############################         Dietembryo seq from sup.tab    ############################## ##############################

st <- as.data.frame(read_excel(path = "C:/Users/alego91/OneDrive - Linköpings universitet/Anita Öst's files - 4_piwi story/TOPSECRET_piwiproject/Figures/Supplementary Material 2.xls", sheet = "TableS9", skip = 3))
rownames(st)<-st$GeneID
st <- st[,-1]

pheno <- data.frame(row.names= colnames(st), Sample_ID = colnames(st))
pheno$diet <- ifelse(grepl("D4", pheno$Sample_ID), "D4", "D6")

pac_master <- make_PAC(counts=st, pheno=pheno)
pac_master@norm$cpm <- pac_master@Counts

diet_embryo_df <- pac_master@norm$cpm
metadata_diet_embryo <- pac_master@Pheno

res_diet_embryo <- as.data.frame(read_excel(path = "C:/Users/alego91/OneDrive - Linköpings universitet/Anita Öst's files - 4_piwi story/TOPSECRET_piwiproject/Figures/Supplementary Material 2.xls", sheet = "TableS11", skip = 3))

rownames(res_diet_embryo)<-res_diet_embryo$GeneID
res_diet_embryo <- res_diet_embryo[,-1]



############################## ############################## ############################## ############################## ##############################
############################## ##############################    Generating all figures      ############################## ##############################
############################## ############################## ############################## ############################## ##############################



### FIG 1A #########################################

##  map pac to mito

pac_mitos <- PAC_summary(pac_mitos, norm="cpm")
#fasta file of mitochondrial genome NC_001709.1
covpi<-PAC_mapper(pac_mitos, input ="~/mito.fa",
                  mismatches = 3, threads=4)
covpirna<-PAC_covplot(pac_mitos, covpi, summary_target = list("cpmMeans_All"))

## clean data
df2 <- covpirna$NC_024511.2$data
df2$chr<-h
df2$end_column <- df2$Position
df2$start_column <- (as.numeric(df2$Position)-1)
df22 <- df2[,c(5,6,4,2,3)]
df22 <- df22[,c(3,2,1,5)]
colnames(df22)<-c("chr", "start", "end", "value1")
df22$start <- as.integer(df22$start)w
df22$end <- as.integer(df22$end)
df22$value1 <- log(df22$value1)
df22[df22$value1 %in% "-Inf",]$value1 <- 0
df22$value1 <- as.numeric(df22$value1)

df_chr = data.frame(
  name  = c("mt:tRNA:Ile-GAT",  "mt:tRNA:Gln-TTG",    "mt:tRNA:Met-CAT",
            "mt:ND2", "mt:tRNA:Trp-TCA", "mt:tRNA:Cys-GCA|mt:tRNA:Trp-TCA",
            "mt:tRNA:Tyr-GTA", "mt:CoI", "Mt:tRNA:Leu-TAA",
            "mt:CoII", "mt:tRNA:Lys-CTT", "mt:tRNA:Asp-GTC",
            "mt:ATP8", "mt:ATP6", "mt:CoIII",
            "mt:tRNA:Gly-TCC", "mt:ND3", "mt:tRNA:Ala-TGC",
            "mt:tRNA:Arg-TCG", "mt:tRNA:Asn-GTT", "mt:tRNA:Ser-GCT",
            "mt:tRNA:Glu-TTC", "mt:tRNA:Phe-GAA", "mt:ND5",
            "mt:tRNA:His-GTG", "mt:ND4", "mt:ND4L",
            "mt:tRNA:Thr-TGT", "mt:tRNA:Pro-TGG", "mt:ND6",
            "mt:Cyt-b", "mt:tRNA:Ser-TGA", "mt:ND1",
            "mt:tRNA:Leu-TAG", "mt:lrRNA", "mt:tRNA:Val-TAC",
            "mt:srRNA", "end"),
  start = c(1, 97, 171,
            240, 1264, 1322,
            1403, 1474, 3012,
            3083, 3768, 3840,
            3907, 4062, 4736,
            5543, 5608, 5983,
            6062, 6126, 6191,
            6259, 6344, 6409,
            8141, 8207, 9545,
            9838, 9904, 9971,
            10499, 11638, 11721,
            12670, 12735, 14058,
            14131, 14917),
  end   = c(65, 165, 239,
            1263, 1329, 1383,
            1468, 3009, 3077,
            3767, 3838, 3906,
            4068, 4736, 5524,
            5607, 5961, 6047,
            6125, 6190, 6258,
            6325, 6408, 8125,
            8206, 9545, 9835,
            9903, 9964, 10495,
            11635, 11703, 12659,
            12734, 14058, 14130,
            14916, 14918))
df_chr$name<-factor(df_chr$name, ordered=TRUE, levels = df_chr$name)
df24 <- df22[1:14918,]

##  plot circos
circos.clear()
circos.par("track.height"=0.50)
circos.genomicInitialize(df_chr,  
                         tickLabelsStartFromZero=FALSE)

circos.track(ylim = c(0, 1),
             bg.col = c("#f4dd91", "#f4dd91", "#f4dd91",
                        "#b1cecf",  "#f4dd91", "#f4dd91", 
                        "#f4dd91",  "#b1cecf",  "#f4dd91",
                        "#b1cecf",  "#f4dd91", "#f4dd91",
                        "#b1cecf","#b1cecf","#b1cecf",
                        "#f4dd91", "#b1cecf", "#f4dd91",
                        "#f4dd91","#f4dd91", "#f4dd91",
                        "#f4dd91", "#f4dd91","#b1cecf",
                        "#f4dd91","#b1cecf","#b1cecf",
                        "#f4dd91", "#f4dd91", "#b1cecf",
                        "#b1cecf","#f4dd91", "#b1cecf",
                        "#f4dd91","#81B185", "#f4dd91",
                        "#81B185","#b2b2b2"),
             bg.border = NA, track.height = 0.05)

circos.genomicTrack(df24, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, area=TRUE,
                                          type="l", col=c("#d09a9a"), lwd=0.5)
                    })



#################### FIG 1B ####################################################

pac_gfp <- PAC_filter(pac_mitos, pheno_target = list("Geno", "GFP"))
pac_gfp <- PAC_summary(pac_gfp, norm="cpm")
PAC_stackbar(pac_gfp, anno_target = list("mt_class"), norm = "cpm", summary_target = list("cpmMeans_All"))


#################### FIG 1C ####################################################
pac_pro <- PAC_filter(pac_mitos, anno_target = list("mt_class", "protein_coding"))
pac_tRNA <- PAC_filter(pac_mitos, anno_target = list("mt_class", "tRNA"))
pac_rRNA <- PAC_filter(pac_mitos, anno_target = list("mt_class", "rRNA"))


sizedist_plots_pro <- PAC_sizedist(pac_pro, 
                               anno_target=list("mt_class"),
                               pheno_target = list("Geno"),
                               summary_target=list("cpmMeans_Geno"),
                               ymax =3000)

sizedist_plots_tRNA <- PAC_sizedist(pac_tRNA, 
                               anno_target=list("mt_class"),
                               pheno_target = list("Geno"),
                               summary_target=list("cpmMeans_Geno"),
                               ymax =3000)

sizedist_plots_rRNA <- PAC_sizedist(pac_rRNA, 
                               anno_target=list("mt_class"),
                               pheno_target = list("Geno"),
                               summary_target=list("cpmMeans_Geno"),
                               ymax =3000)

PAC_pie(pac_pro, anno_target = list("Biotypes_mis3"),summary_target=list("cpmMeans_Geno") )
PAC_pie(pac_tRNA, anno_target = list("Biotypes_mis3"),summary_target=list("cpmMeans_Geno") )
PAC_pie(pac_rRNA, anno_target = list("Biotypes_mis3"),summary_target=list("cpmMeans_Geno") )


#################### FIG 1E and 1SA #############################################################

pac_mitos@Anno <- Anno_df %>%
  column_to_rownames("sequence")
cpm <- pac_mitos@norm$cpm
Anno <- pac_mitos@Anno
metadata <- pac_mitos@Pheno
metadata_sperm <- metadata

meta_gfp <- metadata %>%
  filter(Geno_diet == "GFP_D6")

samples <- intersect(colnames(cpm), rownames(meta_gfp))
cpm_gfp <- cpm[, samples, drop = FALSE]

cpm_gfp_df <- as.data.frame(cpm_gfp) %>%
  rownames_to_column("sequence") %>%
  left_join(
    Anno %>%
      rownames_to_column("sequence") %>%
      dplyr::select(sequence, Gene_target, Biotypes_mis3),
    by = "sequence"
  )

cpm_long <- cpm_gfp_df %>%
  pivot_longer(
    cols = -c(sequence, Gene_target, Biotypes_mis3),
    names_to = "sample",
    values_to = "cpm"
  ) %>%
  filter(
    !is.na(Gene_target),
    !is.na(Biotypes_mis3),
    !grepl("^mt:tRNA|^mt:srRNA|^mt:lrRNA", Gene_target) #(Figure 1SA)
  ) %>%
  group_by(Gene_target, Biotypes_mis3, sample) %>%
  summarise(sum_cpm = sum(cpm, na.rm = TRUE), .groups = "drop")

gene_order <- cpm_long %>%
  group_by(Gene_target) %>%
  summarise(total = sum(sum_cpm), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(Gene_target)

cpm_long$Gene_target <- factor(cpm_long$Gene_target, levels = gene_order)

plot_df <- cpm_long %>%
  group_by(Gene_target, Biotypes_mis3) %>%
  summarise(
    mean_cpm = mean(sum_cpm, na.rm = TRUE),
    sd_cpm = sd(sum_cpm, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  plot_df,
  aes(
    x = Gene_target,
    y = mean_cpm,
    fill = Biotypes_mis3
  )
) +
  geom_col(
    position = position_dodge(width = 0.9),
    width = 0.8
  ) +
  geom_errorbar(
    aes(
      ymin = mean_cpm - sd_cpm,
      ymax = mean_cpm + sd_cpm
    ),
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  theme_bw() +
  labs(
    x = "Gene target",
    y = "Sum CPM",
    fill = "Biotype (mis3)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


stats_res <- cpm_long %>%
  filter(Biotypes_mis3 %in% c("piRNA", "from_gene")) %>%
  
  group_by(Gene_target) %>%
  filter(n_distinct(Biotypes_mis3) == 2) %>%
  ungroup() %>%
  
  group_by(Gene_target) %>%
  summarise(
    test = list(wilcox.test(sum_cpm ~ Biotypes_mis3)),
    .groups = "drop"
  ) %>%
  mutate(tidy = lapply(test, broom::tidy)) %>%
  unnest(tidy) %>%
  dplyr::select(Gene_target, statistic, p.value, method)

stats_res

############ FIG S1 B ##########################################################

pac_mitos_piRNA <- PAC_filter(pac_mitos,nucleotide_range = c(23,29), anno_target = list("Biotypes_mis3", "piRNA"))

dd1 <- PAC_nbias(pac_mitos_piRNA,position = "1",norm = "cpm",
                 anno_target = list("Biotypes_mis3"), 
                 summary_target = list("cpmMeans_All"),data_only = T)

dd10 <- PAC_nbias(pac_mitos_piRNA,position = "10",norm = "cpm",
                  anno_target = list("Biotypes_mis3"), 
                  summary_target = list("cpmMeans_All"),data_only = T)

d1$cpm   <- g@summary$cpmMeans_Geno$`GFP RNAi`
d10$cpm  <- g@summary$cpmMeans_Geno$`GFP RNAi`

dd1$cpm  <- pac_mitos_piRNA@summary$cpmMeans_Geno$`GFP RNAi`
dd10$cpm <- pac_mitos_piRNA@summary$cpmMeans_Geno$`GFP RNAi`


make_plot_df <- function(df, dataset, bias_col){
  
  df %>%
    mutate(
      dataset = dataset,
      bias = .data[[bias_col]]
    ) %>%
    group_by(dataset, bias) %>%
    summarise(
      total_cpm = sum(cpm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      percent = total_cpm / sum(total_cpm) * 100
    )
}

plot_df <- bind_rows(
  make_plot_df(d1,   "fly_all_pos1",   "1_nuc_bias"),
  make_plot_df(d10,  "fly_all_pos10",  "10_nuc_bias"),
  
  make_plot_df(dd1,  "fly_piRNA_pos1",   "1_nuc_bias"),
  make_plot_df(dd10, "fly_piRNA_pos10",  "10_nuc_bias")
  
  
)
# set plotting order
plot_df$dataset <- factor(
  plot_df$dataset,
  levels = c(
    # fly first
    "fly_all_pos1",
    "fly_piRNA_pos1",
    "fly_all_pos10",
    "fly_piRNA_pos10")
)

ggplot(wplot_df,
       aes(x = dataset,
           y = percent,
           fill = bias)) +
  
  geom_col(width = 0.8) +
  
  scale_y_continuous(
    labels = scales::percent_format(scale = 1)
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    x = "",
    y = "Percent of CPM",
    fill = "Nucleotide bias"
  )


######################### FIG 2 B ##############################################


metadata_sperm <- pac_mitos@Pheno
sperm_df <- pac_mitos@norm$cpm

metadata_sperm_D4 <- metadata_sperm %>%
  filter(diet == "D4")

rownames(metadata_sperm_D4) <- metadata_sperm_D4$Sample_ID

common_samples <- intersect(colnames(sperm_df), rownames(metadata_sperm))

sperm_df_D4 <- sperm_df[, common_samples]

sample_sums <- colSums(sperm_df)

sum_df <- data.frame(
  sample = common_samples,
  genotype = metadata_sperm[common_samples, "Geno"],
  sum_expr = sample_sums[common_samples]
)

sum_df$genotype <- factor(
  sum_df$genotype,
  levels = c("GFP", "Piwi", "Aub")
)

anova_res <- aov(sum_expr ~ genotype, data = sum_df)

tukey_res <- TukeyHSD(anova_res)

tukey_table <- as.data.frame(tukey_res$genotype)

tukey_table$comparison <- rownames(tukey_table)
colnames(tukey_table) <- c("diff", "lwr", "upr", "p_adj", "comparison")

tukey_table <- tukey_table %>%
  separate(comparison, into = c("group1", "group2"), sep = "-") %>%
  mutate(
    stars = case_when(
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(group1, group2) %>%
  mutate(
    y.position = max(sum_df$sum_expr, na.rm = TRUE) *
      (1 + 0.1 * row_number())
  )


ggplot(sum_df, aes(x = genotype, y = sum_expr, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, linewidth = 1.5, color = "black") +
  geom_jitter(width = 0.05, alpha = 0.5, size = 3, color = "black") +
  scale_x_discrete(limits = c("GFP", "Piwi", "Aub")) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 2),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 40),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 14)
  ) +
  guides(fill = guide_legend(keywidth = 4, keyheight = 4)) +
  labs(x = "Genotype", y = "Summed CPM", fill = "Genotype") +
  stat_pvalue_manual(
    tukey_table,
    label = "stars",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    bracket.size = 1.5,
    size = 6
  )


######################### FIG 2 C ########################################################

testes_df <- testes_df[grep("^mt:", rownames(testes_df), ignore.case = FALSE), ]

metadata_testes_D4 <- metadata_testes %>%
  filter(diet == "D4")

common_samples <- intersect(colnames(testes_df), rownames(metadata_testes_D4))

testes_df_D4 <- testes_df[, common_samples]
testes_df_D4 <- as.data.frame(lapply(testes_df_D4, as.numeric))

sample_sums <- colSums(testes_df_D4)

sum_df <- data.frame(
  sample = common_samples,
  genotype = metadata_testes_D4[common_samples, "Geno"],
  sum_expr = sample_sums[common_samples]
)

sum_df$genotype <- factor(
  sum_df$genotype,
  levels = c("GFP", "Piwi", "Aub")
)

anova_res <- aov(sum_expr ~ genotype, data = sum_df)

summary(anova_res)

tukey_res <- TukeyHSD(anova_res)

tukey_table <- as.data.frame(tukey_res$genotype)

tukey_table$comparison <- rownames(tukey_table)

tukey_table <- tukey_table %>%
  separate(comparison, into = c("group1", "group2"), sep = "-") %>%
  rename(
    p_adj = `p adj`
  ) %>%
  mutate(
    stars = case_when(
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(group1, group2) %>%
  mutate(
    y.position = max(sum_df$sum_expr, na.rm = TRUE) *
      (1 + 0.1 * row_number())
  )

ggplot(sum_df, aes(x = genotype, y = sum_expr, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, linewidth = 1.5, color = "black") +
  geom_jitter(width = 0.05, alpha = 0.5, size = 3, color = "black") +
  scale_x_discrete(limits = c("GFP", "Piwi", "Aub")) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 2),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 40),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 14)
  ) +
  guides(fill = guide_legend(keywidth = 4, keyheight = 4)) +
  labs(x = "Genotype", y = "Summed CPM (mt genes)", fill = "Genotype") +
  stat_pvalue_manual(
    tukey_table,
    label = "stars",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    bracket.size = 1.5,
    size = 6
  )

######################### FIG 2 D ########################################################

nuc <- PAC_sizedist(pac_nu, 
                    anno_target=list("Biotypes_mis3"), 
                    pheno_target = list("Geno"), 
                    summary_target=list("cpmMeans_Geno"), 
                    ymax = 55000)

mito <- PAC_sizedist(pac_mitos, 
                     anno_target=list("Biotypes_mis3"), 
                     pheno_target = list("Geno"), 
                     summary_target=list("cpmMeans_Geno"), 
                     ymax = 2700)

cowplot::plot_grid(
  plotlist = c(nuc$Histograms, mito$Histograms),
  ncol = 2
)

######################### FIG 2 E ########################################################

pac_mitos_D4 <- PAC_filter(pac_mitos, pheno_target = list("diet", "D4"))


sperm_df <- pac_mitos_D4@norm$cpm %>%
  as.data.frame() %>%
  rownames_to_column("sRNA")

sperm_anno <- pac_mitos@Anno %>%
  rownames_to_column("sRNA") %>%
  dplyr::select(sRNA, gene = Gene_target)

testes_df <- testes_df[grep("^mt:", rownames(testes_df)), ] %>%
  as.data.frame() %>%
  rownames_to_column("gene")

colnames(metadata_testes)[colnames(metadata_testes) == "Geno"] <- "genotype"

testes_long <- testes_df %>%
  pivot_longer(-gene, names_to = "genotype", values_to = "expression") %>%
  mutate(genotype = metadata_testes$genotype[match(genotype, rownames(metadata_testes))])

testes_long <- testes_long %>%
  mutate(expression = as.numeric(expression))

sperm_long <- sperm_df %>%
  pivot_longer(-sRNA, names_to = "genotype", values_to = "expression") %>%
  left_join(sperm_anno, by = "sRNA") %>%
  filter(!is.na(gene)) %>%
  group_by(gene, genotype) %>%
  summarise(expression = sum(expression), .groups = "drop") %>%
  mutate(genotype = metadata_sperm$Geno[match(genotype, rownames(metadata_sperm))]) %>%
  mutate(genotype = gsub("\\s*RNAi\\s*", "", genotype))

combined_long <- bind_rows(
  mutate(testes_long, tissue = "testes"),
  mutate(sperm_long, tissue = "sperm")
)

plot_df <- combined_long %>%
  group_by(gene, genotype, tissue) %>%
  summarise(
    mean_expression = mean(expression),
    sd = sd(expression),
    sem = sd(expression) / sqrt(n()),
    .groups = "drop"
  ) %>%
  group_by(gene) %>%
  filter(all(c("testes", "sperm") %in% tissue)) %>%
  ungroup()

gene_order <- plot_df %>%
  filter(tissue == "sperm") %>%
  group_by(gene) %>%
  summarise(order_val = max(mean_expression), .groups = "drop") %>%
  arrange(desc(order_val)) %>%
  pull(gene)

plot_df$gene <- factor(plot_df$gene, levels = gene_order)

pd <- position_dodge(width = 0.9)

plot_df$genotype <- factor(
  plot_df$genotype,
  levels = c("GFP", "Piwi", "Aub")
)

plot_df$gene <- factor(
  plot_df$gene,
  levels = c("mt:ND4", "mt:ND1", "mt:ND2", "mt:CoI", "mt:CoIII", "mt:Cyt-b", "mt:lrRNA", "mt:srRNA")
)


ggplot(plot_df, aes(x = gene, y = mean_expression, fill = genotype)) +
  geom_col(position = pd, width = 0.85, color = "black") +
  geom_errorbar(
    aes(
      ymin = mean_expression - sd,
      ymax = mean_expression + sd
    ),
    position = pd,
    width = 0.2
  ) +
  facet_wrap(~tissue, scales = "free_y", ncol = 1) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 24),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  ) +
  labs(
    x = "Gene",
    y = "Mean expression",
    fill = "Genotype"
  )




######################### FIG 2 F ########################################################

wide_expr <- combined_long %>%
  group_by(gene, genotype, tissue) %>%
  summarise(mean_expression = mean(expression), .groups = "drop") %>%
  pivot_wider(names_from = tissue, values_from = mean_expression) %>%
  filter(!is.na(testes), !is.na(sperm)) %>%
  mutate(
    group = ifelse(grepl("^mt:ND", gene), "ND", "Rest"),
    label = ifelse(grepl("^mt:ND1|^mt:ND3|^mt:ND4|^mt:ND5|^mt:ND6", gene), gene, NA)
  )



ggplot(wide_expr, aes(x = testes, y = sperm)) +
  geom_point(
    data = subset(wide_expr, group == "Rest"),
    aes(color = genotype, shape = group),
    alpha = 0.7,
    size = 6
  ) +
  geom_point(
    data = subset(wide_expr, group == "ND"),
    aes(color = genotype, shape = group),
    size = 8
  ) +
  geom_smooth(
    data = subset(wide_expr, group == "ND"),
    aes(group = 1),
    method = "lm",
    se = FALSE,
    fullrange = T
  ) +
  geom_smooth(
    data = subset(wide_expr, group == "Rest"),
    aes(group = 1),
    method = "lm",
    se = FALSE,
    fullrange = F
  )+
  geom_text_repel(
    data = subset(wide_expr, !is.na(label)),
    aes(label = label),
    color = "black",
    size = 5,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.5, "lines")
  ) +
  labs(
    x = "Mean expression (testes)",
    y = "Mean expression (sperm)",
    color = "Genotype",
    shape = "Gene group"
  )+ ylim(c(0,10000))

cor_nd <- with(subset(wide_expr, group == "ND"),
               cor(testes, sperm, method = "pearson"))

cor_rest <- with(subset(wide_expr, group == "Rest"),
                 cor(testes, sperm, method = "pearson"))

cor_nd
cor_rest

test_nd <- with(subset(wide_expr, group == "ND"),
                cor.test(testes, sperm, method = "pearson"))
test_rest <- with(subset(wide_expr, group == "Rest"),
                  cor.test(testes, sperm, method = "pearson"))

test_nd
test_rest


######################### FIG S2 A #############################################

# Import Seurat object created by Gabby Vida and Steve DiNardo

testisv1  <- readRDS("~/FCAloomToSeurat2TFP_Annotations.rds")
gonia3 <- subset(testisv1, idents = c("25", "22", "5", "78", "40", "97", "41", "51", 
                                      "73", "35", "33", "48", "105", "21", "56", "100", 
                                      "3", "106", "16", "64", "45", "13", "55", "63",
                                      "39", "27", "83", "15", "0", "28", "26", "24", 
                                      "18", "50", "17", "96", "12", "30", "37", "75", "66", "10", "46"))
soma3 <- subset(testisv1, idents = c("62","36","58","81","77","65","98","47","104",
                                     "88","94","95","74","67","60","69","79","68","57","80","84","73"))

testisv1$cell_group <- "Other"

testisv1$cell_group[colnames(testisv1) %in% colnames(gonia3)] <- "Gonia3"
testisv1$cell_group[colnames(testisv1) %in% colnames(soma3)]  <- "Soma3"
Idents(testisv1) <- "cell_group"

DimPlot(
  testisv1,
  group.by = "cell_group",
  pt.size = 1.3,
  cols = c(
    "Gonia3" = "#ef8060ff",
    "Soma3"  = "#ecc745ff",
    "Other"  = "grey80"
  )
)


FeaturePlot(testisv1, features = "piwi", pt.size = 1.3)
FeaturePlot(testisv1, features = "aub", pt.size = 1.3)


######################### FIG S2 G #############################################

#  Comparison of bulk testis RNA to scRNAseq from Raz et al.
#  using SingleR



##  load the Seurat object from Raz et al. and transform
FCA_bioHub_ScrapeB <- readRDS(file = "~/FCAloomToSeurat2TFP_Annotations.rds")
scFCA <- as.SingleCellExperiment(FCA_bioHub_ScrapeB)

## using TPM-normalized reads

TPM_testes

pred.grun <- SingleR(test=TPM_testes, ref=scFCA, labels=scFCA$S_annotation, de.method="wilcox")
table(pred.grun$labels)
plotScoreHeatmap(pred.grun, show_colnames = TRUE)
plotDeltaDistribution(pred.grun, ncol = 3)
plotMarkerHeatmap(pred.grun, TPM_testes, 
                  label="early-mid elongation-stage spermatid")




######################### FIG S2 I #############################################


ggplot(plot_df, aes(x = gene, y = mean_expression, fill = tissue)) +
  geom_col(position = pd, width = 0.85, color = "black") +
  geom_errorbar(
    aes(ymin = mean_expression - sd,
        ymax = mean_expression + sd),
    position = pd,
    width = 0.2
  ) +
  facet_wrap(~genotype, scales = "free_y", ncol = 1) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 24),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  ) +
  labs(
    x = "Gene",
    y = "Mean expression",
    fill = "Genotype"
  )


######################### FIG 3C ###############################################




P <- EnhancedVolcano(
  res_Piwi,
  lab = rownames(res_Piwi),
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.1,
  FCcutoff = 0.5,
  labSize = 6,
  colAlpha = 0.75,
  legendPosition = "bottom",
  title = "Volcano Plot Piwi",
  pointSize = 4,
  col = c("grey50", "grey50", "darkslategray3", "dodgerblue3"),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  xlim = c(-3, 3),
  ylim = c(0, 15)
)



A <- EnhancedVolcano(
  res_Aub,
  lab = rownames(res_Aub),
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.1,
  FCcutoff = 0.5,
  labSize = 6,
  colAlpha = 0.75,
  legendPosition = "bottom",
  title = "Volcano Plot Aub",
  pointSize = 4,
  col = c("grey50", "grey50", "darkslategray3", "dodgerblue3"),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  xlim = c(-3, 3),
  ylim = c(0, 15)
)


P | A


######################### FIG 3D ###############################################


res_Piwi <- res_Piwi[!is.na(res_Piwi$log2FoldChange) & !is.na(res_Piwi$padj), ]
res_Aub <- res_Aub[!is.na(res_Aub$log2FoldChange) & !is.na(res_Aub$padj), ]

res_Piwi <- as.data.frame(res_Piwi)
res_Aub <- as.data.frame(res_Aub)

set.seed(123)

ranks_Piwi <- res_Piwi$stat
names(ranks_Piwi) <- rownames(res_Piwi)

ranks_Aub <- res_Aub$stat
names(ranks_Aub) <- rownames(res_Aub)

ranks_Piwi <- sort(ranks_Piwi, decreasing = TRUE)
ranks_Aub  <- sort(ranks_Aub, decreasing = TRUE)

gsea_Piwi <- gseGO(
  geneList = ranks_Piwi,
  OrgDb = org.Dm.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  minGSSize = 20,
  maxGSSize = 200,
  pvalueCutoff = 0.05,
  verbose = FALSE,
  pAdjustMethod = "BH",
  nPermSimple = 10000
)

gsea_Aub <- gseGO(
  geneList = ranks_Aub,
  OrgDb = org.Dm.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  minGSSize = 20,
  maxGSSize = 200,
  pvalueCutoff = 0.05,
  verbose = FALSE,
  pAdjustMethod = "BH",
  nPermSimple = 10000
)

gsea_Piwi_s <- clusterProfiler::simplify(
  gsea_Piwi,
  cutoff = 0.6,
  by = "setSize",
  select_fun = min,
  measure = "Wang"
)

gsea_Aub_s <- clusterProfiler::simplify(
  gsea_Aub,
  cutoff = 0.6,
  by = "setSize",
  select_fun = min,
  measure = "Wang"
)


gsea_Piwi_df <- as.data.frame(gsea_Piwi_s) %>%
  mutate(comparison = "Piwi vs GFP")

gsea_Aub_df <- as.data.frame(gsea_Aub_s) %>%
  mutate(comparison = "Aub vs GFP")

gsea_combined <- bind_rows(
  gsea_Piwi_df,
  gsea_Aub_df
) %>%
  mutate(
    direction = ifelse(NES > 0, "Up", "Down"),
    Description = factor(Description)
  )

top_terms <- gsea_combined %>%
  arrange(p.adjust, desc(abs(NES))) %>%   #arrange check how it works
  group_by(comparison) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  pull(Description) %>%
  unique()


gsea_Piwi_df <- as.data.frame(gsea_Piwi) %>%
  mutate(comparison = "Piwi vs GFP")

gsea_Aub_df <- as.data.frame(gsea_Aub) %>%
  mutate(comparison = "Aub vs GFP")

gsea_co <- bind_rows(
  gsea_Piwi_df,
  gsea_Aub_df
) %>%
  mutate(
    direction = ifelse(NES > 0, "Up", "Down"),
    Description = factor(Description)
  )


plot_df <- gsea_co %>%
  filter(Description %in% top_terms) %>%
  mutate(
    Description = reorder(Description, NES)
  )


go_term_map <- c(
  
  # Metabolism / mitochondrial respiration
  "NADH metabolic process" =
    "Metabolism / mitochondrial respiration",
  
  "NAD metabolic process" =
    "Metabolism / mitochondrial respiration",
  
  "respiratory chain complex I" =
    "Metabolism / mitochondrial respiration",
  
  "mitochondrial electron transport, NADH to ubiquinone" =
    "Metabolism / mitochondrial respiration",
  
  "oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor" =
    "Metabolism / mitochondrial respiration",
  
  "tricarboxylic acid cycle" =
    "Metabolism / mitochondrial respiration",
  
  "lipid catabolic process" =
    "Metabolism / mitochondrial respiration",
  
  "fatty acid catabolic process" =
    "Metabolism / mitochondrial respiration",
  
  "small molecule biosynthetic process" =
    "Metabolism / mitochondrial respiration",
  
  "small molecule catabolic process" =
    "Metabolism / mitochondrial respiration",
  
  
  # Translation / ribosome
  "structural constituent of ribosome" =
    "Translation / ribosome",
  "mitochondrial small ribosomal subunit" =
    "Translation / ribosome",
  "mitochondrial gene expression" =
    "Translation / ribosome",
  "mitochondrial translation" =
    "Translation / ribosome",
  
  # Transport
  "ATPase complex" =
    "Transport",
  "transporter complex" =
    "Transport",
  
  # DNA maintenance
  "DNA recombination" =
    "DNA maintenance",
  "catalytic activity, acting on DNA" =
    "DNA maintenance",
  "ATP-dependent activity, acting on DNA" =
    "DNA maintenance",
  
  # Other
  "ligand-gated monoatomic ion channel activity" =
    "Other"
)

plot_df <- gsea_co %>%
  mutate(
    Description_chr = as.character(Description)
  ) %>%
  filter(Description_chr %in% names(go_term_map)) %>%
  mutate(
    GO_group = unname(go_term_map[Description_chr]),
    
    GO_group = factor(
      GO_group,
      levels = unique(go_term_map)
    ),
    
    direction = factor(direction, levels = c("Up", "Down")),
    
    Description = factor(
      Description_chr,
      levels = rev(names(go_term_map))
    )
  )


ggplot(
  plot_df,
  aes(
    x = NES,
    y = Description,
    color = comparison,
    size = -log10(p.adjust)
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey70"
  ) +
  geom_point(alpha = 0.75) +
  scale_size_continuous(range = c(6, 20)) +
  
  facet_grid(
    rows = vars(GO_group),
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_color_manual(
    values = c(
      "Piwi vs GFP" = "#c8951aff",
      "Aub vs GFP" = "#9c1b53ee"
    )
  ) +
  labs(
    x = "Normalized enrichment score, NES",
    y = NULL,
    color = NULL,
    size = "-log10(adj. p)",
    title = "GO GSEA grouped by functional category",
    subtitle = "Piwi vs GFP and Aub vs GFP"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.line.x = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    axis.line.y = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    axis.ticks.x = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    axis.ticks.y = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    
    axis.text.y = element_text(size = 14),
    strip.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )




######################### FIG 3E ###############################################



aerobic_genes <- plot_df %>%
  dplyr::filter(Description == "respiratory chain complex I") %>%
  dplyr::pull(core_enrichment) %>%
  strsplit("/") %>%
  unlist() %>%
  unique()

aerobic_genes_present <- intersect(
  aerobic_genes,
  rownames(RNAi_embryo_df)
)

aerobic_mat <- RNAi_embryo_df[
  aerobic_genes_present,
  ,
  drop = FALSE
]
aerobic_mat <- as.data.frame(lapply(aerobic_mat, as.numeric))
rownames(aerobic_mat) <- aerobic_genes_present

metadata_embryo_RNAi <- metadata_embryo_RNAi[colnames(aerobic_mat), , drop = FALSE]

metadata_embryo_RNAi$genotype <- factor(
  metadata_embryo_RNAi$genotype,
  levels = c("GFP", "Piwi", "Aub")
)

TPM_avg <- sapply(
  levels(metadata_embryo_RNAi$genotype),
  function(g) {
    samples <- rownames(metadata_embryo_RNAi)[metadata_embryo_RNAi$genotype == g]
    rowMeans(aerobic_mat[, samples, drop = FALSE])
  }
)

TPM_avg_log <- log2(TPM_avg + 1)

annotation_col <- data.frame(
  genotype = factor(
    colnames(TPM_avg_log),
    levels = c("GFP", "Piwi", "Aub")
  )
)

rownames(annotation_col) <- colnames(TPM_avg_log)

my_colors <- colorRampPalette(
  c(
    "#9B83B3",
    "#E8DDEB",
    "#F7F5F2",
    "#F2DFC0",
    "#E4B66A"
  )
)(100)

breaks <- seq(-1, 1, length.out = 101)

pheatmap(
  TPM_avg_log,show_rownames = T,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = my_colors,
  #breaks = breaks,
  border_color = "white"
)


######################### FIG 3G ###############################################

RNAi_embryo_df_no_mt <- RNAi_embryo_df[-grep("^mt:"  , rownames(RNAi_embryo_df), ignore.case = FALSE), ]

RNAi_embryo_df_no_mt[] <- lapply(
  RNAi_embryo_df_no_mt,
  as.numeric
)


modules <- list(
  
  N = c(
    "ND-51", "ND-24", "NdufV3", "ND-75", "ND-18",
    "ND-13A", "ND-B8", "ND-B14", "ND-B14.5A", "ND-B17.2"
  ),
  
  Q = c(
    "ND-49", "ND-30", "ND-20", "ND-23",
    "ND-13B", "ND-39"
  ),
  
  P = c(
    "ND-MWFE", "ND-19", "ND-42", "ND-B16.6", "NdufA3",
    "ND-B14.7", "ND-15", "ND-B14.5B",
    "mt:ND1", "mt:ND2", "mt:ND3", "mt:ND6", "mt:ND4L",
    "ND-MNLL", "NP15.6", "ND-B15", "ND-B22", "ND-PDSW",
    "ND-SGDH", "ND-B17", "ND-ASHI", 
    "ND-B18", "ND-B12", "mt:ND4", "mt:ND5"
  )
)



for (module_name in names(modules)) {
  
  cat("\n====================================\n")
  cat("MODULE:", module_name, "\n")
  cat("====================================\n")
  
  genes_of_interest <- modules[[module_name]]
  
  genes_present <- intersect(
    genes_of_interest,
    rownames(RNAi_embryo_df_no_mt)
  )
  
  if (length(genes_present) == 0) {
    cat("No genes found for this module.\n")
    next
  }
  
  module_expression <- colSums(
    RNAi_embryo_df_no_mt[genes_present, , drop = FALSE]
  )
  
  plot_df <- data.frame(
    sample = colnames(RNAi_embryo_df_no_mt),
    genotype = metadata_embryo_RNAi$genotype,
    expression_sum = module_expression
  )
  
  plot_df$genotype <- factor(
    plot_df$genotype,
    levels = c("GFP", "Piwi", "Aub")
  )
  normality_results <- plot_df %>%
    group_by(genotype) %>%
    shapiro_test(expression_sum)
  
  print(normality_results)
  
  all_normal <- all(normality_results$p > 0.05)
  
  if (all_normal) {
    
    cat("\nUsing ANOVA + TukeyHSD\n")
    
    anova_model <- aov(
      expression_sum ~ genotype,
      data = plot_df
    )
    
    print(summary(anova_model))
    
    pairwise_stats <- plot_df %>%
      tukey_hsd(expression_sum ~ genotype)
    
    print(pairwise_stats)
    
  } else {
    
    cat("\nUsing Kruskal-Wallis + Dunn test\n")
    
    kruskal_results <- kruskal_test(
      plot_df,
      expression_sum ~ genotype
    )
    
    print(kruskal_results)
    
    pairwise_stats <- dunn_test(
      plot_df,
      expression_sum ~ genotype,
      p.adjust.method = "BH"
    )
    
    print(pairwise_stats)
  }

  y_max <- max(plot_df$expression_sum, na.rm = TRUE)
  y_min <- min(plot_df$expression_sum, na.rm = TRUE)
  y_range <- y_max - y_min
  
  if (y_range == 0) {
    y_range <- ifelse(y_max == 0, 1, abs(y_max) * 0.02)
  }
  
  pairwise_stats$y.position <- y_max + seq(
    0.15 * y_range,
    0.45 * y_range,
    length.out = nrow(pairwise_stats)
  )
  
  p <- ggplot(
    plot_df,
    aes(
      x = genotype,
      y = expression_sum,
      fill = genotype
    )
  ) +
    
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8
    ) +
    stat_pvalue_manual(
      pairwise_stats,
      label = "p.adj.signif",
      hide.ns = TRUE
    ) +
    
    labs(
      title = paste0(module_name, " Module Expression"),
      x = "Genotype",
      y = "Summed TPM"
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.10))
    ) +
    
    coord_cartesian(clip = "off") +
    
    theme_bw(base_size = 14) +
    
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(10, 20, 20, 20)
    )
  
  print(p)
}



######################### FIG S3 A #############################################

mt_genes <- grep(
  "^mt:",
  rownames(RNAi_embryo_df),
  value = TRUE,
  ignore.case = FALSE
)

RNAi_embryo_df_mt <- RNAi_embryo_df[mt_genes, , drop = FALSE]

RNAi_embryo_df_mt[] <- lapply(
  RNAi_embryo_df_mt,
  as.numeric
)

df_long <- as.data.frame(RNAi_embryo_df_mt) |>
  rownames_to_column("gene") |>
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression"
  ) |>
  left_join(
    data.frame(
      sample = rownames(metadata_embryo_RNAi),
      genotype = metadata_embryo_RNAi$genotype
    ),
    by = "sample"
  ) |>
  mutate(
    gene = as.character(gene),
    expression = as.numeric(expression),
    genotype = factor(
      genotype,
      levels = c("GFP", "Piwi", "Aub")
    )
  )

df_summary <- df_long %>%
  group_by(gene, genotype) %>%
  summarise(
    mean_expression = mean(expression, na.rm = TRUE),
    sd_expression = sd(expression, na.rm = TRUE),
    n = sum(!is.na(expression)),
    sem = sd_expression / sqrt(n),
    .groups = "drop"
  )

gene_order <- df_summary %>%
  filter(genotype == "GFP") %>%
  arrange(desc(mean_expression)) %>%
  pull(gene)

df_long <- df_long %>%
  mutate(
    gene = factor(
      gene,
      levels = gene_order
    )
  )

normality <- df_long %>%
  group_by(gene, genotype) %>%
  filter(sum(!is.na(expression)) >= 3) %>%
  shapiro_test(expression)

gene_normality <- normality %>%
  group_by(gene) %>%
  summarise(
    all_normal = all(p > 0.05),
    .groups = "drop"
  )

pairwise_results <- list()

for (g in levels(df_long$gene)) {
  
  dat <- df_long %>%
    filter(gene == g)
  
  is_normal <- gene_normality %>%
    filter(gene == g) %>%
    pull(all_normal)
  
  if (length(is_normal) == 0) {
    next
  }
  
  if (is_normal) {
    
    pairwise <- dat %>%
      tukey_hsd(expression ~ genotype) %>%
      mutate(
        method = "ANOVA + Tukey"
      )
    
  } else {
    
    pairwise <- dat %>%
      dunn_test(
        expression ~ genotype,
        p.adjust.method = "BH"
      ) %>%
      mutate(
        method = "Kruskal-Wallis + Dunn"
      )
  }
  
  pairwise_results[[as.character(g)]] <- pairwise
}

pairwise_results <- bind_rows(
  pairwise_results,
  .id = "gene"
) %>%
  mutate(
    gene = as.character(gene),
    significance = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001  ~ "***",
      p.adj <= 0.01   ~ "**",
      p.adj <= 0.05   ~ "*",
      TRUE            ~ "ns"
    )
  )

#To supplementary
#write_xlsx(
#  pairwise_results,
#  "pairwise_results.xlsx"
#)

y_positions <- df_long %>%
  mutate(
    gene = as.character(gene)
  ) %>%
  group_by(gene) %>%
  summarise(
    ymax = max(expression, na.rm = TRUE),
    ymin = min(expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    range = ymax - ymin,
    step = case_when(
      range > 0 ~ range * 0.18,
      ymax > 0 ~ ymax * 0.15,
      TRUE ~ 0.2
    )
  )

pairwise_sig <- pairwise_results %>%
  filter(p.adj <= 0.05) %>%
  left_join(
    y_positions,
    by = "gene"
  ) %>%
  group_by(gene) %>%
  arrange(
    factor(group1, levels = c("GFP", "Piwi", "Aub")),
    factor(group2, levels = c("GFP", "Piwi", "Aub"))
  ) %>%
  mutate(
    y.position = ymax + row_number() * step
  ) %>%
  ungroup() %>%
  mutate(
    gene = factor(
      gene,
      levels = gene_order
    ),
    xmin = case_when(
      group1 == "GFP"  ~ 1,
      group1 == "Piwi" ~ 2,
      group1 == "Aub"  ~ 3
    ),
    xmax = case_when(
      group2 == "GFP"  ~ 1,
      group2 == "Piwi" ~ 2,
      group2 == "Aub"  ~ 3
    )
  )

brackets_horizontal <- pairwise_sig %>%
  transmute(
    gene,
    xmin,
    xmax,
    y = y.position
  )

brackets_vertical <- bind_rows(
  pairwise_sig %>%
    transmute(
      gene,
      x = xmin,
      y = y.position,
      yend = y.position - step * 0.25
    ),
  pairwise_sig %>%
    transmute(
      gene,
      x = xmax,
      y = y.position,
      yend = y.position - step * 0.25
    )
)

star_positions <- pairwise_sig %>%
  mutate(
    x = (xmin + xmax) / 2,
    y = y.position + step * 0.12
  )

ggplot(
  df_long,
  aes(
    x = genotype,
    y = expression,
    fill = genotype
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.5,
    alpha = 0.7
  ) +
  geom_segment(
    data = brackets_horizontal,
    aes(
      x = xmin,
      xend = xmax,
      y = y,
      yend = y
    ),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_segment(
    data = brackets_vertical,
    aes(
      x = x,
      xend = x,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_text(
    data = star_positions,
    aes(
      x = x,
      y = y,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 5
  ) +
  facet_wrap(
    ~ gene,
    scales = "free_y",
    ncol = 6
  ) +
  scale_x_discrete(
    limits = c("GFP", "Piwi", "Aub")
  ) +
  theme_bw() +
  labs(
    x = "Genotype",
    y = "Expression (TPM)",
    fill = "Genotype"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

######################### FIG 4A ########################################################


x <- pac_mitos

x@norm$Log2CPM <- log2(x@norm$cpm + 1)

PCA <- PAC_pca(
  PAC = x,
  norm = "Log2CPM",
  style = "pheno",
  graphs = FALSE,
  pheno_target = list("Geno_diet"),
  anno_target = list("genome")
)

target_colors <- c("red", "blue", "orange")
pheno_target <- list("Geno_diet")
col <- x@Pheno[, pheno_target[[1]]]
col <- as.factor(col)

pca_df <- data.frame(
  PCA1 = PCA$ind$coord[,1],
  PCA2 = PCA$ind$coord[,2],
  group = col
)

hulls <- pca_df %>%
  group_by(group) %>%
  dplyr::slice(chull(PCA1, PCA2)) %>%
  ungroup()

n_groups <- length(levels(col))
target_colors <- RColorBrewer::brewer.pal(max(3, n_groups), "Set1")

p <- fviz_pca_ind(
  PCA,
  geom = "point",
  axes = c(1, 2),
  pointsize = 4,
  ggtheme = theme_classic(),
  title = "PCA on Mitos",
  label = "none"
)

p + 
  geom_polygon(data = hulls, 
               aes(x = PCA1, y = PCA2, fill = group, group = group),
               alpha = 0.25, color = NA) +
  geom_point(data = pca_df,
             aes(x = PCA1, y = PCA2, color = group),
             size = 6, shape = 16, alpha = 0.8) +
  scale_fill_manual(values = target_colors) +
  scale_color_manual(values = target_colors) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  labs(x = paste0("PC1 (", round(PCA$eig[1,2], 1), "%)"),
       y = paste0("PC2 (", round(PCA$eig[2,2], 1), "%)"),
       fill = "Group",
       color = "Group")


######################### FIG 4B ###############################################

testes_df_mt <- testes_df[grep("^mt:", rownames(testes_df), ignore.case = FALSE), ]
testes_df_mt[] <- lapply(
  testes_df_mt,
  as.numeric
)

# PCA
res.pca <- PCA(t(log2(testes_df_mt + 1)), graph = FALSE)
pca_df <- as.data.frame(res.pca$ind$coord)
pca_df$group <- as.factor(metadata_testes$genotype_diet)

hulls <- pca_df %>%
  group_by(group) %>%
  slice(chull(Dim.1, Dim.2)) %>%
  ungroup()

n_groups <- length(levels(pca_df$group))
target_colors <- brewer.pal(max(3, min(n_groups, 9)), "Set1")

ggplot(pca_df, aes(x = Dim.1, y = Dim.2)) +
  
  geom_polygon(
    data = hulls,
    aes(fill = group, group = group),
    alpha = 0.20,
    color = NA
  ) +
  
  geom_point(
    aes(color = group),
    size = 4,
    alpha = 0.9
  ) +
  
  scale_fill_manual(values = target_colors) +
  scale_color_manual(values = target_colors) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  
  labs(
    x = paste0("PC1 (", round(res.pca$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(res.pca$eig[2, 2], 1), "%)"),
    color = "Group",
    fill = "Group"
  )





######################### FIG 4C ########################################################
#mitos

lfc <- pac_mitos@summary$Log2FC_Geno_diet

cols_switch <- c(
  "GFP_D4_vs_GFP_D6",
  "Piwi_D4_vs_Piwi_D6",
  "Aub_D4_vs_Aub_D6")

# Reverse the LFC
lfc[, cols_switch] <- -lfc[, cols_switch]

# Rename the columns
colnames(lfc)[match(cols_switch, colnames(lfc))] <- c(
  "GFP_D6_vs_GFP_D4",
  "Piwi_D6_vs_Piwi_D4",
  "Aub_D6_vs_Aub_D4")


pac_mitos@summary$Log2FC_Geno_diet <- lfc


ji <- PAC_jitter(pac_mitos, summary_target = list("Log2FC_Geno_diet"), anno_target = list("genome"),style = "violin",limits = c(-2, 2))

colnames(pac_mitos@summary$Log2FC_Geno_diet)

cowplot::plot_grid(
  plotlist = c(ji$`GFP_D6_vs_GFP_D4`,
               ji$`Piwi_D6_vs_Piwi_D4`,
               ji$`Aub_D6_vs_Aub_D4`),
  ncol = 1
)

#from processing testes data

mat_testes <- cbind(
  Piwi_diet = res_Piwi_diet[rownames(res_Piwi_diet), "log2FC"],
  GFP_diet  = res_gfp_diet[rownames(res_Piwi_diet), "log2FC"],
  Aub_diet  = res_Aub_diet[rownames(res_Piwi_diet), "log2FC"]
)

rownames(mat_testes) <- rownames(res_Piwi_diet)

mat_testes <- mat_testes[
  grepl("^mt:", rownames(mat_testes)),,drop = FALSE]
plot_data <- as.data.frame(mat_testes) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(
    cols = c(Piwi_diet, GFP_diet, Aub_diet),
    names_to = "genotype",
    values_to = "log2FC"
  )

ggplot(
  plot_data,
  aes(x = genotype, y = log2FC, fill = genotype)
) +
  geom_hline(
    yintercept = 0,
    color = "#707177",
    linewidth = 0.6
  ) +
  coord_cartesian(ylim = c(-2, 2)) +
  geom_violin(
    width = 0.9,
    trim = TRUE,
    scale = "width",
    color = "black"
  ) +
  geom_boxplot(
    width = 0.25,
    fill = "white",
    color = "black",
    alpha = 0.7,
    outlier.shape = NA
  ) +
  labs(
    title = "Mitochondrial genes",
    x = "Genotype",
    y = "log2FC (D6 vs D4)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 13)
  )


######################### FIG 4D ########################################################

#sperm sRNA
CPM <- pac_master@summary$cpmMeans_Geno_diet
CPM_sample <- pac_master@norm$cpm

metadata <- pac_master@Pheno
anno <- pac_master@Anno


plot_df <- as.data.frame(CPM) %>%
  rownames_to_column("sequence") %>%
  pivot_longer(
    -sequence,
    names_to = c("Geno", "diet"),
    names_sep = "_(?=D)",   # split before D4/D6
    values_to = "CPM"
  ) %>%
  pivot_wider(
    names_from = diet,
    values_from = CPM,
    names_prefix = "Mean_"
  ) %>%
  left_join(
    anno %>%
      rownames_to_column("sequence") %>%
      dplyr::select(sequence, genome),
    by = "sequence"
  ) %>%
  filter(Mean_D4 > 0, Mean_D6 > 0)


plot_df$genome <- factor(plot_df$genome, levels = c("nu", "mt"))
ggplot(plot_df,
       aes(x = Mean_D4,
           y = Mean_D6,
           colour = genome)) +
  geom_point(
    data = subset(plot_df, genome == "nu"),
    alpha = 0.5,
    size = 1.2
  ) +
  geom_point(
    data = subset(plot_df, genome == "mt"),
    alpha = 0.8,
    size = 1.2
  ) +
  scale_color_manual(
    values = c(
      "nu" = "grey",
      "mt" = "darkred"
    ),
    name = "Genome origin"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~Geno, nrow = 1) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()







######################### FIG 4E ########################################################

res_gfp_diet_sperm$genome <- anno$genome[
  match(rownames(res_gfp_diet_sperm), rownames(anno))
]

res_piwi_sperm$genome <- anno$genome[
  match(rownames(res_piwi_sperm), rownames(anno))
]

res_aub_sperm$genome <- anno$genome[
  match(rownames(res_aub_sperm), rownames(anno))
]


mt_hsd <- subset(res_gfp_diet_sperm,genome=="mt")
mt_piwi <- subset(res_piwi_sperm,genome=="mt")
mt_aub <- subset(res_aub_sperm,genome=="mt")

common_mt <- Reduce(intersect,list(
  
  rownames(mt_hsd),
  rownames(mt_piwi),
  rownames(mt_aub)
  
))

lfc_matrix <- data.frame(
  
  HSD=mt_hsd[common_mt,"log2FoldChange"],
  Piwi=mt_piwi[common_mt,"log2FoldChange"],
  Aub=mt_aub[common_mt,"log2FoldChange"]
  
)

rownames(lfc_matrix) <- common_mt


lfc_matrix$Aub_neg <- lfc_matrix$Aub < 0

base_theme <- theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    aspect.ratio = 1
  )

cols <- c("FALSE" = "grey80",
          "TRUE"  = "#B22222")

p1 <- ggplot(lfc_matrix, aes(HSD, Piwi)) +
  
  geom_point(
    data = subset(lfc_matrix, !Aub_neg),
    colour = "grey80",
    alpha = 0.5,
    size = 2
  ) +
  
  geom_point(
    data = subset(lfc_matrix, Aub_neg),
    aes(colour = Aub_neg),
    size = 2.8
  ) +
  
  geom_smooth(
    method = "lm",
    colour = "#0072B2",
    fill = "#B3D7F0",
    linewidth = 1
  ) +
  
  stat_cor(
    method = "pearson",
    label.x.npc = 0.02,
    label.y.npc = 0.98,
    size = 5
  ) +
  
  scale_colour_manual(values = cols,
                      labels = c("Aub ≥ 0", "Aub < 0")) +
  
  labs(
    title = "A",
    x = expression(HSD~log[2]~Fold~Change),
    y = expression(Piwi~log[2]~Fold~Change)
  ) +
  
  base_theme

p2 <- ggplot(lfc_matrix, aes(HSD, Aub)) +
  
  geom_point(
    data = subset(lfc_matrix, !Aub_neg),
    colour = "grey80",
    alpha = 0.5,
    size = 2
  ) +
  
  geom_point(
    data = subset(lfc_matrix, Aub_neg),
    aes(colour = Aub_neg),
    size = 2.8
  ) +
  
  geom_hline(yintercept = 0,
             linetype = 2,
             colour = "grey60") +
  
  geom_vline(xintercept = 0,
             linetype = 2,
             colour = "grey60") +
  
  geom_smooth(
    method = "lm",
    colour = "#0072B2",
    fill = "#B3D7F0",
    linewidth = 1
  ) +
  
  stat_cor(
    method = "pearson",
    label.x.npc = 0.02,
    label.y.npc = 0.98,
    size = 5
  ) +
  
  scale_colour_manual(values = cols,
                      labels = c("Aub ≥ 0", "Aub < 0")) +
  
  labs(
    title = "B",
    x = expression(HSD~log[2]~Fold~Change),
    y = expression(Aub~log[2]~Fold~Change)
  ) +
  
  base_theme

p3 <- ggplot(lfc_matrix, aes(Piwi, Aub)) +
  
  geom_point(
    data = subset(lfc_matrix, !Aub_neg),
    colour = "grey80",
    alpha = 0.5,
    size = 2
  ) +
  
  geom_point(
    data = subset(lfc_matrix, Aub_neg),
    aes(colour = Aub_neg),
    size = 2.8
  ) +
  
  geom_hline(yintercept = 0,
             linetype = 2,
             colour = "grey60") +
  
  geom_vline(xintercept = 0,
             linetype = 2,
             colour = "grey60") +
  
  geom_smooth(
    method = "lm",
    colour = "#0072B2",
    fill = "#B3D7F0",
    linewidth = 1
  ) +
  
  stat_cor(
    method = "pearson",
    label.x.npc = 0.02,
    label.y.npc = 0.98,
    size = 5
  ) +
  
  scale_colour_manual(values = cols,
                      labels = c("Aub ≥ 0", "Aub < 0")) +
  
  labs(
    title = "C",
    x = expression(Piwi~log[2]~Fold~Change),
    y = expression(Aub~log[2]~Fold~Change)
  ) +
  
  base_theme

final_plot <- (p3 + p1 + p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

final_plot


######################### FIG S4 A & B ########################################################

loadings <- as.data.frame(PCA$var$coord)
loadings$feature <- rownames(loadings)

top100_PC1 <- loadings %>%
  mutate(abs_PC1 = abs(Dim.1)) %>%
  arrange(desc(abs_PC1)) %>%
  slice_head(n = 100) %>%
  dplyr::select(feature, PC1 = Dim.1, abs_PC1)

top100_PC2 <- loadings %>%
  mutate(abs_PC2 = abs(Dim.2)) %>%
  arrange(desc(abs_PC2)) %>%
  slice_head(n = 100) %>%
  dplyr::select(feature, PC2 = Dim.2, abs_PC2)

top100_PC1
top100_PC2



x@Anno$PC1_top100 <- rownames(x@Anno) %in% rownames(top100_PC1)
x@Anno$PC2_top100 <- rownames(x@Anno) %in% rownames(top100_PC2)

sum(x@Anno$PC1_top100)
sum(x@Anno$PC2_top100)
sum(x@Anno$PC1_top100 & x@Anno$PC2_top100)


x@Anno$Biotype_2 <- ifelse(
  x@Anno$Biotypes_mis3 == "piRNA",
  "piRNA",
  x@Anno$mt_class
)

pac_pc1 <- PAC_filter(x, anno_target = list("PC1_top100", TRUE))
pac_pc2 <- PAC_filter(x, anno_target = list("PC2_top100", TRUE))



mean_CPM <- pac_master@summary$cpmMeans_Geno_diet

GFP_D4_CPM <- mean_CPM[, "GFP RNAi_D4"]





x@Anno$GFP_D4_CPM <- GFP_D4_CPM[
  match(rownames(x@Anno), rownames(mean_CPM))
]

df_pc1 <- x@Anno %>%
  filter(PC1_top100) %>%
  group_by(Gene_target, Biotype_2) %>%
  summarise(
    CPM = sum(GFP_D4_CPM, na.rm = TRUE),
    .groups = "drop"
  )


df_pc2 <- x@Anno %>%
  filter(PC2_top100) %>%
  group_by(Gene_target, Biotype_2) %>%
  summarise(
    CPM = sum(GFP_D4_CPM, na.rm = TRUE),
    .groups = "drop"
  )


p_pc1 <- ggplot(
  df_pc1,
  aes(x = Gene_target, y = CPM, fill = Biotype_2)
) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 14) +
  labs(
    title = "Top 100 PC1 drivers",
    x = "Gene target",
    y = "GFP D4 CPM",
    fill = "Biotype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


p_pc2 <- ggplot(
  df_pc2,
  aes(x = Gene_target, y = CPM, fill = Biotype_2)
) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 14) +
  labs(
    title = "Top 100 PC2 drivers",
    x = "Gene target",
    y = "GFP D4 CPM",
    fill = "Biotype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_pc1
p_pc2


######################### FIG S4 C, D and E ########################################################

res_piwi_plot <- res_piwi_sperm[order(res_piwi_sperm$genome == "mt"), ]

keyvals_piwi <- ifelse(
  res_piwi_plot$genome == "mt",
  "Mitochondrial",
  "Nuclear"
)

names(keyvals_piwi) <- rownames(res_piwi_plot)

EnhancedVolcano(
  res_piwi_plot,
  lab = rep("", nrow(res_piwi_plot)),
  x = "log2FoldChange",
  y = "padj",
  selectLab = NULL,
  colCustom = c(
    "Mitochondrial" = "red",
    "Nuclear" = "grey70"
  )[keyvals_piwi],
  pCutoff = 0.1,
  FCcutoff = 0.5,
  labSize = 0,
  colAlpha = 0.75,
  legendPosition = "bottom",
  title = "Piwi vs GFP",
  pointSize = 3,
  drawConnectors = FALSE,
  ylim = c(0, 9),
  xlim = c(-4.5, 4.5)
)

res_aub_plot <- res_aub_sperm[order(res_aub_sperm$genome == "mt"), ]

keyvals_aub <- ifelse(
  res_aub_plot$genome == "mt",
  "Mitochondrial",
  "Nuclear"
)

names(keyvals_aub) <- rownames(res_aub_plot)

EnhancedVolcano(
  res_aub_plot,
  lab = rep("", nrow(res_aub_plot)),
  x = "log2FoldChange",
  y = "padj",
  selectLab = NULL,
  colCustom = c(
    "Mitochondrial" = "red",
    "Nuclear" = "grey70"
  )[keyvals_aub],
  pCutoff = 0.1,
  FCcutoff = 0.5,
  labSize = 0,
  colAlpha = 0.75,
  legendPosition = "bottom",
  title = "aub vs GFP",
  pointSize = 3,
  drawConnectors = FALSE,
  ylim = c(0, 9),
  xlim = c(-4.5, 4.5)
)


res_gfp_plot <- res_gfp_diet_sperm[order(res_gfp_diet_sperm$genome == "mt"), ]

keyvals_gfp <- ifelse(
  res_gfp_plot$genome == "mt",
  "Mitochondrial",
  "Nuclear"
)

names(keyvals_gfp) <- rownames(res_gfp_plot)

EnhancedVolcano(
  res_gfp_plot,
  lab = rep("", nrow(res_gfp_plot)),
  x = "log2FoldChange",
  y = "padj",
  selectLab = NULL,
  colCustom = c(
    "Mitochondrial" = "red",
    "Nuclear" = "grey70"
  )[keyvals_gfp],
  pCutoff = 0.1,
  FCcutoff = 0.5,
  labSize = 0,
  colAlpha = 0.75,
  legendPosition = "bottom",
  title = "gfp d6 vs GFP d4",
  pointSize = 3,
  drawConnectors = FALSE,
  ylim = c(0, 9),
  xlim = c(-4.5, 4.5)
)



######################### FIG S 4F ########################################################

Aub_neg_genes <- rownames(lfc_matrix)[lfc_matrix$Aub_neg]
top100_PC1
top100_PC2



venn.plot <- venn.diagram(
  x = list(
    Orange = Aub_neg_genes,
    PC1 = rownames(top100_PC1),
    PC2 = rownames(top100_PC2)
  ),
  filename = NULL,
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20, 0),
  cat.dist = c(0.05, 0.05, 0.05),
  margin = 0.1
)


venn.plot



######################### FIG S 4G ########################################################

#from testes data preprocessing

P <- EnhancedVolcano(res_Piwi_diet_testes,
                     lab = rownames(res_Piwi_diet_testes),
                     x = 'log2FoldChange',
                     y = 'padj',
                     pCutoff = 0.1,
                     FCcutoff = 0.5,
                     labSize = 4,
                     colAlpha = 0.75,
                     legendPosition = "bottom",
                     title = "Piwi",
                     pointSize = 3,
                     col = c("grey50", "grey50", "darkblue", "darkgreen"),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5)




A <- EnhancedVolcano(res_Aub_diet_testes,
                     lab = rownames(res_Aub_diet_testes),
                     x = 'log2FoldChange',
                     y = 'padj',
                     pCutoff = 0.1,
                     FCcutoff = 0.5,
                     labSize = 4,
                     colAlpha = 0.75,legendPosition = "none",
                     title = "Aub",
                     pointSize = 3,
                     col = c("grey50", "grey50", "darkblue", "darkgreen"),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5)

G <- EnhancedVolcano(res_gfp_diet_testes,
                     lab = rownames(res_gfp_diet_testes),
                     x = 'log2FoldChange',
                     y = 'padj',
                     pCutoff = 0.1,
                     FCcutoff = 0.5,
                     labSize = 4,
                     colAlpha = 0.75,
                     legendPosition = "none",
                     title = "GFP",
                     pointSize = 3,
                     col = c("grey50", "grey50", "darkblue", "darkgreen"),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5)



G
P
A




### FIG 5B ########################################################


res_diet_embryo <- res_diet_embryo[!is.na(res_diet_embryo$padj), ]


deg_genes <- res_diet_embryo[
  res_diet_embryo$padj < 0.1 & abs(res_diet_embryo$log2FoldChange) > 0.5,
]

deg_gene_names <- rownames(deg_genes)





keyvals <- ifelse(
  res_diet_embryo$log2FoldChange < -1 & 
    res_diet_embryo$padj < 0.1, 'darkblue',
  ifelse(res_diet_embryo$log2FoldChange > 1  & 
           res_diet_embryo$padj < 0.1, 'orange',
         'grey70'))
keyvals[is.na(keyvals)] <- 'grey70'
names(keyvals)[keyvals == 'orange'] <- 'Upregulated'
names(keyvals)[keyvals == 'grey70'] <- 'mid'
names(keyvals)[keyvals == 'darkblue'] <- 'Downregulated'



EnhancedVolcano(res_diet_embryo,
                lab = rownames(res_diet_embryo),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = rownames(res_diet_embryo)[which(names(keyvals) %in% c('Upregulated', 'Downregulated'))],
                colCustom = keyvals,
                pCutoff = 0.1,
                FCcutoff = 1,
                labSize = 4,
                colAlpha = 0.75,
                legendPosition = "bottom",
                title = "Volcano Plot diet",
                pointSize = 3,
                drawConnectors = TRUE,
                widthConnectors = 0.5, 
                ylim = c(0, 4),
                xlim=c(-4.5, 4.5))



######################### FIG 5C ########################################################

high_D6 <- rownames(res_diet_embryo[res_diet_embryo$log2FoldChange > 1 & 
                               res_diet_embryo$padj < 0.1, ])

low_D6  <- rownames(res_diet_embryo[res_diet_embryo$log2FoldChange < -1 & 
                                      res_diet_embryo$padj < 0.1, ])



gse_down <- enrichGO(
  gene          = low_D6 ,
  universe = rownames(res_diet),
  OrgDb         = organism,
  keyType       = "SYMBOL",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.2
)

gse_up <- enrichGO(
  gene          = high_D6 ,
  universe = rownames(res_diet),
  OrgDb         = organism,
  keyType       = "SYMBOL",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.2
)

up <- barplot(gse_up, showCategory = 5) +
  scale_fill_gradient(
    low = "orange",
    high = "grey50"
  ) +
  ggtitle("Upregulated by D6")+
  theme(axis.text.y = element_text(size = 20, hjust = 1),
        axis.text.x = element_text(size = 20, hjust = 1),
        plot.title = element_text(size = 22, face = "bold")
  )

down <- barplot(gse_down, showCategory = 5) +
  scale_fill_gradient(
    low = "darkblue",
    high = "grey50"
  ) +
  ggtitle("Downregulated by D6")+
  theme(axis.text.y = element_text(size = 20, hjust = 1),
        axis.text.x = element_text(size = 20, hjust = 1),
        plot.title = element_text(size = 22, face = "bold")
  )


p <- cowplot::plot_grid(
  up,
  down,
  ncol = 1
)

p


######################### FIG 5D ########################################################


diet_embryo_df_no_mt <- diet_embryo_df[-grep("^mt:"  , rownames(diet_embryo_df), ignore.case = FALSE), ]


modules <- list(
  
  N = c(
    "ND-51", "ND-24", "NdufV3", "ND-75", "ND-18",
    "ND-13A", "ND-B8", "ND-B14", "ND-B14.5A", "ND-B17.2"
  ),
  
  Q = c(
    "ND-49", "ND-30", "ND-20", "ND-23",
    "ND-13B", "ND-39"
  ),
  
  P = c(
    "ND-MWFE", "ND-19", "ND-42", "ND-B16.6", "NdufA3",
    "ND-B14.7", "ND-15", "ND-B14.5B",
    "mt:ND1", "mt:ND2", "mt:ND3", "mt:ND6", "mt:ND4L",
    "ND-MNLL", "NP15.6", "ND-B15", "ND-B22", "ND-PDSW",
    "ND-SGDH", "ND-B17", "ND-ASHI", 
    "ND-B18", "ND-B12", "mt:ND4", "mt:ND5"
  )
)



for (module_name in names(modules)) {
  
  cat("\n====================================\n")
  cat("MODULE:", module_name, "\n")
  cat("====================================\n")
  
  genes_of_interest <- modules[[module_name]]
  
  genes_present <- intersect(
    genes_of_interest,
    rownames(diet_embryo_df_no_mt)
  )
  
  if (length(genes_present) == 0) {
    cat("No genes found for this module.\n")
    next
  }
  
  module_expression <- colSums(
    diet_embryo_df_no_mt[genes_present, , drop = FALSE]
  )
  
  plot_df <- data.frame(
    sample = colnames(diet_embryo_df_no_mt),
    diet = metadata_diet_embryo$diet,
    expression_sum = module_expression
  )
  
  plot_df$genotype <- factor(
    plot_df$diet,
    levels = c("D4", "D6")
  )
  normality_results <- plot_df %>%
    group_by(diet) %>%
    shapiro_test(expression_sum)
  
  print(normality_results)
  
  all_normal <- all(normality_results$p > 0.05)
  
  if (all_normal) {
    
    cat("\nUsing ANOVA + TukeyHSD\n")
    
    anova_model <- aov(
      expression_sum ~ diet,
      data = plot_df
    )
    
    print(summary(anova_model))
    
    pairwise_stats <- plot_df %>%
      tukey_hsd(expression_sum ~ diet)
    
    print(pairwise_stats)
    
  } else {
    
    cat("\nUsing Kruskal-Wallis + Dunn test\n")
    
    kruskal_results <- kruskal_test(
      plot_df,
      expression_sum ~ diet
    )
    
    print(kruskal_results)
    
    pairwise_stats <- dunn_test(
      plot_df,
      expression_sum ~ diet,
      p.adjust.method = "BH"
    )
    
    print(pairwise_stats)
  }
  y_max <- max(plot_df$expression_sum, na.rm = TRUE)
  y_min <- min(plot_df$expression_sum, na.rm = TRUE)
  y_range <- y_max - y_min
  
  if (y_range == 0) {
    y_range <- ifelse(y_max == 0, 1, abs(y_max) * 0.02)
  }
  
  pairwise_stats$y.position <- y_max + seq(
    0.15 * y_range,
    0.45 * y_range,
    length.out = nrow(pairwise_stats)
  )
  
  p <- ggplot(
    plot_df,
    aes(
      x = diet,
      y = expression_sum,
      fill = genotype
    )
  ) +
    
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8
    ) +
    stat_pvalue_manual(
      pairwise_stats,
      label = "p.adj.signif",
      hide.ns = TRUE
    ) +
    
    labs(
      title = paste0(module_name, " Module Expression"),
      x = "Genotype",
      y = "Summed TPM"
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.10))
    ) +
    
    coord_cartesian(clip = "off") +
    
    theme_bw(base_size = 14) +
    
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(10, 20, 20, 20)
    )
  
  print(p)
}
