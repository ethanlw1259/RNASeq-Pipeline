library(tximport)
library(DESeq2)
library(data.table)
library(GenomicFeatures)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(AnnotationDbi)
library(dplyr)
library(tidyr)

#Load file names
metaFile="/scratch/eleiter-weintraub/rnaSeqProj/dataFiles/samples/KBaseWithGroups.tsv"
salmonFilePath = "/scratch/eleiter-weintraub/rnaSeqProj/dataFiles/salmonResults"
starFilePath = "/scratch/eleiter-weintraub/rnaSeqProj/dataFiles/starResults"
salmonFiles = list.files(path = salmonFilePath, pattern = "quant.sf", full.names = TRUE, recursive = TRUE)
starFiles = list.files(path = starFilePath, pattern = "ReadsPerGene.out.tab", full.names = TRUE, recursive = TRUE)
names(salmonFiles) <- list.files(salmonFilePath)
txdbFile = "/scratch/eleiter-weintraub/rnaSeqProj/dataFiles/refGenomes/gencode_v40.lncipedia_v5_2_hc.annotation.sqlite"
#gtfFile = "/scratch/eleiter-weintraub/rnaSeqProj/dataFiles/refGenomes/gencode_v40.lncipedia_v5_2_hc.annotation.gtf"
rdsFile = "/home/eleiter-weintraub/R/salmonTxImport.rds"
biomartPath = '/scratch/eleiter-weintraub/rnaSeqProj/dataFiles/mart_export.txt'

#Read in biomart file that will match ensembl IDs to gene names
biomart = fread(biomartPath)
biomart <- biomart %>%
  dplyr::select(`Gene stable ID version`, `Gene name`)

#Create function for updating genes with the biomart file
updateGenes <- function(x) {
  #Merge the results and the biomart table
  resultsDF = as.data.frame(x, row.names = rownames(x))
  resultsDF$gene = rownames(resultsDF)
  resultsDF <- unique(left_join(x = resultsDF, y = biomart, join_by(x$gene == y$`Gene stable ID version`)), )
  rownames(resultsDF) <- rownames(x)
  
  #Replace ensemble IDs with the gene names and remove the merged column to remove redundancy
  resultsDF <- resultsDF %>%
    mutate(gene = ifelse(`Gene name` %in% "" | `Gene name` %in% NA, gene, `Gene name`)) %>%
    dplyr::select(!`Gene name`) %>%
  return(resultsDF)
}

#This is for creating the salmon tximport file and saving it#
# -----------------------------------------------------------#
# txdb = makeTxDbFromGFF(gtfFile)
# saveDb(txdb, txdbFile)
# txdb <- loadDb(txdbFile)
# key = keys(txdb, keytype = "TXNAME")
# tx2gene = AnnotationDbi::select(txdb, key, "GENEID", "TXNAME")
# 
# salmonFiles = list.files(path = salmonFilePath, pattern = "quant.sf", full.names = TRUE, recursive = TRUE)
# names(salmonFiles) <- list.files(salmonFilePath)
# salmonData.g <- tximport(salmonFiles, type = "salmon", tx2gene = tx2gene, txIn = TRUE, txOut = FALSE, countsFromAbundance = "no")
# saveRDS(salmonData.g, rdsFile)
#------------------------------------------------------------#

#Load salmon tximport object
salmonData.g <- readRDS(rdsFile)
exclude <- grep('^lnc|^LINC|^rRNA|^RNA5SP|^Y_RNA|^RPL|^RPS', updateGenes(salmonData.g$counts)$gene, ignore.case = TRUE)
salmonData.g$counts <- salmonData.g$counts[-exclude,]
salmonData.g$abundance <- salmonData.g$abundance[-exclude,]
salmonData.g$length <- salmonData.g$length[-exclude,]

#Read list of pseudogenes and filter from salmon
pseudoGenes <- as.character(fread('/scratch/eleiter-weintraub/rnaSeqProj/dataFiles/refGenomes/ensemblPseudoGenes.txt', header = FALSE)$V1)
salmonData.g$counts <- salmonData.g$counts[!(row.names(salmonData.g$counts) %in% pseudoGenes),]
salmonData.g$abundance <- salmonData.g$abundance[!(row.names(salmonData.g$abundance) %in% pseudoGenes),]
salmonData.g$length <- salmonData.g$length[!(row.names(salmonData.g$length) %in% pseudoGenes),]
  
#Read in star gene count files and format to proper matrix format
names(starFiles) <- list.files(starFilePath)
starCounts = data.frame(fread(starFiles[1]))[c(1,4)]
for(i in 2:length(starFiles))
{
  starCounts = cbind(starCounts, data.frame(fread(starFiles[i]))[4])
}
colnames(starCounts) <- c('Genes', names(starFiles))
starCounts = starCounts[5:nrow(starCounts),]
rownames(starCounts) <- starCounts[,1]
starCounts <- starCounts[,2:ncol(starCounts)]

#Filter star counts for lnc, ribosomal rna, and pseudogenes
starCounts <- updateGenes(starCounts)
exclude <- grep('^lnc|^LINC|rRNA|^RNA5SP|^Y_RNA|^RPL|^RPS', starCounts$gene, ignore.case = TRUE)
starCounts <- starCounts[-exclude,]
starCounts <- starCounts[!(row.names(starCounts) %in% pseudoGenes),]
starCountMatrix <- as.matrix(starCounts[,1:86])


#Read in and format metadata
metadata = fread(metaFile)
str(metadata)
metadata <- metadata %>%
  dplyr::mutate(id = paste0(`DCL Patient ID`, '_', `DCL Sample ID`)) %>%
  dplyr::rename(patient =`DCL Patient ID`, sample = `DCL Sample ID`)
rownames(metadata) <- metadata$id
metadata <- metadata %>%
  dplyr::select(id, patient, sample, Group, Gene, Mutation, Genotype, `Subject Details`) %>%
  dplyr::filter(rownames(metadata) %in% names(salmonFiles))
rownames(metadata)
metadata$Group <- factor(metadata$Group, levels = c('Control', 'MTFMT_Unaffected_HET', 'MT_ATPASE_Control', 'MTFMT_Affcted', 'MT_ATPASE_Affcted', 'ATP6V1A_Affected', 'SURF1_Affected', 'DNM1_Affected'))

#Create DESeq object for Salmon, design has group as factor with control as intercept
#And remove rowmean counts less than 20
salmon.dds <- DESeqDataSetFromTximport(salmonData.g, colData =  metadata, design =  ~ Group)
salmon.dds <- salmon.dds[which(rowMeans(counts(salmon.dds)) > 10),]

#Plot variance of data to help determine good cutoff value and filter at that cutoff
salmon.vars = as.data.frame(rowVars(counts(salmon.dds)))
colnames(salmon.vars) <- c('var')
str(salmon.vars)
salmon.vars <- salmon.vars[which(salmon.vars$var < 20000),]
salmon.vars <- as.data.frame(salmon.vars)
colnames(salmon.vars) <- c('var')
table(salmon.vars$var < 500)
ggplot(data = salmon.vars, aes(x = var)) +
  geom_histogram(bins = 100)

salmon.dds <- salmon.dds[which(rowVars(counts(salmon.dds)) > 500),]


#Finalize creation of salmon deseq object
salmon.dds <- DESeq(salmon.dds)

#Create DESeq object for STAR, design has group as factor with control as intercept
#And remove rowmean counts less than 20
star.dds = DESeqDataSetFromMatrix(starCountMatrix, colData = metadata, design = ~ Group)
star.dds <- star.dds[which(rowMeans(counts(star.dds)) > 10),]

#Plot variance of data to help determine good cutoff value and filter at that cutoff
star.vars = as.data.frame(rowVars(counts(star.dds)))
colnames(star.vars) <- c('var')
str(star.vars)
star.vars <- star.vars[which(star.vars$var < 10000),]
star.vars <- as.data.frame(star.vars)
colnames(star.vars) <- c('var')
table(star.vars$var < 500)
ggplot(data = star.vars, aes(x = var)) +
  geom_histogram(bins = 100)

star.dds <- star.dds[which(rowVars(counts(star.dds)) > 500),]

star.dds <- DESeq(star.dds)

#Create results from salmon
resultsNames(salmon.dds)
salmon.dds.control.mtfmt = updateGenes(as.data.frame(results(salmon.dds, name = 'Group_MTFMT_Affcted_vs_Control', alpha = 0.05)))
salmon.dds.control.ATP6 = updateGenes(as.data.frame(results(salmon.dds, name = 'Group_MT_ATPASE_Affcted_vs_Control', alpha = 0.05)))
salmon.dds.control.dnm1 = updateGenes(as.data.frame(results(salmon.dds, name = 'Group_DNM1_Affected_vs_Control', alpha = 0.05)))
salmon.dds.control.surf1 = updateGenes(as.data.frame(results(salmon.dds, name = 'Group_SURF1_Affected_vs_Control', alpha = 0.05)))
salmon.dds.control.atp6 = updateGenes(as.data.frame(results(salmon.dds, name = 'Group_ATP6V1A_Affected_vs_Control', alpha = 0.05)))
salmon.dds.control.unaff = updateGenes(as.data.frame(results(salmon.dds, name = 'Group_MTFMT_Unaffected_HET_vs_Control', alpha = 0.05)))
salmon.dds.mtfmt = updateGenes(as.data.frame(results(salmon.dds, contrast = c("Group", "MTFMT_Affcted", "MTFMT_Unaffected_HET"), alpha = 0.05)))
salmon.dds.ATP6 = updateGenes(as.data.frame(results(salmon.dds, contrast = c("Group", "MT_ATPASE_Affcted", "MT_ATPASE_Control"), alpha = 0.05)))

#Create results from star
resultsNames(star.dds)
star.dds.control.mtfmt = updateGenes(as.data.frame(results(star.dds, name = 'Group_MTFMT_Affcted_vs_Control', alpha = 0.05)))
star.dds.control.ATP6 = updateGenes(as.data.frame(results(star.dds, name = 'Group_MT_ATPASE_Affcted_vs_Control', alpha = 0.05)))
star.dds.control.dnm1 = updateGenes(as.data.frame(results(star.dds, name = 'Group_DNM1_Affected_vs_Control', alpha = 0.05)))
star.dds.control.surf1 = updateGenes(as.data.frame(results(star.dds, name = 'Group_SURF1_Affected_vs_Control', alpha = 0.05)))
star.dds.control.atp6 = updateGenes(as.data.frame(results(star.dds, name = 'Group_ATP6V1A_Affected_vs_Control', alpha = 0.05)))
star.dds.control.unaff = updateGenes(as.data.frame(results(star.dds, name = 'Group_MTFMT_Unaffected_HET_vs_Control', alpha = 0.05)))
star.dds.mtfmt = updateGenes(as.data.frame(results(star.dds, contrast = c("Group", "MTFMT_Affcted", "MTFMT_Unaffected_HET"), alpha = 0.05)))
star.dds.ATP6 = updateGenes(as.data.frame(results(star.dds, contrast = c("Group", "MT_ATPASE_Affcted", "MT_ATPASE_Control"), alpha = 0.05)))

#Salmon volcano plots
EnhancedVolcano(salmon.dds.control.mtfmt,
                lab = salmon.dds.control.mtfmt$gene,
                title = 'MTFMT Homozygous Affected vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/Salmon/salmon.dds.control.mtfmt.png", width = 9.2, height = 11.9, create.dir = TRUE)

EnhancedVolcano(salmon.dds.control.ATP6,
                lab = salmon.dds.control.ATP6$gene,
                title = 'MT-ATP6 Affected vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                pCutoff = 1e-15, pointSize = 3) +
  scale_y_continuous(limits = c(0, 100))+
  theme(
    legend.position = "none",
    axis.text.x = element_text(family = "Arial", size = 18),
    axis.text.y = element_text(size = 30, family = "Arial"),
    axis.title.x = element_text(size = 35),
    axis.title.y = element_text(size = 35),
    plot.title = element_text(hjust = 0.5, size = 30, family = "Arial", face = "plain")
  )
ggsave("Plots/Salmon/salmon.dds.control.ATP6.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.control.dnm1,
                lab = salmon.dds.control.dnm1$gene,
                title = 'DNM1 Affected vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                pCutoff = 1e-7)
ggsave("Plots/Salmon/salmon.dds.control.dnm1.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.control.surf1,
                lab = salmon.dds.control.surf1$gene,
                title = 'SURF1 Affected vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                pCutoff = 1e-6)
ggsave("Plots/Salmon/salmon.dds.control.surf1.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.control.atp6,
                lab = salmon.dds.control.atp6$gene,
                title = 'ATP6V1A Affected vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                pCutoff = 1e-8)
ggsave("Plots/Salmon/salmon.dds.control.ATP6V1A.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.control.unaff,
                lab = salmon.dds.control.unaff$gene,
                title = 'MTFMT Unaffected Heterozygous vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/Salmon/salmon.dds.control.unaff.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.mtfmt,
                lab = salmon.dds.mtfmt$gene,
                title = 'MTFMT Unaffected Heterozygous vs Affected Homozygous',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                pCutoff = 1e-2)
ggsave("Plots/Salmon/salmon.dds.mtfmt.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.ATP6,
                lab = salmon.dds.ATP6$gene,
                title = 'MT-ATP6 Affected vs Unaffected',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                pCutoff=1e-9)
ggsave("Plots/Salmon/salmon.dds.ATP6.png", width = 9.2, height = 11.9)

#STAR volcano plots
EnhancedVolcano(star.dds.control.mtfmt,
                lab = star.dds.control.mtfmt$gene,
                title = 'MTFMT Homozygous Affected vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.mtfmt.png", width = 9.2, height = 11.9, create.dir = TRUE)

EnhancedVolcano(star.dds.control.ATP6,
                lab = star.dds.control.ATP6$gene,
                title = 'ATP6 Homozygous Affected vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.ATP6.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.control.dnm1,
                lab = star.dds.control.dnm1$gene,
                title = 'DNM1 Affected vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.dnm1.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.control.surf1,
                lab = star.dds.control.surf1$gene,
                title = 'SURF1 Affected vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.surf1.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.control.atp6,
                lab = star.dds.control.atp6$gene,
                title = 'ATP6V1A Homozygous Affected vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.atp6.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.control.unaff,
                star.dds.control.unaff$gene,
                title = 'MTFMT Unaffected Heterozygous vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.unaff.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.mtfmt,
                lab = star.dds.mtfmt$gene,
                title = 'MTFMT Unaffected Heterozygous vs Affected Homozygous',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.mtfmt.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.ATP6,
                lab = star.dds.ATP6$gene,
                title = 'MT-ATP6 Affected vs Controls',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.ATP6.png", width = 9.2, height = 11.9)


#Lets make tables of the significantly differentially expressed genes
getSignificant <- function(x, p, f) {
  significant.df <- x[which(x$padj <= p), ]
  significant.df <- significant.df[which(abs(significant.df$log2FoldChange) >= f),]
  return(significant.df)
}

#Pull out significantly differentially expressed genes from both
salmon.dds.control.mtfmt.sig1 <- getSignificant(salmon.dds.control.mtfmt, 0.05, 2)
salmon.dds.control.atp6.sig1 <- getSignificant(salmon.dds.control.atp6, 0.05, 2)
salmon.dds.control.dnm1.sig1 <- getSignificant(salmon.dds.control.dnm1, 0.05, 2)
salmon.dds.control.ATP6.sig1 <- getSignificant(salmon.dds.control.ATP6, 0.05, 2)
salmon.dds.control.surf1.sig1 <- getSignificant(salmon.dds.control.surf1, 0.05, 2)
salmon.dds.control.unaff.sig1 <- getSignificant(salmon.dds.control.unaff, 0.05, 2)
salmon.dds.ATP6.sig1 <- getSignificant(salmon.dds.ATP6, 0.05, 2)
salmon.dds.mtfmt.sig1 <- getSignificant(salmon.dds.mtfmt, 0.05, 2)

star.dds.control.mtfmt.sig1 <- getSignificant(star.dds.control.mtfmt, 0.05, 2)
star.dds.control.atp6.sig1 <- getSignificant(star.dds.control.atp6, 0.05, 2)
star.dds.control.ATP6.sig1 <- getSignificant(star.dds.control.ATP6, 0.05, 2)
star.dds.control.dnm1.sig1 <- getSignificant(star.dds.control.dnm1, 0.05, 2)
star.dds.control.surf1.sig1 <- getSignificant(star.dds.control.surf1, 0.05, 2)
star.dds.control.unaff.sig1 <- getSignificant(star.dds.control.unaff, 0.05, 2)
star.dds.ATP6.sig1 <- getSignificant(star.dds.ATP6, 0.05, 2)
star.dds.mtfmt.sig1 <- getSignificant(star.dds.mtfmt, 0.05, 2)

salmon.dds.control.mtfmt.sig2 <- getSignificant(salmon.dds.control.mtfmt, 0.01, 1.2)
salmon.dds.control.atp6.sig2 <- getSignificant(salmon.dds.control.atp6, 0.01, 1.2)
salmon.dds.control.ATP6.sig2 <- getSignificant(salmon.dds.control.ATP6, 0.01, 1.2)
salmon.dds.control.dnm1.sig2 <- getSignificant(salmon.dds.control.dnm1, 0.01, 1.2)
salmon.dds.control.surf1.sig2 <- getSignificant(salmon.dds.control.surf1, 0.01, 1.2)
salmon.dds.control.unaff.sig2 <- getSignificant(salmon.dds.control.unaff, 0.01, 1.2)
salmon.dds.ATP6.sig2 <- getSignificant(salmon.dds.ATP6, 0.01, 1.2)
salmon.dds.mtfmt.sig2 <- getSignificant(salmon.dds.mtfmt, 0.01, 1.2)

star.dds.control.mtfmt.sig2 <- getSignificant(star.dds.control.mtfmt, 0.01, 1.2)
star.dds.control.atp6.sig2 <- getSignificant(star.dds.control.atp6, 0.01, 1.2)
star.dds.control.ATP6.sig2 <- getSignificant(star.dds.control.ATP6, 0.01, 1.2)
star.dds.control.dnm1.sig2 <- getSignificant(star.dds.control.dnm1, 0.01, 1.2)
star.dds.control.surf1.sig2 <- getSignificant(star.dds.control.surf1, 0.01, 1.2)
star.dds.control.unaff.sig2 <- getSignificant(star.dds.control.unaff, 0.01, 1.2)
star.dds.ATP6.sig2 <- getSignificant(star.dds.ATP6, 0.01, 1.2)
star.dds.mtfmt.sig2 <- getSignificant(star.dds.mtfmt, 0.01, 1.2)


library(UpSetR)
GeneList1 <- list(salmon_control_mtfmt = c(salmon.dds.control.mtfmt.sig1$gene), 
                 salmon_control_atp6 = c(salmon.dds.control.atp6.sig1$gene), 
                 salmon_control_ATP6 = c(salmon.dds.control.ATP6.sig1$gene), 
                 salmon_control_dnm1 = c(salmon.dds.control.dnm1.sig1$gene), 
                 salmon_control_surf1 = c(salmon.dds.control.surf1.sig1$gene),
                 salmon_control_mtfmtunaff = c(salmon.dds.control.unaff.sig1$gene),
                 salmon_ATP6 = c(salmon.dds.ATP6.sig1$gene), 
                 salmon_mtfmt = c(salmon.dds.mtfmt.sig1$gene),  
                 star_control_mtfmt = c(star.dds.control.mtfmt.sig1$gene), 
                 star_control_atp6 = c(star.dds.control.atp6.sig1$gene), 
                 star_control_ATP6 = c(star.dds.control.ATP6.sig1$gene), 
                 star_control_dnm1 = c(star.dds.control.dnm1.sig1$gene), 
                 star_control_surf1 = c(star.dds.control.surf1.sig1$gene), 
                 star_control_mtfmtunaff = c(star.dds.control.unaff.sig1$gene),
                 star_ATP6 = c(star.dds.ATP6.sig1$gene), 
                 star_mtfmt = c(star.dds.mtfmt.sig1$gene))
upsetGenes <- fromList(GeneList1)
upset(upsetGenes, nsets = 16)

GeneList1 <- list(`Salmon control vs MTFMT Affected` = c(salmon.dds.control.mtfmt.sig2$gene), 
                  `Salmon control vs ATP6V1A` = c(salmon.dds.control.atp6.sig2$gene), 
                  `Salmon control vs MT-ATP6` = c(salmon.dds.control.ATP6.sig2$gene), 
                  `Salmon control vs DNM1` = c(salmon.dds.control.dnm1.sig2$gene), 
                  `Salmon control vs SURF1` = c(salmon.dds.control.surf1.sig2$gene),
                  `Salmon control vs MTFMT Unaffected` = c(salmon.dds.control.unaff.sig2$gene),
                  `Salmon MT-ATP6 Affected vs Unaffected` = c(salmon.dds.ATP6.sig2$gene), 
                  `Salmon MTFMT Affected vs Unaffected` = c(salmon.dds.mtfmt.sig2$gene),  
                  `STAR control vs MTFMT Affected` = c(star.dds.control.mtfmt.sig2$gene), 
                  `STAR control vs ATP6V1A` = c(star.dds.control.atp6.sig2$gene), 
                  `STAR control vs MT-ATP6` = c(star.dds.control.ATP6.sig2$gene), 
                  `STAR control vs DNM1` = c(star.dds.control.dnm1.sig2$gene), 
                  `STAR control vs SURF1` = c(star.dds.control.surf1.sig2$gene),
                  `STAR control vs MTFMT Unaffected` = c(star.dds.control.unaff.sig2$gene),
                  `STAR MT-ATP6 Affected vs Unaffected` = c(star.dds.ATP6.sig2$gene), 
                  `STAR MTFMT Affected vs Unaffected` = c(star.dds.mtfmt.sig2$gene))
upsetGenes <- fromList(GeneList1)
upset(upsetGenes, nsets = 16, text.scale = 7, point.size = 10)
ggsave("Plots/Salmon/UpsetPlot.png", width = 9.2, height = 11.9)

#Load in mitocarta mitochondrial genes
mitocartaPath <- '/scratch/eleiter-weintraub/rnaSeqProj/dataFiles/Human.MitoCarta3.0.xls'
mitocarta <- readxl::read_xls(mitocartaPath, sheet = 'A Human MitoCarta3.0')
mitocartaGeneNames <- mitocarta$Symbol
mitocartaGeneNames <- c(mitocartaGeneNames, unlist(strsplit(mitocarta$Synonyms, '\\|')))

getMitoGenes <- function(x) {
  row_indices <- which(x$gene %in% mitocartaGeneNames)
  x.mito <- x[row_indices, ]
  return(x.mito)
}

#Pull out mitochondrial-associated genes from both
salmon.dds.control.mtfmt.mito <- getMitoGenes(updateGenes(salmon.dds.control.mtfmt))
salmon.dds.control.atp6.mito <- getMitoGenes(updateGenes(salmon.dds.control.atp6))
salmon.dds.control.ATP6.mito <- getMitoGenes(updateGenes(salmon.dds.control.ATP6))
salmon.dds.control.dnm1.mito <- getMitoGenes(updateGenes(salmon.dds.control.dnm1))
salmon.dds.control.surf1.mito <- getMitoGenes(updateGenes(salmon.dds.control.surf1))
salmon.dds.control.unaff.mito <- getMitoGenes(updateGenes(salmon.dds.control.unaff))
salmon.dds.ATP6.mito <- getMitoGenes(updateGenes(salmon.dds.ATP6))
salmon.dds.mtfmt.mito <- getMitoGenes(updateGenes(salmon.dds.mtfmt))

star.dds.control.mtfmt.mito <- getMitoGenes(star.dds.control.mtfmt)
star.dds.control.atp6.mito <- getMitoGenes(star.dds.control.atp6)
star.dds.control.ATP6.mito <- getMitoGenes(star.dds.control.ATP6)
star.dds.control.dnm1.mito <- getMitoGenes(star.dds.control.dnm1)
star.dds.control.surf1.mito <- getMitoGenes(star.dds.control.surf1)
star.dds.control.unaff.mito <- getMitoGenes(star.dds.control.unaff)
star.dds.ATP6.mito <- getMitoGenes(star.dds.ATP6)
star.dds.mtfmt.mito <- getMitoGenes(star.dds.mtfmt)

#Plot volcano plots of mitochodrial genes for both
EnhancedVolcano(salmon.dds.control.mtfmt.mito,
                lab = salmon.dds.control.mtfmt.mito$gene,
                title = 'MTFMT Homozygous Affected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/Salmon/salmon.dds.control.mtfmt.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.control.ATP6.mito,
                lab = salmon.dds.control.ATP6.mito$gene,
                title = 'MT-ATP6 Affected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2, pointSize = 3) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(family = "Arial", size = 18),
    axis.text.y = element_text(size = 30, family = "Arial"),
    axis.title.x = element_text(size = 35),
    axis.title.y = element_text(size = 35),
    plot.title = element_text(hjust = 0.5, size = 25, family = "Arial", face = "plain")
  )
ggsave("Plots/Salmon/salmon.dds.control.ATP6.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.control.dnm1.mito,
                lab = salmon.dds.control.dnm1.mito$gene,
                title = 'DNM1 Affected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/Salmon/salmon.dds.control.dnm1.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.control.surf1.mito,
                lab = salmon.dds.control.surf1.mito$gene,
                title = 'SURF1 Affected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/Salmon/salmon.dds.control.surf1.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.control.atp6.mito,
                lab = salmon.dds.control.atp6.mito$gene,
                title = 'ATP6V1A Affected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/Salmon/salmon.dds.control.atp6.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.control.unaff.mito,
                lab = salmon.dds.control.unaff.mito$gene,
                title = 'MTFMT Unaffected Heterozygous vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/Salmon/salmon.dds.control.unaff.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.mtfmt.mito,
                lab = salmon.dds.mtfmt.mito$gene,
                title = 'MTFMT Unaffected Heterozygous vs Affected Homozygous Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/Salmon/salmon.dds.mtfmt.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(salmon.dds.ATP6.mito,
                lab = salmon.dds.ATP6.mito$gene,
                title = 'MT-ATP6 Affected vs Unaffected Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/Salmon/salmon.dds.ATP6.mito.png", width = 9.2, height = 11.9)


#STAR volcano plots
EnhancedVolcano(star.dds.control.mtfmt.mito,
                lab = star.dds.control.mtfmt.mito$gene,
                title = 'MTFMT Homozygous Affected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.mtfmt.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.control.ATP6.mito,
                lab = star.dds.control.ATP6.mito$gene,
                title = 'ATP6 Homozygous Affected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.ATP6.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.control.dnm1.mito,
                lab = star.dds.control.dnm1.mito$gene,
                title = 'DNM1 Affected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.dnm1.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.control.surf1.mito,
                lab = star.dds.control.surf1.mito$gene,
                title = 'SURF1 Affected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.surf1.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.control.atp6.mito,
                lab = star.dds.control.atp6.mito$gene,
                title = 'ATP6V1A Homozygous Affected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.atp6.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.control.unaff.mito,
                lab = star.dds.control.unaff.mito$gene,
                title = 'MTFMT Heterozygous Unaffected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.control.unaff.mito.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.mtfmt.mito,
                lab = star.dds.mtfmt.mito$gene,
                title = 'MTFMT Unaffected Heterozygous vs Affected Homozygous Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.mtfmt.mito.png.png", width = 9.2, height = 11.9)

EnhancedVolcano(star.dds.ATP6.mito,
                lab = star.dds.ATP6.mito$gene,
                title = 'MT-ATP6 Affected vs Controls Mitochondrial Genes',
                x = 'log2FoldChange',
                y = 'padj',
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 0.2)
ggsave("Plots/STAR/star.dds.ATP6.mito.png", width = 9.2, height = 11.9)

#Function to pull out genes for pathway analysis
forPathway <- function(x) {
  pathway <- x[which(x$padj < 0.25),]
  pathway <- pathway[order(pathway$padj),]
}

#Pull out the genes for pathway analysis and take only top 2000 if there are more than 2000
star.dds.control.mtfmt.pathway <- forPathway(star.dds.control.mtfmt)
star.dds.mtfmt.pathway <- forPathway(star.dds.mtfmt)
star.dds.control.dnm1.pathway <- forPathway(star.dds.control.dnm1)
star.dds.control.dnm1.pathway <- star.dds.control.dnm1.pathway[1:2000,]
star.dds.control.atp6.pathway <- forPathway(star.dds.control.atp6)
star.dds.control.atp6.pathway <- star.dds.control.atp6.pathway[1:2000,]
star.dds.control.ATP6.pathway <- forPathway(star.dds.control.ATP6)
star.dds.control.ATP6.pathway <- star.dds.control.ATP6.pathway [1:2000,]

salmon.dds.control.mtfmt.pathway <- forPathway(salmon.dds.control.mtfmt)
salmon.dds.mtfmt.pathway <- forPathway(salmon.dds.mtfmt)
salmon.dds.control.dnm1.pathway <- forPathway(salmon.dds.control.dnm1)
salmon.dds.control.dnm1.pathway <- salmon.dds.control.dnm1.pathway[1:2000,]
salmon.dds.control.atp6.pathway <- forPathway(salmon.dds.control.atp6)
salmon.dds.control.atp6.pathway <- salmon.dds.control.atp6.pathway[1:2000,]
salmon.dds.control.ATP6.pathway <- forPathway(salmon.dds.control.ATP6)
salmon.dds.control.ATP6.pathway <- salmon.dds.control.ATP6.pathway[1:2000,]
salmon.dds.control.surf1.pathway <- forPathway(salmon.dds.control.surf1)[1:2000,]

write.csv(star.dds.control.mtfmt.pathway$gene, 'pathwayList/star/control.mtfmt.csv', quote = FALSE)
write.csv(star.dds.mtfmt.pathway$gene, 'pathwayList/star/mtfmt.csv', quote = FALSE)
write.csv(star.dds.control.dnm1.pathway$gene, 'pathwayList/star/control.dnm1.csv', quote = FALSE)
write.csv(star.dds.control.atp6.pathway$gene, 'pathwayList/star/control.atp6.csv', quote = FALSE)
write.csv(star.dds.control.ATP6.pathway$gene, 'pathwayList/star/control.ATP6.csv', quote = FALSE)

write.csv(salmon.dds.control.mtfmt.pathway$gene, 'pathwayList/salmon/control.mtfmt.csv', quote = FALSE)
write.csv(salmon.dds.mtfmt.pathway$gene, 'pathwayList/salmon/mtfmt.csv', quote = FALSE)
write.csv(salmon.dds.control.dnm1.pathway$gene, 'pathwayList/salmon/control.dnm1.csv', quote = FALSE)
write.csv(salmon.dds.control.atp6.pathway$gene, 'pathwayList/salmon/control.atp6.csv', quote = FALSE)
write.csv(salmon.dds.control.ATP6.pathway$gene, 'pathwayList/salmon/control.ATP6.csv', quote = FALSE)
write.csv(salmon.dds.control.surf1.pathway$gene, 'pathwayList/salmon/control.surf1.csv', quote = FALSE)

#Import genes that are part of oxidative phosphorylation
oxphos <- names(as.list(fread("OxPhos.csv")))

#Pull out oxphos genes for each salmon analysis group
oxphosDiffExp <- data.frame(sort(oxphos))
colnames(oxphosDiffExp) <- c('gene')
row.names(oxphosDiffExp) <- oxphosDiffExp$gene
oxphosDiffExp <- left_join(x = oxphosDiffExp, y = salmon.dds.ATP6[,c('log2FoldChange', 'gene')], join_by(x$gene == y$gene))
oxphosDiffExp <- left_join(x = oxphosDiffExp, y = salmon.dds.control.ATP6[,c('log2FoldChange', 'gene')], join_by(x$gene == y$gene))
oxphosDiffExp <- left_join(x = oxphosDiffExp, y = salmon.dds.mtfmt[,c('log2FoldChange', 'gene')], join_by(x$gene == y$gene))
oxphosDiffExp <- left_join(x = oxphosDiffExp, y = salmon.dds.control.mtfmt[,c('log2FoldChange', 'gene')], join_by(x$gene == y$gene))
oxphosDiffExp <- left_join(x = oxphosDiffExp, y = salmon.dds.control.unaff[,c('log2FoldChange', 'gene')], join_by(x$gene == y$gene))
oxphosDiffExp <- left_join(x = oxphosDiffExp, y = salmon.dds.control.atp6[,c('log2FoldChange', 'gene')], join_by(x$gene == y$gene))
oxphosDiffExp <- left_join(x = oxphosDiffExp, y = salmon.dds.control.dnm1[,c('log2FoldChange', 'gene')], join_by(x$gene == y$gene))
oxphosDiffExp <- left_join(x = oxphosDiffExp, y = salmon.dds.control.surf1[,c('log2FoldChange', 'gene')], join_by(x$gene == y$gene))
colnames(oxphosDiffExp) <- c('gene', 'MT-ATP6 Unaffected vs Affected', 'Controls vs MT-ATP6 Affected', 'MTFMT Unaffected vs Affected',
                             'Controls vs MTFMT Affected', 'Controls vs MTFMT Unaffected', 'Controls vs ATP6V1A', 'Controls vs DNM1', 'Controls vs SURF1')
oxphosDiffExp <- drop_na(oxphosDiffExp)

#Reformat all of them for heatmap
library(reshape)
tileDF <- reshape::melt(oxphosDiffExp)

#Create heatmap of oxidative phosphorylation gene fold changes
ggplot(tileDF, aes(x = variable, y = gene, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = "#0000FF",
                                     mid = "#FFFFFF",
                                     high = "#FF0000") +
  theme( 
    axis.text.x = element_text(angle = 25, hjust = 1)
  )
ggsave('Plots/Salmon/oxphosHeatMap.png', units = "in", width = 15, height = 30)


#Create DF for number of mapped reads per patient
salmonCountPlot <- as.data.frame(t(counts(salmon.dds)))
salmonCountPlot$identifier <- rownames(salmonCountPlot)  
salmonCountPlot <- salmonCountPlot %>%
  separate(col = "identifier", sep = '_', into = c("ProjectID", "Patient", "LibraryID"), remove = FALSE) %>%
  unite(col = "PatientID", sep = '_', ProjectID, Patient, remove = TRUE)
salmonCountPlot$readCounts <- rowSums(salmonCountPlot[,1:14747])/1000000

#Boxplots of number of mapped reads per patient
ggplot(data = salmonCountPlot) + 
  geom_boxplot(aes(x = PatientID, y = readCounts, fill = PatientID)) + 
  ggtitle('Number of aligned read counts for each sample by patient') +
  ylab("Mapped read counts (10^6)") +
  xlab("Patient ID") + 
  geom_point(aes(x = PatientID, y = readCounts), position = position_dodge2(width = 0.2)) + 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family = "Arial", size = 18),
    axis.text.y = element_text(size = 30, family = "Arial"),
    axis.title.x = element_text(size = 35),
    axis.title.y = element_text(size = 35),
    plot.title = element_text(hjust = 0.5, size = 30, family = "Arial")
  )
ggsave("Plots/Salmon/ReadCounts.png", width = 13, height = 13)

salmon.dds.control.mtfmt.mito.sig <- getSignificant(salmon.dds.control.mtfmt.mito, 0.01, 1.2)
salmon.dds.control.atp6.mito.sig <- getSignificant(salmon.dds.control.atp6.mito, 0.01, 1.2)
salmon.dds.control.ATP6.mito.sig <- getSignificant(salmon.dds.control.ATP6.mito, 0.01, 1.2)
salmon.dds.control.dnm1.mito.sig <- getSignificant(salmon.dds.control.dnm1.mito, 0.01, 1.2)
salmon.dds.control.surf1.mito.sig <- getSignificant(salmon.dds.control.surf1.mito, 0.01, 1.2)
salmon.dds.control.unaff.mito.sig <- getSignificant(salmon.dds.control.unaff.mito, 0.01, 1.2)
salmon.dds.ATP6.mito.sig <- getSignificant(salmon.dds.ATP6.mito, 0.01, 1.2)
salmon.dds.mtfmt.mito.sig <- getSignificant(salmon.dds.mtfmt.mito, 0.01, 1.2)

star.dds.control.mtfmt.mito.sig <- getSignificant(star.dds.control.mtfmt.mito, 0.01, 1.2)
star.dds.control.atp6.mito.sig <- getSignificant(star.dds.control.atp6.mito, 0.01, 1.2)
star.dds.control.ATP6.mito.sig <- getSignificant(star.dds.control.ATP6.mito, 0.01, 1.2)
star.dds.control.dnm1.mito.sig <- getSignificant(star.dds.control.dnm1.mito, 0.01, 1.2)
star.dds.control.surf1.mito.sig <- getSignificant(star.dds.control.surf1.mito, 0.01, 1.2)
star.dds.control.unaff.mito.sig <- getSignificant(star.dds.control.unaff.mito, 0.01, 1.2)
star.dds.ATP6.mito.sig <- getSignificant(star.dds.ATP6.mito, 0.01, 1.2)
star.dds.mtfmt.mito.sig <- getSignificant(star.dds.mtfmt.mito, 0.01, 1.2)

GeneList1 <- list(`Salmon control vs MTFMT Affected` = c(salmon.dds.control.mtfmt.mito.sig$gene), 
                  `Salmon control vs ATP6V1A` = c(salmon.dds.control.atp6.mito.sig$gene), 
                  `Salmon control vs MT-ATP6` = c(salmon.dds.control.ATP6.mito.sig$gene), 
                  `Salmon control vs DNM1` = c(salmon.dds.control.dnm1.mito.sig$gene), 
                  `Salmon control vs SURF1` = c(salmon.dds.control.surf1.mito.sig$gene),
                  `Salmon control vs MTFMT Unaffected` = c(salmon.dds.control.unaff.mito.sig$gene),
                  `Salmon MT-ATP6 Affected vs Unaffected` = c(salmon.dds.ATP6.mito.sig$gene), 
                  `Salmon MTFMT Affected vs Unaffected` = c(salmon.dds.mtfmt.mito.sig$gene),  
                  `STAR control vs MTFMT Affected` = c(star.dds.control.mtfmt.mito.sig$gene), 
                  `STAR control vs ATP6V1A` = c(star.dds.control.atp6.mito.sig$gene), 
                  `STAR control vs MT-ATP6` = c(star.dds.control.ATP6.mito.sig$gene), 
                  `STAR control vs DNM1` = c(star.dds.control.dnm1.mito.sig$gene), 
                  `STAR control vs SURF1` = c(star.dds.control.surf1.mito.sig$gene),
                  `STAR control vs MTFMT Unaffected` = c(star.dds.control.unaff.mito.sig$gene),
                  `STAR MT-ATP6 Affected vs Unaffected` = c(star.dds.ATP6.mito.sig$gene), 
                  `STAR MTFMT Affected vs Unaffected` = c(star.dds.mtfmt.mito.sig$gene))
upsetGenes <- fromList(GeneList1)
upset(upsetGenes, nsets = 16, text.scale = 7, point.size = 10)

control.mtfmt.mito.overlap <- salmon.dds.control.mtfmt[intersect(rownames(salmon.dds.control.mtfmt.mito.sig), rownames(star.dds.control.mtfmt.mito.sig)), ]
control.atp6.mito.overlap <- salmon.dds.control.atp6[intersect(rownames(salmon.dds.control.atp6.mito.sig), rownames(star.dds.control.atp6.mito.sig)), ]
control.ATP6.mito.overlap <- salmon.dds.control.ATP6[intersect(rownames(salmon.dds.control.ATP6.mito.sig), rownames(star.dds.control.ATP6.mito.sig)), ]
control.dnm1.mito.overlap <- salmon.dds.control.dnm1[intersect(rownames(salmon.dds.control.dnm1.mito.sig), rownames(star.dds.control.dnm1.mito.sig)), ]
control.surf1.mito.overlap <- salmon.dds.control.surf1[intersect(rownames(salmon.dds.control.surf1.mito.sig), rownames(star.dds.control.surf1.mito.sig)), ]
control.unaff.mito.overlap <- salmon.dds.control.unaff[intersect(rownames(salmon.dds.control.unaff.mito.sig), rownames(star.dds.control.unaff.mito.sig)), ]
ATP6.mito.overlap <- salmon.dds.ATP6[intersect(rownames(salmon.dds.ATP6.mito.sig), rownames(star.dds.ATP6.mito.sig)), ]
mtfmt.mito.overlap <- salmon.dds.mtfmt[intersect(rownames(salmon.dds.mtfmt.mito.sig), rownames(star.dds.mtfmt.mito.sig)), ]