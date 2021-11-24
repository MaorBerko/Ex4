# Ex4

---
title: "Ex4 answers"
author: "Maor Berkovich"
date: "11/18/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### Functional enrichment analysis

Installation of the new required packages:
```{r eval=F, echo=T}
install.packages("gProfileR")
install.packages("knitr")
nBiocManager::install("gage")
```

Loading the required packages:
```{r message = F, warning=FALSE}
library(compGenomRData)
library(DESeq2)
library(gProfileR)
library(gage)
```

Loading the data:
```{r message = F, warning=FALSE}
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                           package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv",
                            package = "compGenomRData")
counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))

#remove the 'width' column
countData <- as.matrix(subset(counts, select = c(-width)))

#define the experimental setup
colData <- read.table(coldata_file, header = T, sep = '\t',
                      stringsAsFactors = TRUE)

#define the design formula
designFormula <- "~ group"

#create a DESeq dataset object from the count matrix and the colData
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = as.formula(designFormula))
dds <- DESeq(dds)

DEresults = results(dds, contrast = c("group", 'CASE', 'CTRL'))
```

Let’s select the genes that are significantly differentially expressed between the case and control samples. Let’s extract genes that have an adjusted p-value below 0.1 and that show a 2-fold change (either negative or positive) in the case compared to control. We will then feed this gene set into the gProfileR function.

```{r message = F, warning=FALSE}
#remove genes with NA values
DE <- DEresults[!is.na(DEresults$padj),]

#select genes with adjusted p-values below 0.1
DE <- DE[DE$padj < 0.1,]

#select genes with absolute log2 fold change above 1 (two-fold change)
DE <- DE[abs(DE$log2FoldChange) > 1,]

#get the list of genes of interest
genesOfInterest <- rownames(DE)

#calculate enriched GO terms
goResults <- gprofiler(query = genesOfInterest,
                       organism = 'hsapiens',
                       src_filter = 'GO',
                       hier_filtering = 'moderate')
```

1. Re-run gProfileR, this time using pathway annotations such as KEGG, REACTOME, and protein complex databases such as CORUM, in addition to the GO terms. Sort the resulting tables by columns `precision` and/or `recall`. How do the top GO terms change when sorted for `precision`, `recall`, or `p.value`? hint: use `order()` for sorting. [Difficulty: **Beginner**]
```{r}
# Biological pathways
goResultsKEGG <- gprofiler(query = genesOfInterest,
                       organism = 'hsapiens',
                       src_filter = 'KEGG' ,
                        hier_filtering = 'moderate')
dfgoResultsKEGG<-as.data.frame(goResultsKEGG)
goResultsKEGGOrderByPrecision <-dfgoResultsKEGG[order(dfgoResultsKEGG$precision, decreasing = TRUE),9]
head(goResultsKEGGOrderByPrecision)

goResultsKEGGOrderByRecall <- dfgoResultsKEGG[order(dfgoResultsKEGG$recall, decreasing = TRUE),9]
head(goResultsKEGGOrderByRecall)

goResultsKEGGOrderByPvalue <- dfgoResultsKEGG[order(dfgoResultsKEGG$p.value, decreasing = TRUE),9]
head(goResultsKEGGOrderByPvalue)

# We don't get any results when using the 'REACTOME' filter.
goResultsREACTOME <- gprofiler(query = genesOfInterest,
                       organism = 'hsapiens',
                       src_filter = 'REACTOME' ,
                       hier_filtering = 'moderate')
head(goResultsREACTOME)

# Protein databases
# We don't get any results when using the 'CORUM' filter. 
goResultsCORUM <- gprofiler(query = genesOfInterest,
                       organism = 'hsapiens',
                       src_filter = 'CORUM' ,
                       hier_filtering = 'moderate')
head(goResultsCORUM)
```
#### Gene set enrichment analysis

We use the bioconductor package gage to demonstrate how to do GSEA using normalized expression data of the samples as input.

```{r}
#Let's define the first gene set as the list of genes from one of the
#significant GO terms found in the GO analysis. order go results by pvalue
goResults <- goResults[order(goResults$p.value),]

#restrict the terms that have at most 100 genes overlapping with the query
go <- goResults[goResults$overlap.size < 100,]

# use the top term from this table to create a gene set
geneSet1 <- unlist(strsplit(go[1,]$intersection, ','))

#Define another gene set by just randomly selecting 25 genes from the counts
#table get normalized counts from DESeq2 results
normalizedCounts <- DESeq2::counts(dds, normalized = TRUE)

geneSet2 <- sample(rownames(normalizedCounts), 25)

geneSets <- list('top_GO_term' = geneSet1,
                 'random_set' = geneSet2)

# Using the defined gene sets, we’d like to do a group comparison between the case
# samples with respect to the control samples.

#Use the normalized counts to carry out a GSEA.
gseaResults <- gage(exprs = log2(normalizedCounts+1),
                    ref = match(rownames(colData[colData$group =='CTRL',]),
                                colnames(normalizedCounts)),
                    samp = match(rownames(colData[colData$group == 'CASE',]),
                                 colnames(normalizedCounts)),
                    gsets = geneSets, compare = 'as.group')
```

2. Repeat the gene set enrichment analysis by trying different options for the `compare` argument of the `GAGE:gage`
function. How do the results differ? [Difficulty: **Beginner**]
```{r}
# compare argument is 'unpaired'
gseaResults2 <- gage(exprs = log2(normalizedCounts+1),
                    ref = match(rownames(colData[colData$group =='CTRL',]),
                                colnames(normalizedCounts)),
                    samp = match(rownames(colData[colData$group == 'CASE',]),
                                 colnames(normalizedCounts)),
                    gsets = geneSets, compare = 'unpaired')
gseaResults2

# compare argument is '1ongroup'
gseaResults3 <- gage(exprs = log2(normalizedCounts+1),
                    ref = match(rownames(colData[colData$group =='CTRL',]),
                                colnames(normalizedCounts)),
                    samp = match(rownames(colData[colData$group == 'CASE',]),
                                 colnames(normalizedCounts)),
                    gsets = geneSets, compare = '1ongroup')
gseaResults3
```

3. Make a scatter plot of GO term sizes and obtained p-values by setting the `gProfiler::gprofiler` argument `significant = FALSE`. Is there a correlation of term sizes and p-values? (Hint: Take -log10 of p-values). If so, how can this bias be mitigated? [Difficulty: **Intermediate**]
```{r}
goResultsPlot <-gprofiler(query = genesOfInterest,
                       organism = 'hsapiens',
                      significant = FALSE,
                       src_filter = 'GO',
                       hier_filtering = 'moderate')
goResultsLog10PV <- -log10(goResultsPlot$p.value)

goResultsTermSize <- (goResultsPlot$term.size)
GOScatterPlot <- plot(goResultsTermSize,goResultsLog10PV, xlab='Term size', ylab='-log10 P-value' )

# From the scatter plot we can see there is a positive correlation between the term size and the p-value. As the term size increase so as the p-value. 
```

4. Do a gene-set enrichment analysis using gene sets from top 10 GO terms. [Difficulty: **Intermediate**]
```{r}
# The meaning of 'top 10 GO terms' is the first 10 gene sets with the lowest p-value.
goResults <- goResults[order(goResults$p.value),]
top10 <- goResults[1:10,]

# Generating a separated gene-sets:
geneSet1Top10 <- unlist(strsplit(top10[1,]$intersection, ','))
geneSet2Top10 <- unlist(strsplit(top10[2,]$intersection, ','))
geneSet3Top10 <- unlist(strsplit(top10[3,]$intersection, ','))
geneSet4Top10 <- unlist(strsplit(top10[5,]$intersection, ','))
geneSet5Top10 <- unlist(strsplit(top10[5,]$intersection, ','))
geneSet6Top10 <- unlist(strsplit(top10[6,]$intersection, ','))
geneSet7Top10 <- unlist(strsplit(top10[7,]$intersection, ','))
geneSet8Top10 <- unlist(strsplit(top10[8,]$intersection, ','))
geneSet9Top10 <- unlist(strsplit(top10[9,]$intersection, ','))
geneSet10Top10 <- unlist(strsplit(top10[10,]$intersection, ','))

geneSet2Random <- sample(rownames(normalizedCounts), 25)

# GSEA for gene-set number 1:
GeneSets1Top10 <- list('top_GO_term' = geneSet1Top10,
                 'random_set' = geneSet2Random)
gseaResults1Top10 <- gage(exprs = log2(normalizedCounts+1),
                    ref = match(rownames(colData[colData$group =='CTRL',]),
                                colnames(normalizedCounts)),
                    samp = match(rownames(colData[colData$group == 'CASE',]),
                                 colnames(normalizedCounts)),
                    gsets = GeneSets1Top10, compare = 'as.group')
gseaResults1Top10

# GSEA for gene-set number 2:
GeneSets2Top10 <- list('top_GO_term' = geneSet2Top10,
                 'random_set' = geneSet2Random)
gseaResults2Top10 <- gage(exprs = log2(normalizedCounts+1),
                    ref = match(rownames(colData[colData$group =='CTRL',]),
                                colnames(normalizedCounts)),
                    samp = match(rownames(colData[colData$group == 'CASE',]),
                                 colnames(normalizedCounts)),
                    gsets = GeneSets2Top10, compare = 'as.group')

# For some reason the results do not contain the 'top_GO_term' analysis, but only for the 'random_set'. 
```
5. What are the other available R packages that can carry out gene set enrichment analysis for RNA-seq datasets? [Difficulty: **Intermediate**]
```
**Answer**: SeqGSEA, clusterProfiler, msigdbr. 
```

6.  Use the topGO package (https://bioconductor.org/packages/release/bioc/html/topGO.html) to re-do the GO term analysis. Compare and contrast the results with what has been obtained using the `gProfileR` package. Which tool is faster, `gProfileR` or topGO? Why? [Difficulty: **Advanced**]
```{r}
library(topGO)
library(GO.db)
library("org.Mm.eg.db")

DEresultsDf1 <- as.data.frame(DEresults)

geneListGO <- DEresultsDf1$pvalue
names(geneListGO) <-  rownames(DEresultsDf1)

# function that returns TRUE/FALSE for p-values<0.05
selection <- function(allScore){ return(allScore < 0.05)} 

# Create topGOData object
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")
TopGOdata <- new("topGOdata",
             ontology="BP",
              allGenes=geneListGO,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
                geneSel=selection)
TopGOdata

# Fisher test:
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(TopGOdata, test.stat)
geneData(resultWeight)

showSigOfNodes(TopGOdata, score(resultWeight), firstSigNodes = 3, useInfo = 'def')


tab <- GenTable(TopGOdata, raw.p.value = resultWeight, topNodes = length(resultWeight@score), numChar = 120)
head(tab)

# `gProfileR` tool is faster compering to topGO. `gProfileR` including one step, and topGO is a multiple step test. 
```

7. Given a gene set annotated for human, how can it be utilized to work on _C. elegans_ data? (Hint: See `biomaRt::getLDS`). [Difficulty: **Advanced**]
```
We have to use biomart library and the function getLDS. In this function there ia an argument
called "values" - this is the place we need to write the gene set of intrest.  

library(biomaRt)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

c.elegans = useMart("ensembl", dataset = "celegans_gene_ensembl") 

getLDS(attributes = c("hgnc_symbol","chromosome_name", "start_position"), 
    filters = "hgnc_symbol", values = **"GENE_SET"**, mart = human, 
    attributesL = c("chromosome_name","start_position"), martL = c.elegans)
```

8. Import curated pathway gene sets with Entrez identifiers from the [MSIGDB database](http://software.broadinstitute.org/gsea/msigdb/collections.jsp) and re-do the GSEA for all curated gene sets. [Difficulty: **Advanced**]
```
** I had a problem when i tried to knit. I recived an error:
"Error in function (type, msg, asError = TRUE)  : 
  Recv failure: Connection was reset"
But at the first time the code did succeeded to run. 
I wrote the code not in a chunk in order to upload the file.**

library("msigdbr")
library(gProfileR)
library(gage)

msigdbr_df <- msigdbr(species = "human", category = "C2")

msigdbr_GeneOfIntrest <- msigdbr_df$gene_symbol

# GO term analysis 
msigdbr_GoResults <- gprofiler(query = msigdbr_GeneOfIntrest,
                       organism = 'hsapiens',
                       src_filter = 'GO',
                       hier_filtering = 'moderate')

# Order go results by pvalue
msigdbr_go <- msigdbr_GoResults[order(msigdbr_GoResults$p.value),]

# When I did this sort I did not get any results, so I choose to skip this step in order to complete the analysis. 
# msigdbr_go <- msigdbr_GoResults[msigdbr_GoResults$overlap.size < 100,]

# Creating a gene sets
msigdbr_geneSet1 <- unlist(strsplit(msigdbr_go[1,]$intersection, ','))
normalizedCounts <- DESeq2::counts(dds, normalized = TRUE)
msigdbr_geneSet2 <- sample(rownames(normalizedCounts), 25)
msigdbr_geneSets <- list('top_GO_term' = msigdbr_geneSet1,
                 'random_set' = msigdbr_geneSet2)
msigdbr_geneSets

#GSEA
msigdbr_GSEAResults <- gage(exprs = log2(normalizedCounts+1),
                    ref = match(rownames(colData[colData$group =='CTRL',]),
                                colnames(normalizedCounts)),
                    samp = match(rownames(colData[colData$group == 'CASE',]),
                                 colnames(normalizedCounts)),
                    gsets = msigdbr_geneSets, compare = 'as.group')
msigdbr_GSEAResults

# Same as question number 4 - I did not get any results for the top_GO_term list. 
```

