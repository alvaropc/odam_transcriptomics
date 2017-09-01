####
#Microarrays Practise
#Álvaro Ponce Cabrera
#####################
## Published article Blood coagulation protein fibrinogen promotes autoimmunity and demyelination 
## via chemokine release and antigen presentation has associated GEO Series GSE71084 that contains 
## two datasets, one for mouse and one for rat, of fibrin(ogen) stimulation in two different tissues. 
## Download the data and analyze it following the steps and methods seen in the course:

## First we have download the data from GEO (Gene Expression Omnibus).
## Packages GEOquery and oligo will be used for doing that:

library(GEOquery) 
library(oligo)
getwd()
setwd('/home/biocloud')
gseid<-"GSE71084"


gse <- getGEO(gseid,GSEMatrix=TRUE) 

GEO<-getGEOSuppFiles(gseid,makeDirectory=TRUE) #to download raw data and all other 

## We get the uncompressed files. Then we read the data:
untar(file.path(getwd(),"GSE71084_RAW.tar"), exdir = file.path(getwd(),gseid))
setwd(file.path(getwd(),"GSE71084"))

## We have to read them separetely because the arrays of mouse and rat are different.
## We use for that the function read.celfiles (that is in package oligo).
celFiles.mouse = list.files( pattern = "Mouse") #We get Cel files we have in the folder wich Mouse.
GEOFS.mouse <- read.celfiles(celFiles.mouse) #Read the files

celFiles.Rat = list.files( pattern = "Rat") #Same but rat

GEOFS.rat <- read.celfiles(celFiles.Rat)


#Exploration of the data

GEOFS.mouse
GEOFS.rat
class(GEOFS.rat)
class(GEOFS.mouse)


### 1. Quality assessment.

## We can see the arrays with function image to confirm that there is not
## physical damage:

image(GEOFS.mouse)
image(GEOFS.rat)
#Everything ok


## Data exploration and quality assessment:

colos<-rainbow(4)


## Density plots:

hist (GEOFS.mouse, col=colos)
hist (GEOFS.rat, col=colos)

## Box plots: 

boxplot (GEOFS.mouse, col=colos, las=3)
boxplot (GEOFS.rat, col=colos, las=3)

## MA plots:

MAplot (GEOFS.mouse)
MAplot (GEOFS.mouse, pairs=TRUE)
MAplot (GEOFS.rat)
MAplot (GEOFS.rat, pairs=TRUE)

## All of the values are near 0, that is the ideal situation.

## RLE and NUSE for mouse:

GEOFS.fit.mouse <- fitProbeLevelModel (GEOFS.mouse)
RLE (GEOFS.fit.mouse, col=colos)
NUSE (GEOFS.fit.mouse, col=colos)

## RLE is near 0 and NUSE is near 1, so everything is ok.

## RLE and NUSE for rat:

GEOFS.fit.rat <- fitProbeLevelModel (GEOFS.rat)
RLE (GEOFS.fit.rat, col=colos)
NUSE (GEOFS.fit.rat, col=colos)

## Here, we also see that RLE is near 0 and NUSE is near 1.


## QC report

#We can obtain a QC report with arrayQuaityMetrics package

library(arrayQualityMetrics)
arrayQualityMetrics(GEOFS.mouse, force = TRUE, do.logtransform = TRUE)
arrayQualityMetrics(GEOFS.rat, force = TRUE, do.logtransform = TRUE)



###2. Normalization.

#Function rma from oligo package will be use:

genePS.rat <- rma(GEOFS.rat)
slotNames(genePS.rat)
genePS.mouse <- rma(GEOFS.mouse )
slotNames(genePS.mouse)

## Now we can see how the normalization went:

## Density plots:

hist (genePS.rat, col=colos)
hist (genePS.mouse, col=colos)

## Box plots: 

boxplot (genePS.rat, col=colos, las=3)
boxplot (genePS.mouse, col=colos, las=3)



## We also can see the data aggregation in order to check for Batch effect:

use.cor="pairwise.complete.obs"

## Mouse:

x.mouse <- exprs (genePS.mouse)
clust.cor.ward.mouse <- hclust (as.dist (1-cor (x.mouse, use=use.cor)), method="ward.D2")
plot (clust.cor.ward.mouse, main="Hierarchical clustering mouse", hang=-1,cex=0.6)

## The two conditions are grouped, so it is ok.

## Rat:

x.rat <- exprs (genePS.rat)
clust.cor.ward.rat <- hclust (as.dist (1-cor (x.rat, use=use.cor)), method="ward.D2")
plot (clust.cor.ward.rat, main="Hierarchical clustering rat", hang=-1,cex=0.6)

## Same as for mouse.



### 3. Detection of DEG, justify the criteria for the selection of DEG.

## Package limma will be used
library(limma)

## For rat
cond<-as.factor(c(rep("control",2),rep("fibrin",2)))

#Design matrix
design.rat<-model.matrix(~0+cond)

colnames(design.rat)<-gsub("cond","",colnames(design))
rownames(design.rat)<-sampleNames(genePS.rat)
head(design.rat)


#We adjust the model to the matrix using the function lmFit.
#With makeContrats we specify the comparison. 
#With eBayes we obtain the bayesian adjustment. 

fit.rat<-lmFit(genePS.rat,design.rat) #model could also be fitted to the batch corrected data


contrast.matrix.rat<-makeContrasts(fibrin-control,levels=c("control","fibrin"))
head(contrast.matrix.rat)
fit2.rat<-contrasts.fit(fit.rat,contrast.matrix.rat)
fite.rat<-eBayes(fit2.rat)


#With topTable we get the logFC, P.Value and adj.P.Val among other values.

top.table.rat<-topTable(fite.rat,coef=1,number=Inf,adjust="BH",sort.by="logFC")
head(top.table.rat)

#Creation of a summary table of results
Results_summary_rat = data.frame(Adjusted = table(decideTests(fite.rat)),
                               Non_Adjusted = table(decideTests(fite.rat, adjust.method="none")))

Results_summary_rat


## We can try to obtain the DEGs by using the adjusted p-value (previously
## adjusted by Benjamini & Hochberg approach). This strategy is the most
## accurate in a statistical way and will be the choosen one for this reason.
results.adj.p.rat <- top.table.rat[top.table.rat$adj.P.Val < 0.05,] # Subsetting the significance threshold
head(results.adj.p.rat)
dim(results.adj.p.rat)


## Mouse

#Design matrix
design.mouse<-model.matrix(~0+cond)
colnames(design.mouse)<-gsub("cond","",colnames(design.mouse))
rownames(design.mouse)<-sampleNames(genePS.mouse)
fit.mouse<-lmFit(genePS.mouse,design.mouse) 
contrast.matrix.mouse<-makeContrasts(fibrin-control,levels=c("control","fibrin"))
head(contrast.matrix.mouse)
fit2.mouse<-contrasts.fit(fit.mouse,contrast.matrix.mouse)
fite.mouse<-eBayes(fit2.mouse)


top.table.mouse<-topTable(fite,coef=1,number=Inf,adjust="BH",sort.by="logFC")
head(top.table.mouse)
results.adj.p.mouse <- top.table.mouse[top.table.mouse$adj.P.Val < 0.05,] # Subsetting the significance threshold
head(results.adj.p.mouse)

## Since we do not find any results we can try to use just the p-value. 

results.p.mouse <- top.table.mouse[top.table.mouse$P.Value < 0.05,]
head(results.p.mouse)
dim(results.p.mouse)

Results_summary.mouse = data.frame(Adjusted = table(decideTests(fite.mouse)),
                                   Non_Adjusted = table(decideTests(fite.mouse, adjust.method="none")))



### 4. Generate a Volcano plot and a heat map of results.

## We will use for that the package limma and its function volcanoplot and
## the package gplots and its function heatmap.2.

library (gplots)

## Volcano plot for mouse:

volcanoplot (fite.mouse, coef=1, highlight=10, names=fite.mouse$genes$NAME, 
             main="Volcano plot mouse")

## Volcano plot for rat:

volcanoplot (fite.rat ,coef=1, highlight=10, names=fite.rat$genes$NAME, 
             main="Volcano plot rat")



## Heat map for mouse:

data.clus.mouse <- exprs (genePS.mouse [rownames (results.p.mouse),])
clust.rows <- hclust (as.dist (1-cor (t (data.clus.mouse))), method="ward.D2")
clust.cols <- hclust (as.dist (1-cor (data.clus.mouse)), method="average")  
heatcol <- colorRampPalette (c("green", "Black","red"), space = "rgb")

heatm.mouse <- heatmap.2 (as.matrix (data.clus.mouse), col = heatcol(256), 
                    dendrogram="column", Colv=as.dendrogram (clust.cols), 
                    Rowv=as.dendrogram (clust.rows), 
                    scale="row", cexRow=0.1, cexCol=0.5,  
                    main="", key=TRUE, keysize=1, density.info="none", trace="none")

## Heat map for rat:

data.clus.rat <- exprs(genePS.rat[rownames(results.adj.p.rat),])
clust.rows <- hclust(as.dist (1-cor (t (data.clus.rat))) ,method="ward.D2")
clust.cols <- hclust(as.dist (1-cor (data.clus.rat)), method="average")  

heatm.rat <- heatmap.2 (as.matrix (data.clus.rat), col = heatcol(256), 
                    dendrogram="column", Colv=as.dendrogram (clust.cols), 
                    Rowv=as.dendrogram (clust.rows), 
                    scale="row", cexRow=0.1, cexCol=0.5, 
                    main="", key=TRUE, keysize=1, density.info="none", trace="none")

## 5. Annotation of results

# 
## We will use the package annotate for that and then the one specific for each organism. 

library (annotate)

## For the mouse:

library (mogene10sttranscriptcluster.db)

## We get all the values we are interested in to build the annotation:

dat.mouse <- exprs (genePS.mouse) [rownames (results.p.mouse),]
logFC.mouse <- results.p.mouse$logFC
pval.mouse <- results.p.mouse$P.Value
adj.pval.mouse <- results.p.mouse$adj.P.Val

sym.mouse <- unlist (mget (rownames (results.p.mouse), env=mogene10sttranscriptclusterSYMBOL))
name.mouse <- unlist (mget (rownames (results.p.mouse), env=mogene10sttranscriptclusterGENENAME))
chr.mouse <- unlist (mget (rownames (results.p.mouse), env=mogene10sttranscriptclusterCHR))
length (chr.mouse) #2285. It must match the 2283 DEGs we got.
chr.mouse.1 <- match (rownames (results.p.mouse), names(chr.mouse))
chr.mouse.2 <- chr.mouse[chr.mouse.1]
length (chr.mouse.2) #2283. Now it is well. 

## We decide the characteristics that our report file will have, such as its name:

affyids.mouse <- rownames (results.p.mouse)
genelist.mouse <- list (affyids.mouse)
filename.mouse <- "Results.mouse.html"
title.mouse <- "Differentially expressed Mouse. Fibrinogen vs control"
othernames.mouse <- list (sym.mouse, name.mouse, chr.mouse.2, round (logFC.mouse, 1), 
                          round(pval.mouse, 4), round (adj.pval.mouse, 4), round(dat.mouse, 2)) 
head.mouse <- c("Probe ID", "Symbol", "Gene Name", "Chr","logFC", "p-value","adj.p-value",sampleNames(genePS.mouse))
repository.mouse <- list("affy")

## We finally create our .html report with the DEGs. 

htmlpage (genelist.mouse, filename.mouse, title.mouse, othernames.mouse, head.mouse, repository = repository.mouse)

## For the rat:

library(ragene10sttranscriptcluster.db)

## We get all the values we are interested in to build the annotation:

dat.rat <- exprs (genePS.rat) [rownames (results.adj.p.rat),]
logFC.rat <- results.adj.p.rat$logFC
pval.rat <- results.adj.p.rat$P.Value
adj.pval.rat <- results.adj.p.rat$adj.P.Val

sym.rat <- unlist (mget (rownames (results.adj.p.rat), env=ragene10sttranscriptclusterSYMBOL))
name.rat <- unlist (mget (rownames (results.adj.p.rat), env=ragene10sttranscriptclusterGENENAME))
chr.rat <- unlist (mget (rownames (results.adj.p.rat), env=ragene10sttranscriptclusterCHR))

## We decide the characteristics that our report file will have, such as its name:

affyids.rat <- rownames (results.adj.p.rat)
genelist.rat <- list (affyids.rat)
filename.rat <- "Results.rat.html"
title.rat <- "Differentially expressed Rat. Fibrin vs control"
othernames.rat <- list(sym.rat ,name.rat ,chr.rat, round(logFC.rat, 1), 
                       round(pval.rat, 4), round(adj.pval.rat, 4), round(dat.rat, 2)) 
head.rat <- c("Probe ID", "Symbol", "Gene Name", "Chr","logFC", "p-value","adj.p-value", sampleNames(genePS.rat))
repository.rat <- list("affy")

## We finally create our .html report with the DEGs. 

htmlpage (genelist.rat, filename.rat, title.rat, othernames.rat, head.rat, repository = repository.rat)

### 6. Functional analysis using GOstats.

library (GOstats)
library (GSEABase)

## For the mouse:

library (org.Mm.eg.db)

frame = toTable (org.Mm.egGO)
goframeData = data.frame (frame$go_id, frame$Evidence, frame$gene_id)
head (goframeData)

goFrame=GOFrame (goframeData, organism="Mus musculus")
goAllFrame=GOAllFrame (goFrame)

gsc.GO <- GeneSetCollection (goAllFrame, setType = GOCollection())

## We will use entrez ids.

entrez <- unlist (mget (featureNames (genePS.mouse), 
                        env=mogene10sttranscriptclusterENTREZID))
length (entrez)

## Our universe:

length (unique (entrez))
entrezuniverse <- unique (entrez)

## Results: 

results.entrez <- unique (unlist (mget (rownames (results.p.mouse), 
                                        env=mogene10sttranscriptclusterENTREZID)))

## For biological processes:

cutoff=0.001
GO.params.bp <- GSEAGOHyperGParams (name="GOstats",  geneSetCollection=gsc.GO, 
                                    geneIds = results.entrez, universeGeneIds = entrezuniverse,  
                                    ontology = "BP", pvalueCutoff = cutoff,  conditional = FALSE, 
                                    testDirection = "over")

## We create our report in .html:

GO.results.bp <- hyperGTest (GO.params.bp)
head (summary (GO.results.bp))
GO.results.bp
htmlReport (GO.results.bp, "GO.results.BiologicalProcessMouse.html")

## We now, chanching the argument ontology of the function GSEAGOHyperGParams,
## can obtain the other two ontologies, that are molecular function (MF) and
## cellular component (CC).So we do:

## For molecular function:

GO.params.bp <- GSEAGOHyperGParams (name="GOstats",  geneSetCollection=gsc.GO, 
                                    geneIds = results.entrez, universeGeneIds = entrezuniverse,  
                                    ontology = "MF", pvalueCutoff = cutoff,  conditional = FALSE, 
                                    testDirection = "over")

GO.results.bp <- hyperGTest (GO.params.bp)
head (summary (GO.results.bp))
GO.results.bp
htmlReport (GO.results.bp, "GO.results.MolecularFunctionMouse.html")

## For cellular components:

GO.params.bp <- GSEAGOHyperGParams (name="GOstats",  geneSetCollection=gsc.GO, 
                                    geneIds = results.entrez, universeGeneIds = entrezuniverse,  
                                    ontology = "CC", pvalueCutoff = cutoff,  conditional = FALSE, 
                                    testDirection = "over")

GO.results.bp <- hyperGTest (GO.params.bp)
head (summary (GO.results.bp))
GO.results.bp
htmlReport (GO.results.bp, "GO.results.CellularComponentMouse.html")

## For the rat:

library (org.Rn.eg.db)

frame = toTable (org.Rn.egGO)
goframeData = data.frame (frame$go_id, frame$Evidence, frame$gene_id)
head (goframeData)

goFrame=GOFrame (goframeData, organism="Rattus norvegicus")
goAllFrame=GOAllFrame (goFrame)

gsc.GO <- GeneSetCollection (goAllFrame, setType = GOCollection())

## We will use entrez ids.

entrez <- unlist (mget (featureNames (GEOFS.rma.rat), 
                        env=ragene10sttranscriptclusterENTREZID))
length (entrez)

## Our universe:

length (unique (entrez))
entrezuniverse <- unique (entrez)

## Results:

results.entrez <- unique (unlist (mget (rownames (results.adj.p.rat), 
                                        env=ragene10sttranscriptclusterENTREZID)))

## For biological processes:

cutoff=0.001
GO.params.bp <- GSEAGOHyperGParams(name="GOstats",  geneSetCollection=gsc.GO, 
                                   geneIds = results.entrez, universeGeneIds = entrezuniverse,  
                                   ontology = "BP", pvalueCutoff = cutoff,  conditional = FALSE, 
                                   testDirection = "over")

## We create our report in .html:

GO.results.bp<-hyperGTest(GO.params.bp)
head(summary(GO.results.bp))
GO.results.bp
htmlReport(GO.results.bp, "GO.results.BiologicalProcessRat.html")

## For molecular function:

GO.params.bp <- GSEAGOHyperGParams (name="GOstats",  geneSetCollection=gsc.GO, 
                                    geneIds = results.entrez, universeGeneIds = entrezuniverse,  
                                    ontology = "MF", pvalueCutoff = cutoff,  conditional = FALSE, 
                                    testDirection = "over")

GO.results.bp <- hyperGTest (GO.params.bp)
head (summary (GO.results.bp))
GO.results.bp
htmlReport (GO.results.bp, "GO.results.MolecularFunctionRat.html")

## For cellular components:

GO.params.bp <- GSEAGOHyperGParams (name="GOstats",  geneSetCollection=gsc.GO, 
                                    geneIds = results.entrez, universeGeneIds = entrezuniverse,  
                                    ontology = "CC", pvalueCutoff = cutoff,  conditional = FALSE, 
                                    testDirection = "over")

GO.results.bp <- hyperGTest (GO.params.bp)
head (summary (GO.results.bp))
GO.results.bp
htmlReport (GO.results.bp, "GO.results.CellularComponentRat.html")

## 7. Compare homolog results obtained for Mouse and rat, despite the different tissues and 
## conditions, using biomaRt.

library (biomaRt)

listMarts(host = "www.ensembl.org")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
listDatasets(mart)[1:69,]

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="www.ensembl.org")
listAttributes(mart)[1:10,]
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="rnorvegicus_gene_ensembl", host="www.ensembl.org")
listAttributes(mart)[1:180,]

listFilters(mart)

mouse <- useMart (biomart="ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",
                  host = "www.ensembl.org")

rat <- useMart (biomart="ENSEMBL_MART_ENSEMBL", dataset = "rnorvegicus_gene_ensembl",
                host = "www.ensembl.org")

## We use the function getLDS to obtain the information that we want (the attributes)
## to compare the homolog results.

my.results <- getLDS (attributes = c("ensembl_gene_id", "chromosome_name"),
                      #filters = "", 
                      values = sym.rat, 
                      mart = rat,
                      attributesL = c("ensembl_gene_id", "chromosome_name"),
                      martL = mouse)

head (my.results)