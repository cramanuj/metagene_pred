library(cellrangerRkit);library(data.table)

source("~/Research/scripts/r_scripts/plotfns.R")
source("~/Research/scripts/r_scripts/useful_functions.R")

genome <- "GRCh38"
brca = load_cellranger_matrix("~/Research/BRCA/Loi_BCTCells_NatMed_071018/",genome=genome)
dim(exprs(brca))

brca_dat = as.matrix(exprs(brca))
brca_dat = data.frame(brca_dat[,1:5174],check.names=F)

probes = fread("~/Research/pathways/Homo_sapiens.GRCh38.79_table.txt",header=F,data.table=F)
colnames(probes)=c("ENSEMBL","Genes","Coordinates","Strand")
probes$Chr = gsub('\\:.*', '', probes$Coordinates)
probes=probes[c(which(probes$Chr %in% c(1:22)),grep("^X$",probes$Chr),grep("^Y$",probes$Chr)),]
probes = probes[-grep("^MIR",probes$Genes),]
probes = probes[-grep("^RPS",probes$Genes),]
probes = probes[-grep("^RPL",probes$Genes),]

keep_probes = intersect(rownames(brca_dat),probes$ENSEMBL)
seu_exprs = brca_dat[match(keep_probes,rownames(brca_dat)),]
probes_exprs = probes[match(keep_probes,probes$ENSEMBL),]

seu_exprs$Genes = as.factor(as.character(probes_exprs$Genes))
dat = setDT(seu_exprs)[, lapply(.SD, median), by = Genes]
dat = as.data.frame(dat)
rownames(dat) = dat[,1]
dat = dat[,-1]

pdata = data.frame(barcode=colnames(dat))
rownames(pdata) = pdata$barcode

##########
# MONOCLE
##########
library(monocle)

gene_df = data.frame(gene_short_name=rownames(dat))
rownames(gene_df) = rownames(dat)
fd = new("AnnotatedDataFrame", data = gene_df)
mycds = newCellDataSet(as.matrix(dat),phenoData = new("AnnotatedDataFrame", data = pdata), featureData = fd, expressionFamily=negbinomial.size(),lowerDetectionLimit=0.5)
mycds = estimateSizeFactors(mycds)
mycds = estimateDispersions(mycds)
mycds = detectGenes(mycds, min_expr = 0.1)
expressed_genes = row.names(subset(fData(mycds),num_cells_expressed >= 10))
mycds = mycds[expressed_genes,]
x = pData(mycds)$num_genes_expressed
x_1 = (x - mean(x)) / sd(x)
df <- data.frame(x = x_1)
ggplot(df, aes(x)) + geom_histogram(bins = 50) + geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
disp_table <- dispersionTable(mycds)
pData(mycds)$UMI <- Matrix::colSums(exprs(mycds))
ggplot(pData(mycds), aes(num_genes_expressed, UMI)) + geom_point()
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
mycds <- setOrderingFilter(mycds, unsup_clustering_genes$gene_id)
mon_exprs = t(t(exprs(mycds)) /  pData(mycds)[, 'Size_Factor'])
dat_vst = vstExprs(mycds,expr_matrix=mon_exprs)

## ssGSEA

library(GSVA)
imm = read.gmt.file("~/Research/pathways/immune_system.gmt")
gsva_temp = gsva(as.matrix(dat_vst),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Poisson",min.sz=1,max.sz=1000)
rownames(gsva_temp)[20] = c("CD8_CD103")

## Dim Reduction + Density Clustering

FM = gsva_temp

## tSNE
max_components = 2; num_dim = 3; num_clusters = 10; perp = 30
library(Rtsne); library(densityClust); library(plyr); library(dplyr); library(qvalue)
set.seed(123)
irlba_res <- prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1), center = TRUE, scale. = TRUE)
irlba_pca_res <- irlba_res$x
topDim_pca <- irlba_pca_res
set.seed(123)
tsne_res <- Rtsne(as.matrix(topDim_pca), dims = max_components,pca = F,perplexity=perp)
tsne_data <- tsne_res$Y[, 1:max_components]
colnames(tsne_data) = c("Component_1", "Component_2")
rownames(tsne_data) = colnames(FM)
dataDist <- dist(tsne_data)
set.seed(123)
dataClust <- densityClust(dataDist, gaussian = T)
delta_rho_df <- data.frame(delta = dataClust$delta, rho = dataClust$rho)
rho_threshold <- 0
delta_threshold <- sort(delta_rho_df$delta, decreasing = T)[num_clusters] - .Machine$double.eps
dataClust <- densityClust::findClusters(dataClust, rho = rho_threshold, delta = delta_threshold)
tsne_data = data.frame(tsne_data)
tsne_data$Cluster = as.factor(dataClust$cluster)
tcenters = tsne_data[,c(1:3)] %>% dplyr::group_by(Cluster) %>% summarize_all(funs(median))
ggplot(data.frame(tsne_data),aes(x=Component_1,y=Component_2,colour=Cluster,label=Cluster)) + geom_point(size=1) + geom_point(data = tcenters, mapping = aes(x = Component_1, y = Component_2), size = 0, alpha = 0) +
      geom_text(data=tcenters,mapping = aes(label = Cluster), colour="black",size = 6) + theme_pander(boxes=T)

tsne_data1 = cbind(tsne_data,apply(gsva_temp,1,scale))
ave_out = aggregate(. ~ Cluster,tsne_data1,"mean")[,-c(1:3)]
corr_out = cor(ave_out)
corr_out1 = cor(t(ave_out))
centers = tsne_data1[,c(1:3)] %>% dplyr::group_by(Cluster) %>% summarize_all(funs(median))

min_cyt = which.min(aggregate(TCELL_CYTOLYTIC_ACT ~ Cluster,tsne_data1,"mean")[,2])
max_cyt = which.max(aggregate(TCELL_CYTOLYTIC_ACT ~ Cluster,tsne_data1,"mean")[,2])
samples_up = rownames(tsne_data1[tsne_data1$Cluster %in% max_cyt,])
samples_down = rownames(tsne_data1[tsne_data1$Cluster %in% min_cyt,])
samples = c(samples_up,samples_down)

# dat_new = cbind(dat_vst[,match(samples_up,colnames(dat_vst))],dat_vst[,match(samples_down,colnames(dat_vst))])
# dim(dat_new)
# grp1 = length(samples_up); grp2 = length(samples_down)
# pheno_cyt = c(rep(0,grp1),rep(1,grp2))
# dat_new = cbind(dat_vst[,match(samples_up,colnames(dat_vst))],dat_vst[,match(samples_down,colnames(dat_vst))])
# dim(dat_new)
# grp1 = length(samples_up); grp2 = length(samples_down)
# pheno_cyt = c(rep(0,grp1),rep(1,grp2))
# pval_out = laply(1:nrow(dat_new),.progress="time",function(i) wilcox.test(as.numeric(dat_new[i,])~pheno_cyt)$p.val)
# pval_out1 = data.frame("Genes"=rownames(dat_new),"pval"=pval_out,"pval_adj"=p.adjust(pval_out),"qval"=qvalue(pval_out)$qval)
# pval_sig = pval_out1[pval_out1$qval<=0.01,]
# pval_sig = pval_sig[order(pval_sig$pval),]
# pval_sig = na.omit(pval_sig)
# dim(pval_sig)

dat_new = dat[,match(samples,colnames(dat))]
library(Seurat)
pdata = data.frame(barcode=colnames(dat_new))
rownames(pdata) = pdata$barcode
pdata$pheno = pheno_cyt
seed = 123
set.seed(seed)
seu = CreateSeuratObject(raw.data = dat_new, min.cells=3, min.genes=200,is.expr = 0.5, project = "BRCA", meta.data = pdata)
seu = FilterCells(object = seu,subset.names = c("nUMI","nGene"), high.thresholds = c(7500,2500),low.thresholds=c(-Inf,500))
VlnPlot(object=seu,features.plot=c("nGene", "nUMI"),nCol=2)
GenePlot(object = seu, gene1 = "nUMI", gene2 = "nGene")
seu = NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu = FindVariableGenes(seu, do.plot = T, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.015, x.high.cutoff = 3.2, y.cutoff = 0.5)
seu = ScaleData(object = seu)
hist(colSums(seu@data),breaks = 100,main = "Total expression after normalisation",xlab = "Sum of expression")
seu = RunPCA(seu, pcs.print = 0,pc.genes = seu@var.genes)
PCAPlot(object = seu, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seu, dim.1 = 1, dim.2 = 2,group.by="pheno")

PCElbowPlot(seu)
seu = FindClusters(object = seu, reduction.type = "pca", dims.use = 1:10, resolution = 1, print.output = 0, save.SNN = TRUE)
seu = JackStraw(object = seu, num.replicate = 100)
JackStrawPlot(object = seu, PCs = 1:20)
PCHeatmap(object = seu, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE)
seu = RunTSNE(seu, dims.use = 1:15,do.fast=T)
TSNEPlot(seu, do.label = TRUE, pt.size = 0.5)
TSNEPlot(seu, do.label = TRUE, pt.size = 0.5,group.by="pheno")
FeaturePlot(seu, features.plot = c("PRF1","GZMA"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
markers <- FindAllMarkers(object = seu,only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
top10 = markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = seu, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,cex.row = 5,group.cex = 5)

dat_cyt = data.frame(as.matrix(seu@data),check.names=F)
dat_cyt = dat_cyt[match(seu@var.genes,rownames(dat_cyt)),]

########################
### Validation datasets
################################################
## TCGA data

clinical = fread("~/Research/BRCA/tcga_brca_clinical_data.txt",header=TRUE,data.table=F)
dim(clinical)
clinical = clinical[clinical$history_of_neoadjuvant_treatment!=1,]
clinical = clinical[clinical$history_other_malignancy!=1,]
clinical$subtype = as.character(clinical$NEW_SUBTYPE)
clinical$os_year = as.numeric(as.character(clinical$os_year))
clinical$os_Ind = as.numeric(as.character(clinical$os_Ind))
clinical$rfs_year = as.numeric(as.character(clinical$rfs_year))
clinical$rfs_Ind = as.numeric(as.character(clinical$rfs_Ind))
clinical$pfs_year = as.numeric(as.character(clinical$pfs_year))
clinical$pfs_Ind = as.numeric(as.character(clinical$pfs_Ind))
clinical$Age = as.numeric(as.character(clinical$Age))
clinical = clinical[,-c(18:19)]
clinical = clinical[!is.na(clinical$NEW_SUBTYPE),]
samples_clinical = as.character(clinical$bcr_patient_barcode)

datTCGA = fread("~/Research/BRCA/breast_immune/TCGA_BR_Exp.txt",header=T,data.table=F);
rownames(datTCGA) = datTCGA[,1];
datTCGA = datTCGA[,-1];

rownames(clinical)=clinical$bcr_patient_barcode
clinTN = clinical[clinical$subtype=="TN",]
clinTN = clinTN[clinTN$PR != 1,]
clinTN = clinTN[!is.na(clinTN$subtype),]
keep_samples = intersect(clinTN$bcr_patient_barcode,colnames(datTCGA))
tcgaTN = datTCGA[,match(keep_samples,colnames(datTCGA))]
clinTN = clinTN[match(keep_samples,clinTN$bcr_patient_barcode),]

clinER = clinical[clinical$subtype=="ER+",]
keep_ER = intersect(clinER$bcr_patient_barcode,colnames(datTCGA))
tcgaER = datTCGA[,match(keep_ER,colnames(datTCGA))]

## GSVA/ssGSEA
tcgaTN_gsva = gsva(as.matrix(tcgaTN),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(tcgaTN_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/tcga_brca_ssgsea.txt",sep="\t",col.names=NA,quote=F)

tcgaER_gsva = gsva(as.matrix(tcgaER),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(tcgaER_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/tcga_er_brca_ssgsea.txt",sep="\t",col.names=NA,quote=F)

###############################
## METABRIC datasets
## 1) Discovery 2) Validation
###############################
library(impute)
# 1) MB Discovery
metaB = fread("~/Research/Metabric/_ega-box-04_discovery_ExpressionMatrix.txt",header=F,data.table=F)
dim(metaB)
rownames(metaB) = metaB[,1]; metaB = metaB[,-1]
system("head -1 ~/Research/Metabric/_ega-box-04_discovery_ExpressionMatrix.txt > metaB_samples.txt")
metaB_samples = scan("metaB_samples.txt",what="")
colnames(metaB) = metaB_samples

temp = read.delim("~/Research/Metabric/brca_metabric/data_clinical_supp_patient.txt",header=T)
temp = temp[match(colnames(metaB),temp$PATIENT_ID),]
temp1 = read.delim("~/Research/Metabric/brca_metabric/data_clinical_supp_sample.txt",header=T)
temp1 = temp1[match(colnames(metaB),temp1$PATIENT_ID),]
temp1$ER_STATUS = as.character(temp1$ER_STATUS)
temp1$HER2_STATUS = as.character(temp1$HER2_STATUS)
temp1$PR_STATUS = as.character(temp1$PR_STATUS)
temp1$TN = laply(1:nrow(temp1),function(i) length(grep("-",temp1[i,6:8])))
temp2 = cbind(temp,temp1[,-c(1:2)])
clinTN_metaB = temp2[temp2$TN=="3",]
clinTN_metaB$OS_year = clinTN_metaB$OS_MONTHS/12
clinTN_metaB$OS_Ind = ifelse(clinTN_metaB$OS_STATUS=="DECEASED",1,0)
metaB_dss = read.csv("~/Research/Metabric/brca_metabric/Clinical_Overall_Survival_Data_from_METABRIC.txt",header=T)
metaB_dss1 = metaB_dss[match(clinTN_metaB$PATIENT_ID,metaB_dss$SampleID),]
clinTN_metaB$DSS_days = metaB_dss1$time
clinTN_metaB$DSS_event = metaB_dss1$status
metaB_TN = metaB[,match(clinTN_metaB$PATIENT_ID,colnames(metaB))]
annot = read.delim("~/Illumina/HumanHT-12_V4_0_R2_15002873_B.txt",header=T)
annot = annot[,c(14,5,12)]
annot$Symbol = as.character(annot$Symbol)
illumina_probes = intersect(rownames(metaB_TN),annot$Probe_Id)
annot = annot[match(illumina_probes,annot$Probe_Id),]
metaB_TN = metaB_TN[match(illumina_probes,rownames(metaB_TN)),]
metaB_TN$Genes = as.character(annot$ILMN_Gene)
metaB_TN = setDT(metaB_TN)[, lapply(.SD, mean), by = Genes]
metaB_TN = as.data.frame(metaB_TN)
rownames(metaB_TN) = metaB_TN[,1]
metaB_TN = metaB_TN[,-1]

## GSVA/ssGSEA
metaB_TN_gsva = gsva(as.matrix(metaB_TN),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(metaB_TN_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/metaB_ssgsea.txt",sep="\t",col.names=NA,quote=F)

### 2) MB Validation
metaBV = fread("~/Research/Metabric/_ega-box-04_validation_ExpressionMatrix.txt",header=F,sep=" ",data.table=F)
rownames(metaBV) = metaBV[,1]; metaBV = metaBV[,-1]
system("head -1 ~/Research/Metabric/_ega-box-04_validation_ExpressionMatrix.txt > metaBV_samples.txt")
metaBV_samples = scan("metaBV_samples.txt",what="")
colnames(metaBV) = metaBV_samples

temp = read.delim("~/Research/Metabric/brca_metabric/data_clinical_supp_patient.txt",header=T)
keep_samples = intersect(colnames(metaBV),temp$PATIENT_ID)
temp = temp[match(keep_samples,temp$PATIENT_ID),]
metaBV = metaBV[,match(keep_samples,colnames(metaBV))]
temp1 = read.delim("~/Research/Metabric/brca_metabric/data_clinical_supp_sample.txt",header=T)
temp1 = temp1[match(colnames(metaBV),temp1$PATIENT_ID),]
temp1$ER_STATUS = as.character(temp1$ER_STATUS)
temp1$HER2_STATUS = as.character(temp1$HER2_STATUS)
temp1$PR_STATUS = as.character(temp1$PR_STATUS)
temp1$TN = laply(1:nrow(temp1),function(i) length(grep("-",temp1[i,6:8])))
temp2 = cbind(temp,temp1[,-c(1:2)])
clinTN_metaBV = temp2[temp2$TN=="3",]
clinTN_metaBV$OS_year = clinTN_metaBV$OS_MONTHS/12
clinTN_metaBV$OS_Ind = ifelse(clinTN_metaBV$OS_STATUS=="DECEASED",1,0)
metaBV_dss = read.csv("~/Research/Metabric/brca_metabric/Clinical_Overall_Survival_Data_from_METABRIC.txt",header=T)
metaBV_dss1 = metaBV_dss[match(clinTN_metaBV$PATIENT_ID,metaBV_dss$SampleID),]
clinTN_metaBV$DSS_days = metaBV_dss1$time
clinTN_metaBV$DSS_event = metaBV_dss1$status

metaBV_TN = metaBV[,match(clinTN_metaBV$PATIENT_ID,colnames(metaBV))]
metaBV_TN = metaBV_TN[match(illumina_probes,rownames(metaBV_TN)),]
metaBV_TN$Genes = as.character(annot$ILMN_Gene)
metaBV_TN = setDT(metaBV_TN)[, lapply(.SD, mean), by = Genes]
metaBV_TN = as.data.frame(metaBV_TN)
rownames(metaBV_TN) = metaBV_TN[,1]
metaBV_TN = metaBV_TN[,-1]
temp_imp = impute.knn(as.matrix(metaBV_TN))
metaBV_TN = temp_imp$data

### ssGSEA
metaBV_TN_gsva = gsva(as.matrix(metaBV_TN),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(metaBV_TN_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/metaBV_ssgsea.txt",sep="\t",col.names=NA,quote=F)

##########
## TCGA SKCM/Melanoma
#######################

tcgaMELA = fread("~/Research/SKCM/SKCM_normalized.txt",header=T,data.table=F)
rownames(tcgaMELA) = tcgaMELA[,1]
tcgaMELA = tcgaMELA[,-1]

# ssGSEA
tcga_mela_gsva = gsva(as.matrix(tcgaMELA),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(tcga_mela_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/tcga_mela_ssgsea.txt",sep="\t",col.names=NA,quote=F)

tcga_clinMELA = read.delim("~/Research/SKCM/SKCM_clinical.txt",header=T,row.names=1)
dim(tcga_clinMELA)

######
# LUAD
################

tcgaLUAD = fread("~/Research/LUAD/LUAD_normalized.txt",header=T,data.table=F)
rownames(tcgaLUAD) = tcgaLUAD[,1]
tcgaLUAD = tcgaLUAD[,-1]

# ssGSEA
tcga_luad_gsva = gsva(as.matrix(tcgaLUAD),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(tcga_luad_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/tcga_luad_ssgsea.txt",sep="\t",col.names=NA,quote=F)

tcga_clinLUAD = read.delim("~/Research/LUAD/LUAD_clinical.txt",header=T,row.names=1)
dim(tcga_clinLUAD)

######
# LUSC
################

tcgaLUSC = fread("~/Research/LUSC/LUSC_normalized.txt",header=T,data.table=F)
rownames(tcgaLUSC) = tcgaLUSC[,1]
tcgaLUSC = tcgaLUSC[,-1]

# ssGSEA
tcga_lusc_gsva = gsva(as.matrix(tcgaLUSC),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(tcga_lusc_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/tcga_lusc_ssgsea.txt",sep="\t",col.names=NA,quote=F)

tcga_clinLUSC = read.delim("~/Research/LUSC/LUSC_clinical.txt",header=T,row.names=1)
dim(tcga_clinLUSC)

######
# HNSC
################

tcgaHNSC = fread("~/Research/HNSC/HNSC_normalized.txt",header=T,data.table=F)
rownames(tcgaHNSC) = tcgaHNSC[,1]
tcgaHNSC = tcgaHNSC[,-1]

# ssGSEA
tcga_hnsc_gsva = gsva(as.matrix(tcgaHNSC),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(tcga_hnsc_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/tcga_hnsc_ssgsea.txt",sep="\t",col.names=NA,quote=F)

tcga_clinHNSC = read.delim("~/Research/HNSC/HNSC_clinical.txt",header=T,row.names=1)
dim(tcga_clinHNSC)

######
# COAD
################

tcgaCOAD = fread("~/Research/COAD/COAD_normalized.txt",header=T,data.table=F)
rownames(tcgaCOAD) = tcgaCOAD[,1]
tcgaCOAD = tcgaCOAD[,-1]

# ssGSEA
tcga_coad_gsva = gsva(as.matrix(tcgaCOAD),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(tcga_coad_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/tcga_coad_ssgsea.txt",sep="\t",col.names=NA,quote=F)

tcga_clinCOAD = read.delim("~/Research/COAD/COAD_clinical.txt",header=T,row.names=1)
dim(tcga_clinCOAD)

######
# KIRC
################

tcgaKIRC = fread("~/Research/KIRC/KIRC_normalized.txt",header=T,data.table=F)
rownames(tcgaKIRC) = tcgaKIRC[,1]
tcgaKIRC = tcgaKIRC[,-1]

# ssGSEA
tcga_kirc_gsva = gsva(as.matrix(tcgaKIRC),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(tcga_kirc_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/tcga_kirc_ssgsea.txt",sep="\t",col.names=NA,quote=F)

tcga_clinKIRC = read.delim("~/Research/KIRC/KIRC_clinical.txt",header=T,row.names=1)
dim(tcga_clinKIRC)

#######
## GSE78220
## Pre-treatment melanomas undergoing anti-PD-1 checkpoint inhibition therapy
#######################

clinMELA = read.delim("~/Research/antiPD1_data/GSE78220/sampleInfo.txt",header=T)
clinMELA$OS_yr = clinMELA$Overall_Survival/365
clinMELA$OS_ind = clinMELA$Patient_Status
clinMELA$OS_ind = gsub("Dead","1",clinMELA$OS_ind)
clinMELA$OS_ind = gsub("Alive","0",clinMELA$OS_ind)
clinMELA$OS_ind = as.numeric(clinMELA$OS_ind)

mela = fread("~/Research/antiPD1_data/GSE78220/GSE78220_normalized.txt",header=T,data.table=F)
rownames(mela) = mela[,1]; mela = mela[,-1]
#
### ssGSEA
mela_gsva = gsva(as.matrix(mela),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(mela_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/mela_ssgsea.txt",sep="\t",col.names=NA,quote=F)

###########
## GSE91061
## Immune checkpoint blockade therapy
##############################################

gse91061 = read.delim("~/Research/antiPD1_data/GSE91061/GSE91061_normalized.txt",header=T,row.names=1,check.names=F)
dim(gse91061)
geneid = read.delim("~/Research/antiPD1_data/GSE91061/gene_id_conversion.txt",header=T)
gse91061$Genes = geneid$Gene_Symbols
gse91061 = setDT(gse91061)[, lapply(.SD, median), by = Genes]
gse91061 = data.frame(gse91061,check.names=F)
gse91061 = gse91061[!is.na(gse91061$Genes),]
rownames(gse91061) = gse91061$Genes
gse91061 = gse91061[,-1]

gse91061_CLIN = read.delim("~/Research/antiPD1_data/GSE91061/GSE91061_fullClinical.txt",header=T,row.names=1,check.names=F)
dim(gse91061_CLIN)
gse91061_CLIN$VISIT = as.character(gse91061_CLIN$VISIT)
gse91061_CLIN$VISIT = gsub(" ","",gse91061_CLIN$VISIT)
gse91061_CLIN_pre = gse91061_CLIN[gse91061_CLIN$VISIT=="Pre",]
gse91061_CLIN_pre$SAMPLE_TITLE = as.character(gse91061_CLIN_pre$SAMPLE_TITLE)
gse91061_pre = gse91061[,match(gse91061_CLIN_pre$SAMPLE_TITLE,colnames(gse91061))]

### ssGSEA
gse_gsva = gsva(as.matrix(gse91061_pre),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(gse_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/GSE91061_ssgsea.txt",sep="\t",col.names=NA,quote=F)

#######
## GSE78220
## Pre-treatment melanomas undergoing anti-PD-1 checkpoint inhibition therapy
#######################

clinMELA = read.delim("~/Research/antiPD1_data/GSE78220/sampleInfo.txt",header=T)
clinMELA$OS_yr = clinMELA$Overall_Survival/365
clinMELA$OS_ind = clinMELA$Patient_Status
clinMELA$OS_ind = gsub("Dead","1",clinMELA$OS_ind)
clinMELA$OS_ind = gsub("Alive","0",clinMELA$OS_ind)
clinMELA$OS_ind = as.numeric(clinMELA$OS_ind)

gse78220 = fread("~/Research/antiPD1_data/GSE78220/GSE78220_normalized.txt",header=T,data.table=F)
rownames(gse78220) = gse78220[,1]; gse78220 = gse78220[,-1]
#
### ssGSEA
mela_gsva = gsva(as.matrix(gse78220),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
write.table(mela_gsva,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/mela_ssgsea.txt",sep="\t",col.names=NA,quote=F)

save.image("metagene_predictor_newer.rds")
