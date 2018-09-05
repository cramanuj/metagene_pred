library(monocle); library(cellrangerRkit); library(ggplot2); library(cowplot); library(data.table);
library(GSA); library(GSVA); library(ggpubr); library(ggthemes); library(car); library(qvalue);
library(data.table); library(limma); library(sva); library(randomForest)
library(FactoMineR);library(factoextra); library(survival); library(Seurat);
library(survival); library(survminer); library(maftools);
library(gplots); library(igraph); library(irlba); library(Rtsne); library(densityClust);
library(plyr); library(dplyr); library(caret); library(broom); library(impute);
library(ggcorrplot); library(genefu)

source("~/Research/scripts/r_scripts/plotfns.R")
source("~/Research/scripts/r_scripts/useful_functions.R")

ggplotRegression <- function(fit){
ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     " P =",signif(summary(fit)$coef[2,4], 5))) + xlab("Predicted Probability") + ylab("T-cell Cytolytic Activity")
}

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
fwrite(dat,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/BRCA_scRNAseq_gene_level.txt",sep="\t",col.names=T,quote=F,row.names=T)

pdata = data.frame(barcode=colnames(dat))
rownames(pdata) = pdata$barcode

##########
# MONOCLE
##########
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
plot_ordering_genes(mycds)
plot_pc_variance_explained(mycds, return_all = FALSE)
mycds = reduceDimension(mycds, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = T)
mycds = clusterCells(mycds, num_clusters = 10)
plot_cell_clusters(mycds, 1, 2, color = "Cluster")
mycds = reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds = orderCells(mycds)
plot_cell_trajectory(mycds, color_by = "Cluster")

mon_exprs = t(t(exprs(mycds)) /  pData(mycds)[, 'Size_Factor'])
dat_vst = vstExprs(mycds,expr_matrix=mon_exprs)
fwrite(dat_vst,"BRCA_scRNAseq_VST.txt",sep="\t",col.names=T,row.names=T,quote=F,)
imm = read.gmt.file("~/Research/pathways/immune_system.gmt")
gsva_temp = gsva(as.matrix(dat_vst),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Poisson",min.sz=1,max.sz=1000)
rownames(gsva_temp)[20] = c("CD8_CD103")
write.table(gsva_temp,"BRCA_ssRNAseq_ssgsea.txt",sep="\t",col.names=NA,quote=F)

##

FM = gsva_temp

## tSNE
max_components = 2; num_dim = 3; num_clusters = 10; perp = 30

set.seed(123)
irlba_res <- prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1), center = TRUE, scale. = TRUE)
irlba_pca_res <- irlba_res$x
topDim_pca <- irlba_pca_res
set.seed(123)
tsne_res <- Rtsne(as.matrix(topDim_pca), dims = max_components,pca = F,perplexity=perp)
tsne_data <- tsne_res$Y[, 1:max_components]
colnames(tsne_data) = c("Component_1", "Component_2")
rownames(tsne_data) = colnames(FM)

### Clustering

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
ggsave("tSNE_plot_ssGSEA.pdf",width=10.9,height=6.71)


tsne_data1 = cbind(tsne_data,apply(gsva_temp,1,scale))
ave_out = aggregate(. ~ Cluster,tsne_data1,"mean")[,-c(1:3)]
corr_out = cor(ave_out)
ggcorrplot(corr_out, hc.order = TRUE,type = "lower", lab = TRUE, lab_size = 2, method="circle", colors = c("tomato2", "white", "springgreen3"),
           title="Correlogram of immune metagenes", tl.cex = 6, ggtheme=theme_bw)
ggsave("corr_plot_ssGSEA.pdf",width=10.9,height=6.71)

corr_out1 = cor(t(ave_out))
ggcorrplot(corr_out1, hc.order = TRUE,type = "lower", lab = TRUE, lab_size = 2, method="circle", colors = c("tomato2", "white", "springgreen3"),
           title="Correlogram of clusters", tl.cex = 6, ggtheme=theme_bw)

centers = tsne_data1[,c(1:3)] %>% dplyr::group_by(Cluster) %>% summarize_all(funs(median))

pdf("BRCA_boxplots_ssgsea.pdf")
for(i in 1:nrow(gsva_temp)){
  print(i)
    print(ggboxplot(tsne_data1,x="Cluster",y=rownames(gsva_temp)[i],fill="grey"))
}
dev.off()

print(ggplot(tsne_data1,aes(x=Component_1,y=Component_2)) + geom_point(colour="light grey",size=1) +
      geom_point(data=tsne_data1[tsne_data1$TCELL_CYTOLYTIC_ACT > quantile(tsne_data1$TCELL_CYTOLYTIC_ACT,0.9),],aes(x=Component_1, y=Component_2), colour="blue", size=1) +
      geom_point(data = centers, mapping = aes(x = Component_1, y = Component_2), size = 0, alpha = 0) + geom_text(data=centers,mapping = aes(label = Cluster), colour="black",size = 5)
    )
ggsave("tSNE_plot_CYT_ssGSEA.pdf",width=10.9,height=6.71)

print(ggplot(tsne_data1,aes(x=Component_1,y=Component_2)) + geom_point(colour="light grey",size=1) +
      geom_point(data=tsne_data1[tsne_data1$MACROPHAGE_ACTIVITY > quantile(tsne_data1$MACROPHAGE_ACTIVITY,0.95),],aes(x=Component_1, y=Component_2), colour="blue", size=1) +
      geom_point(data = centers, mapping = aes(x = Component_1, y = Component_2), size = 0, alpha = 0) + geom_text(data=centers,mapping = aes(label = Cluster), colour="black",size = 4)
    )

print(ggplot(tsne_data1,aes(x=Component_1,y=Component_2)) + geom_point(colour="light grey",size=1) +
      geom_point(data=tsne_data1[tsne_data1$MHC_II > quantile(tsne_data1$MHC_II,0.9),],aes(x=Component_1, y=Component_2), colour="blue", size=1) + geom_label(label=rownames(data)))
print(ggplot(tsne_data1,aes(x=Component_1,y=Component_2)) + geom_point(colour="light grey",size=1) +
      geom_point(data=tsne_data1[tsne_data1$CHEMOKINES > quantile(tsne_data1$MHC_II,0.9),],aes(x=Component_1, y=Component_2), colour="blue", size=1))

min_cyt = which.min(aggregate(TCELL_CYTOLYTIC_ACT ~ Cluster,tsne_data1,"mean")[,2])
max_cyt = which.max(aggregate(TCELL_CYTOLYTIC_ACT ~ Cluster,tsne_data1,"mean")[,2])
samples_up = rownames(tsne_data1[tsne_data1$Cluster %in% max_cyt,])
samples_down = rownames(tsne_data1[tsne_data1$Cluster %in% min_cyt,])
dat_new = cbind(dat_vst[,match(samples_up,colnames(dat_vst))],dat_vst[,match(samples_down,colnames(dat_vst))])
dim(dat_new)
grp1 = length(samples_up); grp2 = length(samples_down)
pheno_cyt = c(rep(0,grp1),rep(1,grp2))
pval_out = laply(1:nrow(dat_new),.progress="time",function(i) wilcox.test(as.numeric(dat_new[i,])~pheno_cyt)$p.val)
pval_out1 = data.frame("Genes"=rownames(dat_new),"pval"=pval_out,"pval_adj"=p.adjust(pval_out),"qval"=qvalue(pval_out)$qval)
pval_sig = pval_out1[pval_out1$qval<=0.01,]
pval_sig = pval_sig[order(pval_sig$pval),]
pval_sig = na.omit(pval_sig)
means = do.call("rbind",llply(1:nrow(pval_sig),function(i) aggregate(as.numeric(dat_new[grep(paste("^",pval_sig$Genes[i],"$",sep=""),rownames(dat_new)),]) ~ pheno_cyt,FUN="mean")[,-1],.progress="time"))
colnames(means) = c("Group1_means","Group2_means")
cyt_sig = cbind(pval_sig,means)
dim(cyt_sig)
cyt_sig$Genes = as.character(cyt_sig$Genes)
write.table(cyt_sig,"CYT_DE_sig_genes.txt",sep="\t",col.names=NA,quote=F)
dat_cyt = dat_new[match(cyt_sig$Genes,rownames(dat_new)),]

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
#
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

### 2) MD Validation
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

### Common probes between TCGA/METABRIC and single-cell RNAseq
classifier_genes = Reduce("intersect",list(rownames(dat_cyt),rownames(tcgaTN),rownames(metaB_TN),rownames(metaBV_TN),
                                          rownames(tcgaMELA),rownames(tcgaLUAD),rownames(tcgaLUSC),rownames(tcgaHNSC),rownames(tcgaCOAD),rownames(tcgaKIRC)))


### Training datasets
dat_cyt_train = dat_cyt[match(classifier_genes,rownames(dat_cyt)),]

### Test datasets
# BRCA
tcgaTN_cyt = tcgaTN[match(classifier_genes,rownames(tcgaTN)),]
tcgaER_cyt = tcgaER[match(classifier_genes,rownames(tcgaER)),]
metaB_TN_cyt = metaB_TN[match(classifier_genes,rownames(metaB_TN)),]
metaBV_TN_cyt = metaBV_TN[match(classifier_genes,rownames(metaBV_TN)),]

# TCGA datasets
tcgaMELA_cyt = tcgaMELA[match(classifier_genes,rownames(tcgaMELA)),]
tcgaLUAD_cyt = tcgaLUAD[match(classifier_genes,rownames(tcgaLUAD)),]
tcgaLUSC_cyt = tcgaLUSC[match(classifier_genes,rownames(tcgaLUSC)),]
tcgaHNSC_cyt = tcgaHNSC[match(classifier_genes,rownames(tcgaHNSC)),]
tcgaCOAD_cyt = tcgaCOAD[match(classifier_genes,rownames(tcgaCOAD)),]
tcgaKIRC_cyt = tcgaKIRC[match(classifier_genes,rownames(tcgaKIRC)),]

### Training datasets
seed = 12345
set.seed(seed)
tr_ctrl = trainControl(method="oob",classProbs=TRUE,verboseIter=FALSE,sampling = "down")
mtryGrid = expand.grid(mtry = c(50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,1000))
trainIndex = createDataPartition(pheno_cyt,0.40,list=F,times=1)
trainData = t(dat_cyt_train[,trainIndex])
trainPheno = as.numeric(as.character(pheno_cyt[as.numeric(trainIndex)]))
trainPheno[trainPheno=="0"]<-"UP"
trainPheno[trainPheno=="1"]<-"DOWN"
trainPheno = as.factor(trainPheno)
testData = t(dat_cyt_train[,-trainIndex])
testPheno = as.numeric(as.character(pheno_cyt[-as.numeric(trainIndex)]))
testPheno[testPheno=="0"]<-"UP"
testPheno[testPheno=="1"]<-"DOWN"
testPheno = as.factor(testPheno)

set.seed(seed);
train_out = train(trainData, trainPheno, method = "rf",importance=TRUE,metric="Accuracy",proximity=TRUE,trControl = tr_ctrl,tuneGrid = mtryGrid,ntree=100)

pdf("RF_valid_survCurves.pdf")
### Validate classifier in TCGA Breast Data
## Triple-negative
clinTN$Prob = predict(train_out,apply(tcgaTN_cyt,1,scale),type="prob")[,1]
clinTN$CYT_ACT = as.numeric(tcgaTN_gsva[4,])
clinTN$Preds = predict(train_out,apply(tcgaTN_cyt,1,scale),type="raw")
tcga_brca_pfs = as.numeric(pairwise_survdiff(Surv(pfs_year,pfs_Ind) ~ Preds,clinTN,p.adjust="none")$p.value)
tcga_brca_os = as.numeric(pairwise_survdiff(Surv(os_year,os_Ind) ~ Preds,clinTN,p.adjust="none")$p.value)

fitDSS_TCGA = survfit(Surv(pfs_year,pfs_Ind) ~ Preds,clinTN)
ggsurvplot(fitDSS_TCGA,palette = genespring.colors(2),risk.table=T,pval=T,data=clinTN,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")

fitOS_TCGA = survfit(Surv(os_year,os_Ind) ~ Preds,clinTN)
ggsurvplot(fitOS_TCGA,palette = genespring.colors(2),risk.table=T,pval=T,data=clinTN,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")

ggplotRegression(lm(as.numeric(tcgaTN_gsva[4,]) ~ clinTN$Prob))
summary(coxph(Surv(pfs_year,pfs_Ind) ~ CYT_ACT,clinTN))



## ER positive
clinER$Prob = predict(train_out,apply(tcgaER_cyt,1,scale),type="prob")[,1]
clinER$CYT_ACT = as.numeric(tcgaER_gsva[4,])
clinER$Preds = predict(train_out,apply(tcgaER_cyt,1,scale),type="raw")
tcga_er_brca_pfs = as.numeric(pairwise_survdiff(Surv(pfs_year,pfs_Ind) ~ Preds,clinER,p.adjust="none")$p.value)
tcga_er_brca_os = as.numeric(pairwise_survdiff(Surv(os_year,os_Ind) ~ Preds,clinER,p.adjust="none")$p.value)

fitDSS_ER_TCGA = survfit(Surv(pfs_year,pfs_Ind) ~ Preds,clinER)
ggsurvplot(fitDSS_ER_TCGA,palette = genespring.colors(2),risk.table=T,pval=T,data=clinER,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")
fitOS_ER_TCGA = survfit(Surv(os_year,os_Ind) ~ Preds,clinER)
ggsurvplot(fitOS_ER_TCGA,palette = genespring.colors(2),risk.table=T,pval=T,data=clinER,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")

ggplotRegression(lm(as.numeric(tcgaER_gsva[4,]) ~ clinER$Prob))

### Validate classifier in METABRIC Data
clinTN_metaB$Probs = predict(train_out,apply(metaB_TN_cyt,1,scale),type="prob")[,1]
clinTN_metaB$CYT_ACT = as.numeric(metaB_TN_gsva[4,])
clinTN_metaB$Preds = predict(train_out,apply(metaB_TN_cyt,1,scale))
clinTN_metaBV$Probs = predict(train_out,apply(metaBV_TN_cyt,1,scale),type="prob")[,1]
clinTN_metaBV$CYT_ACT = as.numeric(metaBV_TN_gsva[4,])
clinTN_metaBV$Preds = predict(train_out,apply(metaBV_TN_cyt,1,scale))
metaB_valid = rbind(clinTN_metaB,clinTN_metaBV)
metab_brca_pfs = as.numeric(pairwise_survdiff(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_valid,p.adjust="none")$p.value)
metab_brca_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_Ind) ~ Preds,metaB_valid,p.adjust="none")$p.value)

fitDSS_META = survfit(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_valid)
ggsurvplot(fitDSS_META,palette = genespring.colors(2),risk.table=T,pval=T,data=metaB_valid,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")
fitOS_META = survfit(Surv(OS_year,OS_Ind) ~ Preds,metaB_valid)
ggsurvplot(fitOS_META,palette = genespring.colors(2),risk.table=T,pval=T,data=metaB_valid,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")

ggplotRegression(lm(CYT_ACT ~ Probs,data=metaB_valid))

tcga_clinMELA$Probs = predict(train_out,apply(tcgaMELA_cyt,1,scale),type="prob")[,1]
tcga_clinMELA$CYT_ACT = as.numeric(tcga_mela_gsva[4,])
tcga_clinMELA$Preds = predict(train_out,apply(tcgaMELA_cyt,1,scale),type="raw")
tcga_skcm_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinMELA,p.adjust="none")$p.value)
tcga_skcm_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinMELA,p.adjust="none")$p.value)
fitDSS_SKCM = survfit(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinMELA)
ggsurvplot(fitDSS_SKCM,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinMELA,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")
fitOS_SKCM = survfit(Surv(OS_year,OS_ind) ~ Preds,tcga_clinMELA)
ggsurvplot(fitOS_SKCM,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinMELA,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")

ggplotRegression(lm(CYT_ACT ~ Probs,data=tcga_clinMELA))

# print("TCGA_LUAD")
tcga_clinLUAD$Probs = predict(train_out,apply(tcgaLUAD_cyt,1,scale),type="prob")[,1]
tcga_clinLUAD$CYT_ACT = as.numeric(tcga_luad_gsva[4,])
tcga_clinLUAD$Preds = predict(train_out,apply(tcgaLUAD_cyt,1,scale),type="raw")
tcga_luad_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinLUAD,p.adjust="none")$p.value)
tcga_luad_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinLUAD,p.adjust="none")$p.value)
fitDSS_LUAD = survfit(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinLUAD)
ggsurvplot(fitDSS_LUAD,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinLUAD,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")
fitOS_LUAD = survfit(Surv(OS_year,OS_ind) ~ Preds,tcga_clinMELA)
ggsurvplot(fitOS_LUAD,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinLUAD,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")

ggplotRegression(lm(CYT_ACT ~ Probs,data=tcga_clinLUAD))

# print("TCGA_LUSC")
tcga_clinLUSC$Probs = predict(train_out,apply(tcgaLUSC_cyt,1,scale),type="prob")[,1]
tcga_clinLUSC$CYT_ACT = as.numeric(tcga_lusc_gsva[4,])
tcga_clinLUSC$Preds = predict(train_out,apply(tcgaLUSC_cyt,1,scale),type="raw")
tcga_lusc_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinLUSC,p.adjust="none")$p.value)
tcga_lusc_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinLUSC,p.adjust="none")$p.value)
fitDSS_LUSC = survfit(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinLUSC)
ggsurvplot(fitDSS_LUSC,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinLUSC,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")
fitOS_LUSC = survfit(Surv(OS_year,OS_ind) ~ Preds,tcga_clinLUSC)
ggsurvplot(fitOS_LUSC,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinLUSC,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")

ggplotRegression(lm(CYT_ACT ~ Probs,data=tcga_clinLUSC))

# print("TCGA_HNSC")
tcga_clinHNSC$Probs = predict(train_out,apply(tcgaHNSC_cyt,1,scale),type="prob")[,1]
tcga_clinHNSC$CYT_ACT = as.numeric(tcga_hnsc_gsva[4,])
tcga_clinHNSC$Preds = predict(train_out,apply(tcgaHNSC_cyt,1,scale),type="raw")
tcga_hnsc_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinHNSC,p.adjust="none")$p.value)
tcga_hnsc_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinHNSC,p.adjust="none")$p.value)
fitDSS_HNSC = survfit(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinHNSC)
ggsurvplot(fitDSS_HNSC,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinHNSC,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")
fitOS_HNSC = survfit(Surv(OS_year,OS_ind) ~ Preds,tcga_clinHNSC)
ggsurvplot(fitOS_HNSC,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinHNSC,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")

ggplotRegression(lm(CYT_ACT ~ Probs,data=tcga_clinHNSC))

# print("TCGA_COAD")
tcga_clinCOAD$Probs = predict(train_out,apply(tcgaCOAD_cyt,1,scale),type="prob")[,1]
tcga_clinCOAD$CYT_ACT = as.numeric(tcga_coad_gsva[4,])
tcga_clinCOAD$Preds = predict(train_out,apply(tcgaCOAD_cyt,1,scale),type="raw")
tcga_coad_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinCOAD,p.adjust="none")$p.value)
tcga_coad_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinCOAD,p.adjust="none")$p.value)
fitDSS_COAD = survfit(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinCOAD)
ggsurvplot(fitDSS_COAD,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinCOAD,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")
fitOS_COAD = survfit(Surv(OS_year,OS_ind) ~ Preds,tcga_clinCOAD)
ggsurvplot(fitOS_COAD,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinCOAD,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")

ggplotRegression(lm(CYT_ACT ~ Probs,data=tcga_clinCOAD))

# print("TCGA_KIRC")
tcga_clinKIRC$Probs = predict(train_out,apply(tcgaKIRC_cyt,1,scale),type="prob")[,1]
tcga_clinKIRC$CYT_ACT = as.numeric(tcga_kirc_gsva[4,])
tcga_clinKIRC$Preds = predict(train_out,apply(tcgaKIRC_cyt,1,scale),type="raw")
tcga_kirc_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinKIRC,p.adjust="none")$p.value)
tcga_kirc_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinKIRC,p.adjust="none")$p.value)
fitDSS_KIRC = survfit(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinKIRC)
ggsurvplot(fitDSS_KIRC,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinKIRC,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")
fitOS_KIRC = survfit(Surv(OS_year,OS_ind) ~ Preds,tcga_clinKIRC)
ggsurvplot(fitOS_KIRC,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinKIRC,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")

ggplotRegression(lm(CYT_ACT ~ Probs,data=tcga_clinKIRC))

dev.off()
out = data.frame("TCGA_BR_PFS"=tcga_brca_pfs,"TCGA_BR_OS"=tcga_brca_os,
                "META_BR_PFS"= metab_brca_pfs,"META_BR_OS"=metab_brca_os,
                "TCGA_SKCM_PFS"= tcga_skcm_pfs,"TCGA_SKCM_OS"=tcga_skcm_os,
                "TCGA_LUAD_PFS"= tcga_luad_pfs,"TCGA_LUAD_OS"=tcga_luad_os,
                "TCGA_LUSC_PFS"= tcga_lusc_pfs,"TCGA_LUSC_OS"=tcga_lusc_os,
                "TCGA_HNSC_PFS"= tcga_hnsc_pfs,"TCGA_HNSC_OS"=tcga_hnsc_os,
                "TCGA_COAD_PFS"= tcga_coad_pfs,"TCGA_COAD_OS"=tcga_coad_os,
                "TCGA_KIRC_PFS"= tcga_kirc_pfs,"TCGA_KIRC_OS"=tcga_kirc_os)

fin_mod = train_out$finalModel
pdf("rf_final_model_plots.pdf",width=11.1,height=6.64)
MDSplot(fin_mod, fac=trainPheno, k=2, palette=genespring.colors(2), pch=19)
abline(h=0,v=0,lty=2)
varImpPlot(fin_mod,n.var=25,pch=19,main="Random Forest -- Final Model")
plot(fin_mod,log="y",col=c("black","red","blue"),main="Final Model")
legend("topright", colnames(fin_mod$err.rate),col=c("black","red","blue"),cex=0.8,fill=c("black","red","blue"))
dev.off()

test_prob=data.frame(predict(train_out,testData,type="prob"))
test_prob$true_class = testPheno
test_prob$pred_class = predict(train_out,testData)
test_prob = test_prob[order(test_prob$true_class),]
chisq.test(table(test_prob$true_class,test_prob$pred_class))
pdf("test_set_validation.pdf",width=11.1,height=6.64)
plot(test_prob$DOWN,col=genespring.colors(2)[as.numeric(test_prob$true_class)],pch=19,ylab="Random Forest Probability",xlab="Patients in Training Set")
abline(h=0.5,lty=2)
dev.off()
predictors = data.frame(varImp(train_out,scale=T)$importance)
predictors$Genes = rownames(predictors)
predictors = predictors[order(predictors$DOWN,decreasing=T),]
ggbarplot(predictors[150:1,],x="Genes",y="DOWN",fill="grey",ylab="Scaled Variable Importance",lab.size=1)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=rel(0.5)))+theme(text = element_text(size=5))+coord_flip()
ggsave("top_predictors.pdf",width=11.1,height=6.64,units="in")

cyt_sig_classifier = cyt_sig[match(classifier_genes,cyt_sig$Genes),]
pdf("volcano_plot.pdf")
with(cyt_sig_classifier,plot(logFC,-log10(qval),pch=20,main="",xlim=c(-5,5),xlab="LOG FOLD CHANGE",ylab="-log10 ADJ P VAL"))
with(subset(cyt_sig_classifier, abs(logFC)>2), textxy(logFC, -log10(qval), labs=Genes, cex=.5))
dev.off()

########
#
