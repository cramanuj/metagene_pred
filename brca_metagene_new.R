source("~/Research/scripts/r_scripts/plotfns.R")
source("~/Research/scripts/r_scripts/useful_functions.R")

library(data.table)
# dr_out = fread("/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/BRCA_scRNAseq_gene_level.txt",header=T,data.table=F)
# rownames(dr_out) = dr_out[,1]
# dr_out = dr_out[,-1]
# dim(dr_out)
library(cellrangerRkit)
genome <- "GRCh38"
brca = load_cellranger_matrix("~/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final",genome=genome)
dim(exprs(brca))

brca_dat = data.frame(as.matrix(exprs(brca)),check.names=F)

probes = fread("~/Research/pathways/Homo_sapiens.GRCh38.79_table.txt",header=F,data.table=F)
colnames(probes)=c("ENSEMBL","Genes","Coordinates","Strand")
probes$Chr = gsub('\\:.*', '', probes$Coordinates)
probes = probes[c(which(probes$Chr %in% c(1:22)),grep("^X$",probes$Chr),grep("^Y$",probes$Chr)),]
probes = probes[-grep("^MIR",probes$Genes),]
probes = probes[-grep("^RPS",probes$Genes),]
probes = probes[-grep("^RPL",probes$Genes),]
keep_probes = intersect(rownames(brca_dat),probes$ENSEMBL)
seu_exprs = brca_dat[match(keep_probes,rownames(brca_dat)),]
probes_exprs = probes[match(keep_probes,probes$ENSEMBL),]

seu_exprs$Genes = as.factor(probes_exprs$Genes)
dat = setDT(seu_exprs)[, lapply(.SD, median), by = Genes]
dat = as.data.frame(dat)
rownames(dat) = dat[,1]
dat = dat[,-1]
fwrite(dat,"/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/BRCA_scRNAseq_gene_level.txt",sep="\t",col.names=T,quote=F,row.names=T)
pdata = data.frame(barcode=colnames(dat))
rownames(pdata) = pdata$barcode

##########
# MONOCLE
###############

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
pData(mycds)$Total_mRNAs = Matrix::colSums(exprs(mycds))
mycds = mycds[,pData(mycds)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(mycds)$Total_mRNAs)) + 2*sd(log10(pData(mycds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(mycds)$Total_mRNAs)) - 2*sd(log10(pData(mycds)$Total_mRNAs)))
qplot(Total_mRNAs, data = pData(mycds), color = "red", geom ="density") + geom_vline(xintercept = lower_bound) + geom_vline(xintercept = upper_bound)
mycds <- mycds[,pData(mycds)$Total_mRNAs > lower_bound &   pData(mycds)$Total_mRNAs < upper_bound]
L <- log(exprs(mycds[expressed_genes,]))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab("Standardized log(FPKM)") + ylab("Density")
# x_1 = (x - mean(x)) / sd(x)
# df <- data.frame(x = x_1)
# ggplot(df, aes(x)) + geom_histogram(bins = 50) + geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
# disp_table <- dispersionTable(mycds)
# pData(mycds)$UMI <- Matrix::colSums(exprs(mycds))
# ggplot(pData(mycds), aes(num_genes_expressed, UMI)) + geom_point()
# unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
# mycds <- setOrderingFilter(mycds, unsup_clustering_genes$gene_id)
# plot_ordering_genes(mycds)
# plot_pc_variance_explained(mycds, return_all = FALSE)
# mycds = reduceDimension(mycds, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = T)
# mycds = clusterCells(mycds, num_clusters = 10)
# plot_cell_clusters(mycds, 1, 2, color = "Cluster")
# mycds = reduceDimension(mycds, max_components = 2, method = 'DDRTree')
# mycds = orderCells(mycds)
# plot_cell_trajectory(mycds, color_by = "Cluster")

mon_exprs = t(t(exprs(mycds)) /  pData(mycds)[, 'Size_Factor'])
dat_vst = vstExprs(mycds,expr_matrix=mon_exprs)
fwrite(data.frame(dat_vst,check.names=F),"BRCA_scRNAseq_VST.txt",sep="\t",col.names=T,row.names=T,quote=F)

# cth <- newCellTypeHierarchy()
# cth <- addCellType(cth, "CD8+ T cell",classify_func=function(x) { x["CD8A",] >= 1 | x["CD8B",] >= 1 })
# cth <- addCellType(cth, "CD4+ T cell", classify_func=function(x) {x["CD4",] >= 1})
# mycds <- classifyCells(mycds, cth)
# table(pData(mycds)$CellType)

# x = pData(mycds)$num_genes_expressed
# x_1 = (x - mean(x)) / sd(x)
# df <- data.frame(x = x_1)
# ggplot(df, aes(x)) + geom_histogram(bins = 50) + geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
# disp_table <- dispersionTable(mycds)
# pData(mycds)$UMI_monocle <- Matrix::colSums(exprs(mycds))
# ggplot(pData(mycds), aes(num_genes_expressed, UMI_monocle)) + geom_point()

plot_pc_variance_explained(mycds, return_all = FALSE)
mycds = reduceDimension(mycds, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = T)
mycds = clusterCells(mycds, num_clusters = 10)
plot_cell_clusters(mycds, 1, 2, color = "Cluster")
plot_cell_clusters(mycds, 1, 2, color = "Cluster",markers = c("CD8A","CD8B","CD4"))

# tsne_data = data.frame(t(reducedDimA(mycds)))
# colnames(tsne_data)=c("Component_1","Component_2")
# rownames(tsne_data) = colnames(exprs(mycds))
# datN = log10(Matrix::t(Matrix::t(exprs(mycds))/sizeFactors(mycds))+0.1)
# tsne_data$CD8A = as.numeric(datN[grep("CD8A",rownames(datN)),])
# tsne_data$CD8B = as.numeric(datN[grep("CD8B",rownames(datN)),])
# tsne_data$CD4 = as.numeric(datN[grep("^CD4$",rownames(datN)),])
#
# ggplot(tsne_data,aes(x=Component_1,y=Component_2)) + geom_point(colour="light grey",size=1) +
#       geom_point(data=tsne_data[tsne_data$CD8A > quantile(tsne_data$CD8A,0.75),],aes(x=Component_1, y=Component_2), colour="blue", size=1) +
#       geom_point(data=tsne_data[tsne_data$CD4 > quantile(tsne_data$CD4,0.75),],aes(x=Component_1, y=Component_2), colour="red", size=1)
#
# ### CD8+ subset
# CD8_subset = row.names(subset(pData(mycds), CellType == "CD8+ T cell"))
# mycds_cd8 = mycds[,CD8_subset]
# mycds_cd8 = reduceDimension(mycds_cd8, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = T)
# mycds_cd8 = clusterCells(mycds_cd8, num_clusters = 5)
# plot_cell_clusters(mycds_cd8, 1, 2, color = "Cluster")
# plot_cell_clusters(mycds_cd8, 1, 2, color = "Cluster",markers= c("CD8A","CD8B","CD4"))
# plot_cell_clusters(mycds_cd8, 1, 2, color = "Cluster",markers= c("PRF1","GZMA"))
# mon_exprs = t(t(exprs(mycds_cd8)) /  pData(mycds_cd8)[, 'Size_Factor'])
# datVST_CD8 = vstExprs(mycds_cd8,expr_matrix=mon_exprs)

########
# ssGSEA
################

library(GSVA)
imm = read.gmt.file("~/Research/pathways/immune_system.gmt")
ssgsea_out = gsva(as.matrix(dat_vst),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Poisson",min.sz=1,max.sz=1000)
write.table(ssgsea_out,"BRCA_scRNAseq_ssGSEA.txt",sep="\t",col.names=NA,quote=F)
# ssgsea_cd8_out = gsva(as.matrix(datVST_CD8),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Poisson",min.sz=1,max.sz=1000)

#### ssGSEA global
## Dimension reduction + Density clustering
############################################
library(Rtsne); library(densityClust); library(ggthemes); library(plyr); library(dplyr); library(ggcorrplot); library(qvalue); library(ggpubr)
FM = ssgsea_out
max_components = 2; num_dim = 3; num_clusters = 10; perp = 30
set.seed(123)
irlba_res <- prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1), center = TRUE, scale. = TRUE)
irlba_pca_res <- irlba_res$x
topDim_pca <- irlba_pca_res
tsne_res <- Rtsne(as.matrix(topDim_pca), dims = max_components,pca = F,perplexity=perp)
tsne_data <- tsne_res$Y[, 1:max_components]
colnames(tsne_data) = c("Component_1", "Component_2")
rownames(tsne_data) = colnames(FM)
dataDist <- dist(tsne_data)
dataClust <- densityClust(dataDist, gaussian = T)
delta_rho_df <- data.frame(delta = dataClust$delta, rho = dataClust$rho)
rho_threshold <- 0
delta_threshold <- sort(delta_rho_df$delta, decreasing = T)[num_clusters] - .Machine$double.eps
dataClust <- densityClust::findClusters(dataClust, rho = rho_threshold, delta = delta_threshold)
tsne_data = data.frame(tsne_data)
tsne_data$Cluster = as.factor(dataClust$cluster)
tcenters = tsne_data[,c(1:3)] %>% dplyr::group_by(Cluster) %>% summarize_all(funs(median))
plot1 = ggplot(data.frame(tsne_data),aes(x=Component_1,y=Component_2,colour=Cluster,label=Cluster)) + geom_point(size=1) + geom_point(data = tcenters, mapping = aes(x = Component_1, y = Component_2), size = 0, alpha = 0) +
      geom_text(data=tcenters,mapping = aes(label = Cluster), colour="black",size = 6) + theme_pander(boxes=T) + theme(legend.position="none")
tsne_data1 = cbind(tsne_data,apply(FM,1,scale))
ave_out = aggregate(. ~ Cluster,tsne_data1,"mean")[,-c(1:3)]
corr_out = cor(ave_out)
ggcorrplot(corr_out, hc.order = TRUE,type = "lower", lab = TRUE, lab_size = 2, method="circle", colors = c("tomato2", "white", "springgreen3"),
            title="Correlogram of immune metagenes", tl.cex = 6, ggtheme=theme_bw)
ggsave("corr_plot_ssGSEA.pdf",width=10.9,height=6.71)
centers = tsne_data1[,c(1:3)] %>% dplyr::group_by(Cluster) %>% summarize_all(funs(median))
plot2 = ggplot(tsne_data1,aes(x=Component_1,y=Component_2)) + geom_point(colour="light grey",size=1) +
      geom_point(data=tsne_data1[tsne_data1$TCELL_CYTOLYTIC_ACT > quantile(tsne_data1$TCELL_CYTOLYTIC_ACT,0.9),],aes(x=Component_1, y=Component_2), colour="red", size=1) +
      geom_point(data = centers, mapping = aes(x = Component_1, y = Component_2), size = 0, alpha = 0) + geom_text(data=centers,mapping = aes(label = Cluster), colour="black",size = 5) +
      geom_point(data=tsne_data1[tsne_data1$TCELL_CYTOLYTIC_ACT < quantile(tsne_data1$TCELL_CYTOLYTIC_ACT,0.1),],aes(x=Component_1, y=Component_2), colour="blue", size=1) +
      geom_point(data = centers, mapping = aes(x = Component_1, y = Component_2), size = 0, alpha = 0) + geom_text(data=centers,mapping = aes(label = Cluster), colour="black",size = 5) + theme_pander(boxes=T)
ggarrange(plot1, plot2, nrow=1,ncol = 2)
ggsave("tSNE_plots_ssGSEA.pdf",width=10.9,height=6.71)

pdf("temp_boxplots_ssgsea.pdf")
for(i in 1:nrow(FM)){
  print(i)
  print(ggboxplot(tsne_data1,x="Cluster",y=rownames(FM)[i],fill="grey"))
}
dev.off()

min_cyt = which.min(aggregate(TCELL_CYTOLYTIC_ACT ~ Cluster,tsne_data1,"median")[,2])
max_cyt = which.max(aggregate(TCELL_CYTOLYTIC_ACT ~ Cluster,tsne_data1,"median")[,2])
samples_up = rownames(tsne_data1[tsne_data1$Cluster %in% max_cyt,])
samples_down = rownames(tsne_data1[tsne_data1$Cluster %in% min_cyt,])

# samples_down = rownames(tsne_data1)[which(tsne_data1$TCELL_CYTOLYTIC_ACT <= quantile(tsne_data1$TCELL_CYTOLYTIC_ACT,0.10))]
# samples_up = rownames(tsne_data1)[which(tsne_data1$TCELL_CYTOLYTIC_ACT > quantile(tsne_data1$TCELL_CYTOLYTIC_ACT,0.90))]
# dat_new = cbind(dat_vst[,match(samples_up,colnames(dat_vst))],dat_vst[,match(samples_down,colnames(dat_vst))])
# dat_new = cbind(datVST_CD8[,match(samples_up,colnames(datVST_CD8))],datVST_CD8[,match(samples_down,colnames(datVST_CD8))])

grp1 = length(samples_up); grp2 = length(samples_down)
samples = c(samples_up,samples_down)
pheno_cyt = c(rep(0,grp1),rep(1,grp2))
pheno_cyt[pheno_cyt>0]<- "TCA_DOWN"
pheno_cyt[pheno_cyt!="TCA_DOWN"]<- "TCA_UP"
pheno_cyt = as.factor(pheno_cyt)
pheno_cyt = relevel(pheno_cyt,ref="TCA_DOWN")

# pval_out = laply(1:nrow(dat_new),.progress="time",function(i) wilcox.test(as.numeric(dat_new[i,])~pheno_cyt)$p.val)
# pval_out = data.frame("Genes"=rownames(dat_new),"pval"=pval_out,"pval_adj"=p.adjust(pval_out),"qval"=qvalue(pval_out)$qval)
# pval_out = pval_out[order(pval_out$pval),]
# pval_out = na.omit(pval_out)
# pval_sig = pval_out[pval_out$qval<=0.05,]
# dat_cyt = dat_new[match(pval_sig$Genes,rownames(dat_new)),]

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

# out = FindMarkers(seu,ident.1=0,ident.2=2,only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
seu_new = SubsetData(seu,idents.use=c(0,2))
seu_new@data = data.frame(as.matrix(seu_new@data),check.names=F)
seu_new = NormalizeData(seu_new, normalization.method = "LogNormalize", scale.factor = 10000)
seu_new = FindVariableGenes(seu_new, do.plot = T, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.015, x.high.cutoff = 3.2, y.cutoff = 0.5)
seu_new = RunPCA(seu_new, pcs.print = 0,pc.genes = seu@var.genes)
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2)
PCHeatmap(object = seu_new, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE)

seu_new = FindClusters(object = seu_new, reduction.type = "pca", dims.use = 1:5, resolution = 1, print.output = 0, save.SNN = TRUE)
seu_new = RunTSNE(seu_new, dims.use = 1:5,do.fast=T)
TSNEPlot(seu_new, do.label = TRUE, pt.size = 0.5)
FeaturePlot(seu_new, features.plot = c("PRF1","GZMA"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
markers <- FindAllMarkers(object = seu_new,only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
top10 = markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = seu_new, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,cex.row = 5,group.cex = 5)
seu_new = JackStraw(object = seu_new, num.replicate = 100)
JackStrawPlot(object = seu_new, PCs = 1:20)

dat_new = data.frame(as.matrix(seu_new@data),check.names=F)
pc1_genes = names(which(seu_new@dr$pca@jackstraw@emperical.p.value[,1]<=0.05))

dat_cyt = dat_new[match(pc1_genes,rownames(dat_new)),]
pheno_cyt = seu_new@meta.data$pheno
pheno_cyt = relevel(pheno_cyt,ref="TCA_DOWN")
pc1_genes = names(which(seu@dr$pca@jackstraw@emperical.p.value[,1]<=0.05))
# pc2_genes = names(which(seu@dr$pca@jackstraw@emperical.p.value[,2]<=0.05))
# pc3_genes = names(which(seu@dr$pca@jackstraw@emperical.p.value[,3]<=0.05))
# pc4_genes = names(which(seu@dr$pca@jackstraw@emperical.p.value[,4]<=0.05))
# pc5_genes = names(which(seu@dr$pca@jackstraw@emperical.p.value[,5]<=0.05))

## Gene score
# library(genefu)
# annot = data.frame("EntrezGene.ID"=rownames(mark1))
# x = data.frame("probe"=rownames(mark1),"EntrezGene.ID"=rownames(mark1),"coefficient"=mark1$avg_logFC)
# gene_score = sig.score(x,data=t(dat_cyt),annot)

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

##################
### Common probes between TCGA/METABRIC and single-cell RNAseq
################################################################

classifier_genes = Reduce("intersect",list(pc1_genes,rownames(tcgaTN),rownames(metaB_TN),rownames(metaBV_TN),
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

# annot = data.frame("EntrezGene.ID"=rownames(mark2))
# x = data.frame("probe"=rownames(mark2),"EntrezGene.ID"=rownames(mark2),"coefficient"=mark2$avg_logFC)
# tcga_brca_score = sig.score(x,data=t(tcgaTN_cyt),annot)
# mb2_brca_score = sig.score(x,data=t(metaBV_TN_cyt),annot)
# temp = data.frame("Samples"=clinTN_metaBV$PATIENT_ID,"Score"=mb2_brca_score$score)
# head(temp)
# temp$groups = ifelse(temp$Score<quantile(mb2_brca_score$score,0.25),"TCR_DOWN","TCR_UP")
# head(temp)
# temp$year = clinTN_metaBV$DSS_days
# temp$Ind = clinTN_metaBV$DSS_event
# fitDSS_MB = survfit(Surv(year,Ind) ~ groups,temp)
# ggsurvplot(fitDSS_MB,palette = genespring.colors(2),risk.table=T,pval=T,data=temp,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_pander(),legend.title="")
