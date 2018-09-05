load("metagene_predictor_newer.rds")
library(caret); library(ggplot2); library(doParallel); library(parallel)
library(survminer); library(survival)

cluster = makeCluster(10)
registerDoParallel(cluster)

pc1_genes = names(which(seu@dr$pca@jackstraw@emperical.p.value[,1]<=0.05))
pc2_genes = names(which(seu@dr$pca@jackstraw@emperical.p.value[,2]<=0.05))
pc3_genes = names(which(seu@dr$pca@jackstraw@emperical.p.value[,3]<=0.05))
pc4_genes = names(which(seu@dr$pca@jackstraw@emperical.p.value[,4]<=0.05))
pc5_genes = names(which(seu@dr$pca@jackstraw@emperical.p.value[,5]<=0.05))

pc_gene_list = list(pc1_genes,pc2_genes,pc3_genes,pc4_genes,pc5_genes)
pc_out = list()

for(k in 1:length(pc_gene_list)){
	cat("PC Metagene: ",k,"of",length(pc_gene_list),"\n")
	pc_genes = pc_gene_list[[k]]
	dat_cyt = dat_new[match(pc_genes,rownames(dat_new)),]
	pheno_cyt = seu@meta.data$pheno
	pheno_cyt[pheno_cyt==0]<-"TCA_UP"; pheno_cyt[pheno_cyt==1]<-"TCA_DOWN";
	pheno_cyt = as.factor(pheno_cyt);
	pheno_cyt = relevel(pheno_cyt,ref="TCA_DOWN")
	classifier_genes = Reduce("intersect",list(pc_genes,rownames(tcgaTN),rownames(metaB_TN),rownames(metaBV_TN),
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

	seed = 12345
	tr_ctrl = trainControl(method="oob",classProbs=TRUE,verboseIter=FALSE,sampling = "down", allowParallel=TRUE)
	mtryGrid = expand.grid(mtry = seq(50,length(classifier_genes),50))
	ntrees = c(10,15,25,50,75,100,175,190,200,225,250,275,300,325,350,375,400,450,500)
	prop = c(0.25,0.30,0.35,0.4,0.45,0.5,0.6,0.75,1)
	fin=list()
	for(j in 1:length(prop)){
		cat("Train Prop: ",prop[j]*100,"% \n")
		set.seed(seed)
		trainIndex = createDataPartition(pheno_cyt,prop[j],list=F,times=1)
		trainData = t(dat_cyt_train[,trainIndex])
		trainPheno = pheno_cyt[as.numeric(trainIndex)]
		# testData = t(dat_cyt_train[,-trainIndex])
		# testPheno = pheno_cyt[-as.numeric(trainIndex)]
		out = data.frame(matrix(0,ncol=16,nrow=length(ntrees)))
		for(i in 1:length(ntrees)){
			cat("ntrees: ",ntrees[i],"\n")
			set.seed(seed);
			train_out = train(trainData, trainPheno, method = "rf",importance=TRUE,metric="Accuracy",proximity=TRUE,trControl = tr_ctrl,tuneGrid = mtryGrid,ntree=ntrees[i])

			print("TCGA_BRCA")
			clinTN$Prob = predict(train_out,apply(tcgaTN_cyt,1,scale),type="prob")[,1]
			clinTN$Preds = predict(train_out,apply(tcgaTN_cyt,1,scale),type="raw")
			tcga_brca_pfs = as.numeric(pairwise_survdiff(Surv(pfs_year,pfs_Ind) ~ Preds,clinTN,p.adjust="none")$p.value)
			tcga_brca_pfs = ifelse(length(tcga_brca_pfs)<1,NA,tcga_brca_pfs)
			tcga_brca_os = as.numeric(pairwise_survdiff(Surv(os_year,os_Ind) ~ Preds,clinTN,p.adjust="none")$p.value)
			tcga_brca_os = ifelse(length(tcga_brca_os)<1,NA,tcga_brca_os)

			print("METAB_BRCA")
			clinTN_metaB$Probs = predict(train_out,apply(metaB_TN_cyt,1,scale),type="prob")[,1]
			clinTN_metaB$Preds = predict(train_out,apply(metaB_TN_cyt,1,scale))
			clinTN_metaBV$Probs = predict(train_out,apply(metaBV_TN_cyt,1,scale),type="prob")[,1]
			clinTN_metaBV$Preds = predict(train_out,apply(metaBV_TN_cyt,1,scale))
			metaB_valid = rbind(clinTN_metaB,clinTN_metaBV)
			metab_brca_pfs = as.numeric(pairwise_survdiff(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_valid,p.adjust="none")$p.value)
			metab_brca_pfs = ifelse(length(metab_brca_pfs)<1,NA,metab_brca_pfs)
			metab_brca_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_Ind) ~ Preds,metaB_valid,p.adjust="none")$p.value)
			metab_brca_os = ifelse(length(metab_brca_os)<1,NA,metab_brca_os)

			print("TCGA_SKCM")
			tcga_clinMELA$Probs = predict(train_out,apply(tcgaMELA_cyt,1,scale),type="prob")[,1]
			tcga_clinMELA$Preds = predict(train_out,apply(tcgaMELA_cyt,1,scale),type="raw")
			tcga_skcm_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinMELA,p.adjust="none")$p.value)
			tcga_skcm_pfs = ifelse(length(tcga_skcm_pfs)<1,NA,tcga_skcm_pfs)
			tcga_skcm_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinMELA,p.adjust="none")$p.value)
			tcga_skcm_os = ifelse(length(tcga_skcm_os)<1,NA,tcga_skcm_os)

			print("TCGA_LUAD")
			tcga_clinLUAD$Probs = predict(train_out,apply(tcgaLUAD_cyt,1,scale),type="prob")[,1]
			tcga_clinLUAD$Preds = predict(train_out,apply(tcgaLUAD_cyt,1,scale),type="raw")
			tcga_luad_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinLUAD,p.adjust="none")$p.value)
			tcga_luad_pfs = ifelse(length(tcga_luad_pfs)<1,NA,tcga_luad_pfs)
			tcga_luad_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinLUAD,p.adjust="none")$p.value)
			tcga_luad_os = ifelse(length(tcga_luad_os)<1,NA,tcga_luad_os)

			print("TCGA_LUSC")
			tcga_clinLUSC$Probs = predict(train_out,apply(tcgaLUSC_cyt,1,scale),type="prob")[,1]
			tcga_clinLUSC$Preds = predict(train_out,apply(tcgaLUSC_cyt,1,scale),type="raw")
			tcga_lusc_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinLUSC,p.adjust="none")$p.value)
			tcga_lusc_pfs = ifelse(length(tcga_lusc_pfs)<1,NA,tcga_lusc_pfs)
			tcga_lusc_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinLUSC,p.adjust="none")$p.value)
			tcga_lusc_os = ifelse(length(tcga_lusc_os)<1,NA,tcga_lusc_os)

			print("TCGA_HNSC")
			tcga_clinHNSC$Probs = predict(train_out,apply(tcgaHNSC_cyt,1,scale),type="prob")[,1]
			tcga_clinHNSC$Preds = predict(train_out,apply(tcgaHNSC_cyt,1,scale),type="raw")
			tcga_hnsc_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinHNSC,p.adjust="none")$p.value)
			tcga_hnsc_pfs = ifelse(length(tcga_hnsc_pfs)<1,NA,tcga_hnsc_pfs)
			tcga_hnsc_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinHNSC,p.adjust="none")$p.value)
			tcga_hnsc_os = ifelse(length(tcga_hnsc_os)<1,NA,tcga_hnsc_os)

			print("TCGA_COAD")
			tcga_clinCOAD$Probs = predict(train_out,apply(tcgaCOAD_cyt,1,scale),type="prob")[,1]
			tcga_clinCOAD$Preds = predict(train_out,apply(tcgaCOAD_cyt,1,scale),type="raw")
			tcga_coad_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinCOAD,p.adjust="none")$p.value)
			tcga_coad_pfs = ifelse(length(tcga_coad_pfs)<1,NA,tcga_coad_pfs)
			tcga_coad_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinCOAD,p.adjust="none")$p.value)
			tcga_coad_os = ifelse(length(tcga_coad_os)<1,NA,tcga_coad_os)

			print("TCGA_KIRC")
			tcga_clinKIRC$Probs = predict(train_out,apply(tcgaKIRC_cyt,1,scale),type="prob")[,1]
			tcga_clinKIRC$Preds = predict(train_out,apply(tcgaKIRC_cyt,1,scale),type="raw")
			tcga_kirc_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinKIRC,p.adjust="none")$p.value)
			tcga_kirc_pfs = ifelse(length(tcga_kirc_pfs)<1,NA,tcga_kirc_pfs)
			tcga_kirc_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_ind) ~ Preds,tcga_clinKIRC,p.adjust="none")$p.value)
			tcga_kirc_os = ifelse(length(tcga_kirc_os)<1,NA,tcga_kirc_os)

			out[i,] = data.frame("TCGA_BR_PFS"=tcga_brca_pfs,"TCGA_BR_OS"=tcga_brca_os,
											"META_BR_PFS"=metab_brca_pfs,"META_BR_OS"=metab_brca_os,
											"TCGA_SKCM_PFS"= tcga_skcm_pfs,"TCGA_SKCM_OS"=tcga_skcm_os,
											"TCGA_LUAD_PFS"= tcga_luad_pfs,"TCGA_LUAD_OS"=tcga_luad_os,
											"TCGA_LUSC_PFS"= tcga_lusc_pfs,"TCGA_LUSC_OS"=tcga_lusc_os,
											"TCGA_HNSC_PFS"= tcga_hnsc_pfs,"TCGA_HNSC_OS"=tcga_hnsc_os,
											"TCGA_COAD_PFS"= tcga_coad_pfs,"TCGA_COAD_OS"=tcga_coad_os,
											"TCGA_KIRC_PFS"= tcga_kirc_pfs,"TCGA_KIRC_OS"=tcga_kirc_os)
		}
		rownames(out) = ntrees
		colnames(out) = c("TCGA_BR_PFS","TCGA_BR_OS","META_BR_PFS","META_BR_OS",
											"TCGA_SKCM_PFS","TCGA_SKCM_OS","TCGA_LUAD_PFS","TCGA_LUAD_OS",
											"TCGA_LUSC_PFS","TCGA_LUSC_OS","TCGA_HNSC_PFS","TCGA_HNSC_OS","TCGA_COAD_PFS","TCGA_COAD_OS","TCGA_KIRC_PFS","TCGA_KIRC_OS")

		fin[[j]] = out
	}
	fin_out = do.call("rbind",fin)
	fin_out$train_prop = rep(prop*100,each=length(ntrees))
	fin_out$ntrees = rep(ntrees,length(prop))
	fin_out$PC = k
	rownames(fin_out) = NULL
	pc_out[[k]] = fin_out
}
pc_fin = do.call("rbind",pc_out)
dim(pc_fin)
