# clean working environment

rm(list = ls())

# load necessary libraries

library(variancePartition)
library(car)
library(BiocParallel)
library(NMF)
library(ggplot2)
# 
# library(limma)
# library(edgeR)
# library(ggplot2)
# library(igraph)
# library(plyr)
# library(corrplot)
# library(Hmisc)

# Specify parallel processing parameters
# this is used by dream() to run in parallel

param = SnowParam(5, "SOCK", progressbar=TRUE)

# load preprocessed gene-count data

load("sources_of_variation_in_NPC_RNAseq.RData")

# change the values of the categorial variables that will go into the model into something more practical

phenotypes$Group <- NA
phenotypes$Group[phenotypes$Case_control == "Schizophrenia"] <- "SCZ"
phenotypes$Group[phenotypes$Case_control == "Control"] <- "CON"

phenotypes$Treatment <- NA
phenotypes$Treatment[phenotypes$Treatment_Condition == "Interferon Gamma"] <- "INF-gamma"
phenotypes$Treatment[phenotypes$Treatment_Condition == "Interleukin 1 Beta"] <- "IL1-beta"
phenotypes$Treatment[phenotypes$Treatment_Condition == "Vehicle"] <- "Vehicle"

# before anything else, make sure the reference group are set properly for the variables that will go
# into the model

phenotypes$Group <- as.factor(phenotypes$Group)
phenotypes$Group <- relevel(phenotypes$Group, ref = "CON")

phenotypes$Treatment <- as.factor(phenotypes$Treatment)
phenotypes$Treatment <- relevel(phenotypes$Treatment, ref = "Vehicle")

# define differential expression model

ddexpr_model <- ~ Group*Treatment + Age + Assigned_percent + (1|Individual_ID)

#re-generate voom object using "voomWithDreamWeights" instead

gene_voom <- voomWithDreamWeights(counts = dge_gene, formula = ddexpr_model, data = phenotypes)

#check potential colinearity problems in the model

diag_pheno <- phenotypes
diag_pheno$gene_expr <- gene_voom$E[1,]

diag_lm <- lm(gene_expr ~ Group*Treatment + Age + Assigned_percent, data = diag_pheno)

diag_vif <- vif(diag_lm)

print(paste("The DF-adjusted GVIF value for the variable of interest is", round(diag_vif[1,3], digits = 2), sep = " "))

if(diag_vif[1,3] > 2){
	print("Variance of coefficients inflated beyond acceptable levels")
} else if(diag_vif[1,3] < 2){
	print("No coefficient variance inflation problems with the model")
}

#fit model to gene expression data to estimate coefficients
# use the "dream" function rather than the lmFit function, as the former allows us to model the
# interindividual variability by adding the Individual IDs as a random effect

# We are working with very small sample size , so use the "Kenward-Roger" method 
# rather than the "Satterthwaite" approximation

gene_fit <- dream(gene_voom, formula = ddexpr_model, data = phenotypes, ddf = "Kenward-Roger", BBPARAM = param)

# #before building the contrasts, explore the landscape of coefficients; 
# that will give me an idea of which of them could be added and which not

coeffs <- gene_fit$coefficients
colnames(coeffs)[1] <- "Intercept"

# generate heatmap for the raw coefficients

png("heatmap_of_coefficients.png", width = 300, height = 900)
aheatmap(coeffs[,-1], scale = "column")
dev.off()

coeffs_cormatrix <- cor(coeffs[,-1], method = "pearson")

pdf("coefficient_profiles_cormatrix_plot.pdf", width = 30, height = 16, onefile = FALSE)
aheatmap(coeffs_cormatrix, scale = "none", txt = round(coeffs_cormatrix, digits = 2), labCol = NA, fontsize = 18, Rowv = NA, Colv = NA)
dev.off()

# set what coefficients we want to extract differential gene expression results for

deg_coefs <- colnames(gene_fit$coefficients)[-1]

# generate a vector with the names of the signatures

sign_names <- c("SCZ_vs_CON_in_vehicle", "IL1b_vs_vehicle_in_controls", "INFg_vs_vehicle_in_controls",
                "Age", "Assigned_percent", "IL1b_effect_scz_vs_controls", "INFg_effect_scz_vs_controls")

# create empty lists in which to store differential expression signatures

DEG_list <- list(length = length(deg_coefs))
FC_list <- list(length = length(deg_coefs))

# for visualizations, add rownames of the gene-fit as the GeneSymbols in the "genes" object

gene_fit$genes$GeneSymbol <- rownames(gene_fit$coefficients)

# source in the visualization functions

source("plotting_functions.R")

for(j in 1:length(deg_coefs)){

  # define "i" as "j+1"
  
  i <- j +  1 
  

  #extract DE signature and the whole logFC table for all the genes

  DEG_FCs <- topTable(gene_fit, coef = i, number = Inf, p.value = 1)
  DEGs <- DEG_FCs[DEG_FCs$adj.P.Val < 0.05, ]

    #write the tables out

  write.table(DEGs, paste("./DGE_results/", sign_names[j], "_significant_DGE_hits.txt", sep = ""), col.names = T, row.names = T, sep = "\t", quote = F)
  write.table(DEG_FCs, paste("./DGE_results/", sign_names[j], "_results_all_genes.txt", sep = ""), col.names = T, row.names = T, sep = "\t", quote = F)

  # print out name of signature and number of significant hits
  
  print(paste(nrow(DEGs), "significant genes for", sign_names[j]))
  
  # let's perform some visualizations


  pdf(paste("./DGE_results/", sign_names[j], "_volcanoplot.pdf", sep = ""))
  volcanoplot2(gene_fit, coef = i, highlight = round(nrow(DEGs) * 0.1), names = gene_fit$genes$GeneSymbol, 
               main = paste(sign_names[j], "signature", sep = "_"))
  dev.off()

  #store both Differential expression signatures and 

  DEG_list[[j]] <- DEGs 
  FC_list[[j]] <- DEG_FCs

  #assign the name of the comparison to the signature and the FC list

  assign(paste(sign_names[j], "signature", sep = "_"), DEGs)
  assign(paste(sign_names[j], "logFCs", sep = "_"), DEG_FCs)

}

#now prepare a summary plot to describe the sizes of the signatures and the distribution of up- and down-regulated genes

downregulated_genes <- numeric(length = length(sign_names))
upregulated_genes <- numeric(length = length(sign_names))

for (i in 1:length(sign_names)) {
	degs <- DEG_list[[i]]
	if (nrow(degs[degs$logFC < 0, ]) > 0){
	downregulated_genes[i] <- nrow(degs[degs$logFC < 0, ])
	} else {
		downregulated_genes[i] <- 0
	}
	if (nrow(degs[degs$logFC > 0, ]) > 0){
		upregulated_genes[i] <- nrow(degs[degs$logFC > 0, ])
	} else {
		upregulated_genes[i] <- 0
	}
}

comps <- length(sign_names)

de_genes <- c(downregulated_genes, upregulated_genes)
comparison <- rep(sign_names, times = 2)
deregulation <- c(rep("downregulation", times = comps), rep("upregulation", times = comps))
colors <- rep(c("blue", "red"), times = comps)
de_genes_df <- data.frame(comparison, deregulation, de_genes, colors)

# generate vector with font sizes

fontsizes <- log(de_genes_df$de_genes, 2)
fontsizes[fontsizes < 2] <- 2


pdf("./DGE_results/signature_lengths.pdf", width = 24, height = 16)

g1 <- ggplot(de_genes_df, aes(comparison, de_genes), ylab = "Number of DE genes") + 
  geom_bar(stat = "identity", position = "stack", aes(fill = deregulation)) + 
  scale_fill_manual(values = colors) + geom_text(aes(x = comparison, y = c(de_genes_df[(comps+1):nrow(de_genes_df), 
  "de_genes"], de_genes_df[1:comps, "de_genes"]), label = c(de_genes_df[(comps+1):nrow(de_genes_df), "de_genes"], 
  de_genes_df[1:comps, "de_genes"]), size = fontsizes), hjust = 0.5, vjust = 1.5, position = "stack", 
  fontface = "bold") + theme_bw() + theme(axis.text.x = element_text(face = "bold", family = "Helvetica", 
  size = 18, angle = 60, vjust = 1, hjust = 1), axis.title = element_text(face = "bold", family = "Helvetica", 
  size = 18), title = element_text(face = "bold", family = "Helvetica", size = 24), 
  axis.text.y = element_text(face = "bold", family = "Helvetica", size = 18)) 

print(g1)

dev.off()

# extract expression values for the receptors of IL1beta and INF-gamma

receptor_names <- c("IL1R", "IFNGR")

# get normalized log2-cpm expression matrix from the voom object

norm_expr <- gene_voom$E

# grep receptor names to subset their expression

recep_expr <- norm_expr[unique(unlist(lapply(receptor_names, function(x) grep(x, rownames(norm_expr))))),]

# add those expression values to the phenotypic table

phenos_recep <- merge(phenotypes, as.data.frame(t(recep_expr)), by = 0)

# in a loop, generate a control vs schizophrenia expression distribution boxplot for
# each of the 5 IL1B/INF-g receptor genes

recep_genes <- rownames(recep_expr)

# put this into long format

recep_genes_long <- melt(data = phenos_recep, id.vars = "Group", measure.vars = rownames(recep_expr))
colnames(recep_genes_long) <- c("Group", "Receptor", "log2CPM")

# generate boxplots using the facet_wrap function on the receptors

pdf("./DGE_results/scz_vs_con_expression_of_cytokine_receptors_violin_plots_free_scales.pdf",width = 9)
g2 <- ggplot(data = recep_genes_long, aes(x = Group, y = log2CPM)) +
  geom_violin(aes(fill = Group), alpha = 0.5) + ylab("Gene expression (log2CPMs)") + xlab("Clinical group") +
  geom_jitter(aes(colour = Group)) + theme_bw() + facet_wrap(~ Receptor, scales = "free") 

print(g2)
  
dev.off()


pdf("./DGE_results/scz_vs_con_expression_of_cytokine_receptors_violin_plots_one_scale.pdf", width = 9)
g2 <- ggplot(data = recep_genes_long, aes(x = Group, y = log2CPM)) +
  geom_violin(aes(fill = Group), alpha = 0.5) + ylab("Gene expression (log2CPMs)") + xlab("Clinical group") +
  geom_jitter(aes(colour = Group)) + theme_bw() + facet_wrap(~ Receptor) 

print(g2)

dev.off()

# save objects in working environment into RData file

save.image("differential_expression_from_RNAseq_using_DREAM.RData")
