# simple script to explore source of variation in RNAseq gene expression data

# clean working environment

rm(list = ls())

#upload necessary libraries

library(limma)
library(edgeR)
library(variancePartition)
library(doParallel)
library(FactoMineR)
library(sva)

#upload the phenotype-table and gcount matrix

phenotypes <- read.table("Phenotype_Table_RNASeq.txt", header = T, row.names = 1, sep = "\t", quote = "")
gcounts <- read.table("raw_count_matrix_2", header = T, row.names = 1, sep = "\t")

#check if samples are in the same order

print(paste("Are the samples in the same order?", identical(colnames(gcounts), rownames(phenotypes))))

# Normalize RNA seq count data

samples <- nrow(phenotypes)
dge_gene <- DGEList(counts = gcounts)

# calculate the counts-per-million

cpms <- cpm(dge_gene)

# plot out the cpm distribution

pdf("cpm_density_distribution.pdf")

plot(density(as.numeric(log(cpms, 10))))

dev.off()

dge_gene <- dge_gene[rowSums(cpm(dge_gene) >10^0.6)>=0.01*(3/samples),,keep.lib.sizes = FALSE]
dge_gene <- calcNormFactors(dge_gene)

pdf("filtered_cpm_density_distribution.pdf")

plot(density(as.numeric(log(dge_gene@.Data[[1]], 10))))

dev.off()


#generate model design matrix to be fed to voom and lmFit

global_design <- model.matrix(~ Case_control + Treatment_Condition + Age + Assigned_percent, phenotypes)

#apply voom-transformation on the filtered count matrix and output expression vs. variance plot to check if the distribution looks as it should

pdf("global_voom.pdf") 
voom <- voom(dge_gene, design =  global_design, plot = TRUE, span = 0.7)
dev.off()

norm_expr <- voom$E


##################
## PCA ANALYSIS ##
##################

#perform Principal Component Analysis on the expression data

pca_res <- PCA(X = norm_expr, scale.unit = TRUE, graph = FALSE)

#extract coordinates of sample across principal components

pca_coords <- pca_res$var$coord

#extract percentages of explained variance for the first 3 PCs

pc1_expl <- pca_res$eig[1,2]
pc2_expl <- pca_res$eig[2,2]
pc3_expl <- pca_res$eig[3,2]

#merge sample coordinates with phenotypes

phenos_coords <- merge(phenotypes, pca_coords, by.x = 0, by.y = 0)
rownames(phenos_coords) <- phenos_coords$Row.names
phenos_coords <- phenos_coords[,-1]

#now, generate PC1 vs PC2 and PC2 vs PC3 plots for each phenotype (using the phenotype as coloring variable)

pheno_list <- colnames(phenos_coords)[1:(grep("Dim.1", colnames(phenos_coords))-1)]

#loop through the phenotypes that need to be color-coded in each PCA plot set

#open a for loop

for(i in 1:length(pheno_list)){
  
  #open pdf device
  
  pdf(paste(pheno_list[i], "pca_plots.pdf", sep = "_"), width = 10, height = 8)
  
  #extract column that contains the data of the corresponding phenotype and turn it into a factor
  
  coloring <- phenos_coords[, pheno_list[i]]  
  
  #plot differently depending on whether the coloring vectors is a factor or a numeric vector
  
  if(class(coloring) == "factor") {
    
    print("Plotting a factor as color-code")
    
    #generate the ggplot objects with all the necessary layers
    
    #PC1 vs PC2
    
    g1 <- ggplot(data = phenos_coords, aes(x = Dim.1, y = Dim.2)) + geom_point(aes(colour = coloring), size = 3) + xlab(paste("PC1 ( ", as.character(round(pc1_expl, digits = 2)), "% )", sep = "")) + ylab(paste("PC2 ( ", as.character(round(pc2_expl, digits = 2)), "% )", sep = "")) + ggtitle(paste("PC1 vs PC2 colored by", pheno_list[i], sep = " ")) + stat_ellipse(aes(fill = coloring), geom = "polygon", alpha = 0.25) + theme_bw()
    
    #print ggplot object out
    
    print(g1)
    
    #PC2 vs PC3
    
    g2 <- ggplot(data = phenos_coords, aes(x = Dim.2, y = Dim.3)) + geom_point(aes(colour = coloring), size = 3) + xlab(paste("PC2 ( ", as.character(round(pc2_expl, digits = 2)), "% )", sep = "")) + ylab(paste("PC3 ( ", as.character(round(pc3_expl, digits = 2)), "% )", sep = "")) + ggtitle(paste("PC2 vs PC3 colored by", pheno_list[i], sep = " ")) + stat_ellipse(aes(fill = coloring), geom = "polygon", alpha = 0.25) + theme_bw()
    
    #print ggplot object out
    
    print(g2)
    
  } else if (class(coloring) == "integer" || class(coloring) == "numeric") {
    
    print("Plotting an integer as color-code")
    
    #generate the ggplot objects with all the necessary layers
    
    #PC1 vs PC2
    
    g1 <- ggplot(data = phenos_coords, aes(x = Dim.1, y = Dim.2)) + geom_point(aes(fill = coloring), colour = "black", pch = 21, size = 3) + xlab(paste("PC1 ( ", as.character(round(pc1_expl, digits = 2)), "% )", sep = "")) + ylab(paste("PC2 ( ", as.character(round(pc2_expl, digits = 2)), "% )", sep = "")) + ggtitle(paste("PC1 vs PC2 colored by", pheno_list[i], sep = " ")) + scale_fill_gradient(low = "white", high = "red") + theme_bw()
    
    #print ggplot object out
    
    print(g1)
    
    #PC2 vs PC3
    
    g2 <- ggplot(data = phenos_coords, aes(x = Dim.2, y = Dim.3)) + geom_point(aes(fill = coloring), colour = "black", pch = 21, size = 3) + xlab(paste("PC2 ( ", as.character(round(pc2_expl, digits = 2)), "% )", sep = "")) + ylab(paste("PC3 ( ", as.character(round(pc3_expl, digits = 2)), "% )", sep = "")) + ggtitle(paste("PC2 vs PC3 colored by", pheno_list[i], sep = " ")) + scale_fill_gradient(low = "white", high = "red") + theme_bw()
    
    #print ggplot object out
    
    print(g2) 
    
  }
  
  #close PDF device
  dev.off()
  
  #close the for loop
  
}

##################################
## IDENTIFY SURROGATE VARIABLES ##
##################################

# use the sva package to identify surrogate variables

# build model matrix from know variables that are potential sources of variation of gene expression

mod0 <- model.matrix(~ 1, data = phenotypes)

# ran sva to identify hidden sources of variance as surrogate variables

surro_vars <- sva(as.matrix(norm_expr), mod = global_design, mod0 = mod0, method = "irw", B = 100)

sv_mat <- surro_vars$sv
colnames(sv_mat) <- paste("sva", seq(from = 1, to = ncol(sv_mat)), sep = "_")

# add significant surrogate variables to phenotypic table

phenotypes_2 <- as.data.frame(cbind(phenotypes, sv_mat))

########################
## VARIANCE PARTITION ##
########################

## write down form that describes the model for variancePartition

form = ~ (1|Individual_ID) + (1|Case_control) + (1|Treatment_Condition) + Age + Assigned_percent + sva_1 + sva_2 + sva_3

#apply variance Partition to extract the fraction of gene expression variance explained by each explanatory variable for each gene

varPart <- fitExtractVarPartModel(exprObj = voom, formula = form, data = phenotypes_2)

#sort explanatory variables based on the median variance explained (from high to low); this ordering reflects the explanatory power

sortVarPart <- sortCols(varPart)

#output variance Partition plot

pdf("VarPartPlot.pdf", width = 10, height = 8)
print(plotVarPart(sortVarPart, label.angle = 45) + theme(text=element_text(size=20), axis.text.x=element_text(size = 18)))
dev.off()

## WITHOUT SURROGATE VARIABLES

## write down form that describes the model for variancePartition

form_2 = ~ (1|Individual_ID) + (1|Case_control) + (1|Treatment_Condition) + Age + Assigned_percent

#apply variance Partition to extract the fraction of gene expression variance explained by each explanatory variable for each gene

varPart_2 <- fitExtractVarPartModel(exprObj = voom, formula = form_2, data = phenotypes_2)

#sort explanatory variables based on the median variance explained (from high to low); this ordering reflects the explanatory power

sortVarPart_2 <- sortCols(varPart_2)

#output variance Partition plot

pdf("VarPartPlot_2.pdf", width = 10, height = 8)
print(plotVarPart(sortVarPart_2, label.angle = 45) + theme(text=element_text(size=20), axis.text.x=element_text(size = 18)))
dev.off()

# save objects in .RData file

save.image("sources_of_variation_in_NPC_RNAseq.RData")
