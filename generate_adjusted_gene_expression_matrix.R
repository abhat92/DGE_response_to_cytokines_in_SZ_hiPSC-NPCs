#this is a script to extract residuals from gene expression data using a linear mixed model in limma and evaluate the performance of the correction with variancePartition

#upload necessary libraries

library(limma)
library(edgeR)
library(variancePartition)
library(doParallel)

## set the number of cores to run variancePartition in parallel

# registerDoParallel(28)

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

dge_gene <- dge_gene[rowSums(cpm(dge_gene) >1)>=0.01*(3/samples),,keep.lib.sizes = FALSE]
dge_gene <- calcNormFactors(dge_gene)

#generate model design matrix to be fed to voom and lmFit

global_design <- model.matrix(~ as.factor(IBD_BiopsyRegion) + as.factor(IBD_BiopsyType) + as.factor(disease_type) + as.numeric(combined_age) + as.factor(RNA_mol_gender) + as.factor(CT_RNA_extractionBatch) + as.numeric(CT_RQN), correct_pheno)

#apply voom-transformation on the filtered count matrix and output expression vs. variance plot to check if the distribution looks as it should

pdf("global_voom.pdf") 
voom <- voom(dge_gene, design =  global_design, plot = TRUE, span = 0.7)
dev.off()

norm_expr <- voom$E

#fit mixed linear model for each gene, extract variance of interest an generate adjusted matrices

Fit <- lmFit(voom, design = global_design)

#extract the coefficients from the fit object, sum the technical effects, on one hand, and the technical/demographic effect, on the other, and remove them from the expression data matrix

coeffs <- Fit$coef

tech_effects <- coeffs[,c(16:94)]%*%t(global_design[,c(16:94)])
tech_demo_effects <- coeffs[,c(14:94)]%*%t(global_design[,c(14:94)])
region_tech_demo_effects <- coeffs[,c(2:9,14:94)]%*%t(global_design[,c(2:9,14:94)])

tech_corr_expr <- norm_expr - tech_effects
tech_demo_corr_expr <- norm_expr - tech_demo_effects
region_tech_demo_corr_expr <- norm_expr - region_tech_demo_effects

#round residuals to get a lighter matrix

tech_corr_expr <- round(tech_corr_expr, digits = 2)
tech_demo_corr_expr <- round(tech_demo_corr_expr, digits = 2)
region_tech_demo_corr_expr <- round(region_tech_demo_corr_expr, digits = 2)

#write the adjusted matrices out

write.table(tech_corr_expr, "extraction_batch_RQN_corrected_matrix_with_intercepts", col.names = T, row.names = T, sep = "\t", quote = F)
write.table(tech_demo_corr_expr, "extraction_batch_RQN_gender_age_corrected_matrix_with_intercepts", col.names = T, row.names = T, sep = "\t", quote = F)
write.table(region_tech_demo_corr_expr, "extraction_batch_RQN_gender_age_biopsy_region_corrected_matrix_with_intercepts", col.names = T, row.names = T, sep = "\t", quote = F)

#now apply variance Partition and see how the effect of those variables has been reduced 

#it looks like numeric variables are not accepted as random effects in the version of variancePartition installed in Minerva;
#do as.factor() on them and see if that works

correct_pheno$GRID_2 <- as.factor(correct_pheno$GRID)
correct_pheno$CT_RNA_extractionBatch_2 <- as.factor(correct_pheno$CT_RNA_extractionBatch)

#include the variables of interest in the form to be fed to variancePartition

form = ~ (1|GRID_2) + (1|IBD_BiopsyRegion) + (1|IBD_BiopsyType) + (1|disease_type) + combined_age + (1|RNA_mol_gender) + CT_RQN + (1|CT_RNA_extractionBatch_2) + (1|seq_batch) + total_RNA_ug

#check whether both tables are in the same order in both tables

print(paste("Are the samples in the same order?", identical(colnames(tech_corr_expr), rownames(correct_pheno))))
print(paste("Are the samples in the same order?", identical(colnames(tech_demo_corr_expr), rownames(correct_pheno))))
print(paste("Are the samples in the same order?", identical(colnames(region_tech_demo_corr_expr), rownames(correct_pheno))))

#apply variance Partition to extract the fraction of gene expression variance explained by its explanatory variable for each gene

# tech_corr_expr <- tech_corr_expr[1:200,]
# tech_demo_corr_expr <- tech_demo_corr_expr[1:200,]
# region_tech_demo_corr_expr <- region_tech_demo_corr_expr[1:200,]

varPart_tech <- fitExtractVarPartModel(tech_corr_expr, form, correct_pheno)
varPart_tech_demo <- fitExtractVarPartModel(tech_demo_corr_expr, form, correct_pheno)
varPart_region_tech_demo <- fitExtractVarPartModel(region_tech_demo_corr_expr, form, correct_pheno)

#sort explanatory variables based on the median variance explained (from high to low); this ordering reflects the explanatory power

sortVarPart_tech <- sortCols(varPart_tech)
sortVarPart_tech_demo <- sortCols(varPart_tech_demo)
sortVarPart_region_tech_demo <- sortCols(varPart_region_tech_demo)

#output variance Partition plot

pdf("extract_batch_RQN_corrected_VarPartPlot.pdf", width = 10, height = 8)
print(plotVarPart(sortVarPart_tech , label.angle = 45) + theme(text=element_text(size=20), axis.text.x=element_text(size = 18)))
dev.off()

pdf("extract_batch_RQN_gender_age_corrected_VarPartPlot.pdf", width = 10, height = 8)
print(plotVarPart(sortVarPart_tech_demo, label.angle = 45) + theme(text=element_text(size=20), axis.text.x=element_text(size = 18)))
dev.off()

pdf("extract_batch_RQN_gender_age_biopsy_region_corrected_VarPartPlot.pdf", width = 10, height = 8)
print(plotVarPart(sortVarPart_region_tech_demo, label.angle = 45) + theme(text=element_text(size=20), axis.text.x=element_text(size = 18)))
dev.off()


#save image of object space

save.image("gene_expression_correction.RData")


 

