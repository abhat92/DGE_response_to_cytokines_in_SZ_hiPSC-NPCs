#upload libraries

library(ggplot2)
library(ggrepel)
# 
# #read mapping stats extract from the final STAR log files
# 
# map_stats <- read.table("mapping_stats_table", header = F, row.names = 1, sep = "\t")
# colnames(map_stats) <- c("TotalReads", "UniquelyMapped", "MultiMapped")
# map_stats$TotalMapped <- apply(map_stats[,2:3], 1,sum)
# 
# #generate mapping stat distribution plots
# 
# pdf("mapped_read_count_distributions.pdf")
# plot(density(map_stats$TotalReads), col = "red")
# lines(density(map_stats$UniquelyMapped), col = "blue")
# lines(density(map_stats$MultiMapped), col = "green")
# legend("topright", legend = c("TotalReads", "UniquelyMapped", "MultiMapped"), col = c("red", "blue", "green"), lwd = 1)
# g <- ggplot(data = map_stats, aes(x = TotalReads, y = UniquelyMapped)) + geom_point() + geom_text(aes(label = rownames(map_stats))) + theme_bw()
# print(g)
# dev.off()
# 
# #upload featureCounts counting summary stats table
# 
# count_stats <- read.table("counting_stats_table", header = F, row.names = 1, sep = " ")
# colnames(count_stats) <- c("Assigned", "UnassignedAmbiguity", "UnassignedMultimapping", "UnassignedNoFeatures")
# 
# #compare number of reads aligned by STAR and number of reads assigned by featureCounts
# 
# map_count_stats <- merge(map_stats, count_stats, by.x = 0, by.y = 0)
# colnames(map_count_stats)[1] <- "Sample"
# 
# write.table(map_count_stats, 'alignment_counting_stats.txt', col.names = T, row.names = F, sep = "\t", quote = F)
# 
# pdf("aligned_vs_assigned_reads.pdf")
# g <- ggplot(data = map_count_stats, aes(x = TotalMapped, y = Assigned)) + geom_point() + geom_text(aes(label = Sample)) + theme_bw() + geom_abline(slope= 1, intercept = c(0,0), color = "red")
# print(g)
# dev.off()


#read raw count matrix in

raw_counts <- read.table("raw_count_matrix", header = T, row.names = 1, sep = " ")

#turn column names from paths to sample names

colnames(raw_counts) <- gsub(".*ANB01", "ANB01", colnames(raw_counts))
colnames(raw_counts) <- gsub("Aligned.*", "", colnames(raw_counts))

#reorder columns and row of raw_count matrix

raw_counts <- raw_counts[order(rownames(raw_counts)), order(colnames(raw_counts))]

#write raw count matrix out and read it in again

write.table(raw_counts, "raw_count_matrix_2", col.names = T, row.names = T, sep = "\t")

# raw_counts_2 <- read.table("raw_count_matrix_2", header = T, row.names = 1, sep = "\t")

## reorder raw count matrix columns

# raw_counts_2 <- raw_counts_2[,order(colnames(raw_counts_2))]

#read phenotypic table in 

phenotypes <- read.table("Phenotype_Table_RNASeq.txt", header = T, row.names = 1, sep = "\t")

#chekc what sample names are different in the phenotypic table and count matrix

diff_samples <- setdiff(rownames(phenotypes), colnames(raw_counts))

#re-order both tables and check if all the samples included in the phenotypic table are included in the gcounts matrix

raw_counts <- raw_counts[,order(colnames(raw_counts_2))]
phenotypes <- phenotypes[order(rownames(phenotypes)),]

write.table(raw_counts_2, "gcount_186_samples_april_24th_2017", col.names = T, row.names = T, sep = "\t", quote = F)

#extract counts for samples described in the phenotypic table

gcounts <- raw_counts_2[, rownames(phenotypes)]

print(paste("Are the samples in the same order?", identical(colnames(gcounts), rownames(phenotypes)), sep = " "))

write.table(gcounts, 'gcounts_158_samples_april_24th_2017', col.names = T, row.names = T, sep = "\t", quote = F)


#do some sanity check

#distribution of library sizes

lib_sizes <- colSums(raw_counts)
lib_size_den <- density(lib_sizes)

lib_sizes <- as.data.frame(lib_sizes)
lib_sizes$Sample <- rownames(lib_sizes)
colnames(lib_sizes) <- c("Reads", "Sample")

pdf("library_size_distribution.pdf", height = 8, width = 16)
plot(lib_size_den)
g <- ggplot(data = lib_sizes, aes(x = Sample, y = as.numeric(lib_sizes$Reads))) + geom_col(fill = "white", colour = "black") + geom_text_repel(aes(y = as.numeric(lib_sizes$Reads), label = Sample), angle = 90) + theme_bw() + ylab("Reads assigned to features")
print(g)
dev.off()

#check that samples from the same individual do not have the exact same count profiles

sample_freq <- table(phenotypes$Individual)
double_ind <- unique(names(sample_freq[sample_freq == 2]))

pdf("intraindividual_count_plots.pdf")
for(i in 1:length(double_ind)){
	ind_pheno <- phenotypes[phenotypes$Individual == double_ind[i],]
	plot(log(gcounts[,rownames(ind_pheno)[1]],10), log(gcounts[,rownames(ind_pheno)[2]],10), xlab = paste(ind_pheno[1,"Individual"], ind_pheno[1,"Cell.type"], sep = "_"), ylab = paste(ind_pheno[2,"Individual"], ind_pheno[2,"Cell.type"], sep = "_"))
}
dev.off()

save.image("count_matrix_preparation.RData")
