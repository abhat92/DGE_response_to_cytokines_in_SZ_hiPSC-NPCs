# give appropriate format to gene annotations file

rm(list = ls())

# read in GTF file

gene_gtf <- read.table("Homo_sapiens.GRCh38.99.genes.gtf", header = F, sep = "\t")

# split the ninth column

split_ninth <- lapply(gene_gtf$V9, function(y) unlist(strsplit(x = as.character(y), split = "; ")))

# put all the vectors into a data frame

ninth_df <- as.data.frame(do.call("rbind", split_ninth))

# remove unnecessary columns

ninth_df <- ninth_df[,c("V1", "V3", "V5")]

# rename columnes

colnames(ninth_df) <- c("gene_id", "gene_name", "gene_biotype")

# remove the annotation names from the entries

ninth_df$gene_id <- gsub("gene_id ", "", ninth_df$gene_id)
ninth_df$gene_name <- gsub("gene_name ", "", ninth_df$gene_name)
ninth_df$gene_biotype <- gsub("gene_biotype ", "", ninth_df$gene_biotype)
ninth_df$gene_biotype <- gsub(";", "", ninth_df$gene_biotype)

gene_annots <- as.data.frame(cbind(gene_gtf[,-9], ninth_df))

# re-organize and select columns

gene_annots <- gene_annots[,c("gene_id", "V1", "V4", "V5", "V7", "gene_name", "gene_biotype")]
colnames(gene_annots)[2:5] <- c("Chr", "start", "end", "strand")

# identify the gene-names that have more than one entry

gene_freqs <- table(gene_annots$gene_name)
dupli_genes <- names(gene_freqs)[gene_freqs > 1]

# open for loop

# generate empty list in which to store collapsed annotations

collap_list <- vector("list", length = length(dupli_genes))

for(i in 1:length(dupli_genes)){
  
  # extract annotations for corresponding gene
  
  dupli_annots <- gene_annots[gene_annots$gene_name == dupli_genes[i],]
  
  # create empty vector in which to store collapsed values
  
  collap_vec <- vector("character", length = ncol(dupli_annots))
  
  for(j in 1:(ncol(dupli_annots))){
    
    # extract values for the corresponding column
    
    col_vals <- as.character(unique(dupli_annots[,j]))
    
    if(length(col_vals) == 1){
      
      collap_vec[j] <- col_vals
      
    } else if (length(col_vals) > 1){
      
      collap_vec[j] <- paste(col_vals, collapse = ";")
    }
    
  }
  
  collap_list[[i]] <- collap_vec
}

# collapse into a data frame

collap_df <- as.data.frame(do.call("rbind", collap_list))
colnames(collap_df) <- colnames(gene_annots)

# remove the duplicated annotations and add the collapsed ones

gene_annots_2 <- gene_annots[which(is.na(match(gene_annots$gene_name, dupli_genes))), ]
gene_annots_3 <- as.data.frame(rbind(gene_annots_2, collap_df))

# let's see if there are any duplicated gene-names left

gene_freqs_2 <- table(gene_annots_3$gene_name)
dupli_genes_2 <- names(gene_freqs_2)[gene_freqs_2 > 1]
  
print(paste("There are", length(dupli_genes_2), "duplicated gene names left in the gene annotations"))
  
# add gene- names as rownames

rownames(gene_annots_3) <- gene_annots_3$gene_name

# write table out

write.table(gene_annots_3, "Ensembl_Hg38_v99_gene_annotations.txt", col.names = T, row.names = F, sep = "\t", quote = F)
