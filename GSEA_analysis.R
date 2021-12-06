# clean working environment

rm(list = ls())

# load necessary libraries

library(qusage)
library(org.Hs.eg.db)
library(fgsea)
library(data.table)
library(NMF)
library(ggplot2)

# read in GMT file with gene-sets

gs_list <- read.gmt("./cleaned_CNS_and_immune_gene_sets_2020_04_27.gmt")

# get list of files that contain the differential expression signatures we have generated

sign_fnames <- list.files("./DGEs_3/DGEs_2_two_covariates")[grep("all_genes.txt", list.files("./DGEs_3/DGEs_2_two_covariates"))]

# we are not interested in performing a gene-set enrichment analysis on the genes that are 
# differentially expressed with age or due to technical effects, so we will remove those
# filenames from the list
#sign_fnames <- sign_fnames[-unlist(lapply(c("Assigned"), 
                                          #function(x) grep(x, sign_fnames)))]

# generate signature names from the file names

sign_names <- gsub("_results.*", "", sign_fnames)

# now, in a loop, we will iterate throught the signatures to perform a GSEA analysis

for(i in 1:length(sign_fnames)){
  
  # read in signature
  
  sign <- read.table(paste("./DGEs_3/DGEs_2_two_covariates/", sign_fnames[i], sep = ""), header = T, sep = "\t")
  
  # translate the GeneSymbols into Entrez IDs
  
  sign$EntrezID <- mapIds(x = org.Hs.eg.db, keys = as.character(sign$GeneSymbol),
                          column = 'ENTREZID', keytype = 'SYMBOL')
  
  # remove those entries without an EntrezID
  
  na_entries <- which(is.na(sign$EntrezID))
  print(paste(length(na_entries), "out of", nrow(sign), "genes don't have an EntrezID"))
  print("Removing...")
  
  sign_2 <- sign[-na_entries,]
  
  # re-order the signature by the "z.std" value (this is comparable across genes; not the t-value,
  # in this case as, the way DREAM does the analysis, each gene has different degrees of freedom)

  sign_2 <- sign_2[order(sign_2$z.std, decreasing = TRUE),]
  
  # extract "z.std" values and add EntrezID names to the vector
  
  z_std <- sign_2$z.std
  names(z_std) <- sign_2$EntrezID
  
  # perform GSEA analysis
    
  gsea_res <- fgsea(pathways = gs_list, stats = z_std, nperm = 100000, minSize = 5, nproc = 3)
 
  # extract database and GeneSet names
  
  gsea_res$database <- gsub("_.*", "", gsea_res$pathway)
  gsea_res$GeneSet <- unlist(lapply(gsea_res$pathway, 
            function(x) gsub(paste(paste(unique(gsea_res$database),
                          "_", sep = ""), collapse = "|"), "", x)))
  
  # remove underscores from GeneSet names (they make Cytoscape think there is
  # only word)
  
  gsea_res$GeneSet <- gsub("_", " ", gsea_res$GeneSet)
    
  # calculate -log10p-value
  
  gsea_res$log10padj <- -log(gsea_res$padj, 10)
  
  # subset GSEA results for the significant ones
  
  sig_gsea_res <- gsea_res[gsea_res$padj < 0.05,]
  
  # write results out
  
  fwrite(as.data.frame(gsea_res), paste("./GSEA_results/", sign_names[i],
                              "_GSEA_enrichment_all_results.txt", sep = ""), 
              col.names = TRUE, row.names = F, sep = "\t", quote = F)
  fwrite(as.data.frame(sig_gsea_res), paste("./GSEA_results/", sign_names[i],
                                        "_GSEA_enrichment_significant_results.txt", sep = ""), 
         col.names = TRUE, row.names = F, sep = "\t", quote = F)
  
  # now, calculate pairwise Jaccard similarity between the gene-sets
  # that are significantly enriched in the signature
  
  # extract the list of gene-sets
  
  sig_gs_list <- gs_list[sig_gsea_res$pathway]
  
  # calculate all pairwise Jaccard similarity values and put them in a
  # symmetric matrix
  
  jacc_dist <- do.call("rbind", lapply(sig_gs_list, 
          function(x) unlist(lapply(sig_gs_list, 
          function(y) length(intersect(x,y))/length(union(x,y))))))
  
  # transform the matrix into a 3-column data frame from which we can
  # build a network in Cytoscape
  
  jacc_net <- data.frame(node1 = unlist(lapply(rownames(jacc_dist), 
                function(x) rep(x, times = ncol(jacc_dist)))), 
                         node2 = rep(colnames(jacc_dist), times = nrow(jacc_dist)), 
                jaccard = c(jacc_dist))
  
  # remove the entries that contain the Jaccard values for each GeneSet
  # with itself
  
  jacc_net <- jacc_net[-which(jacc_net$node1 == jacc_net$node2),]
  
  # and set values below 0.5 to 0
  
  jacc_net$jaccard[jacc_net$jaccard <= 0.5] <- 0
  
  # write data frame with network out
  
  write.table(jacc_net, paste("./GSEA_results/", sign_names[i],
                              "_enriched_GeneSet_network.txt", sep = ""), 
              col.names = TRUE, row.names = F, sep = "\t", quote = F)
  
}

#########################################
## CLUSTERING OF SIGNIFICANT GENE-SETS ##
#########################################

# I should have store the GSEA results in lists but I didn't, so I will have
# to read them in

# now do store the results in a list

gsea_res_all <- vector("list", length = length(sign_names))

for(i in 1:length(sign_names)){

  # read in significant GSEA results
  
  gsea_res <- read.table(paste("./GSEA_results/", sign_names[i],
                               "_GSEA_enrichment_all_results.txt", sep = ""),
                         header = T, sep = "\t", quote = "")
  
  # and store them in corresponding slot of list
  
  gsea_res_all[[i]] <- gsea_res
  
  # subset GSEA results for the significant ones
  
  sig_gsea_res <- gsea_res[gsea_res$padj < 0.05,]
  
  # print out number of signifiant GeneSets obtained for the signature
  
  print(paste(nrow(sig_gsea_res), 
              "enriched Genesets identified for", sign_names[i]))
  
  # now, calculate pairwise Jaccard similarity between the gene-sets
  # that are significantly enriched in the signature
  
  # extract the list of gene-sets
  
  sig_gs_list <- gs_list[names(gs_list) %in% sig_gsea_res$pathway]
  
  # calculate all pairwise Jaccard similarity values and put them in a
  # symmetric matrix
  
  jacc_dist <- do.call("rbind", lapply(sig_gs_list, 
                  function(x) unlist(lapply(sig_gs_list, 
                  function(y) length(intersect(x,y))/length(union(x,y))))))
  
  # transform the matrix into a 3-column data frame from which we can
  # build a network in Cytoscape
  
  jacc_net <- data.frame(node1 = unlist(lapply(rownames(jacc_dist), 
                                               function(x) rep(x, times = ncol(jacc_dist)))), 
                         node2 = rep(colnames(jacc_dist), times = nrow(jacc_dist)), 
                         jaccard = c(jacc_dist))
  
  # remove the entries that contain the Jaccard values for each GeneSet
  # with itself
  
  jacc_net <- jacc_net[-which(jacc_net$node1 == jacc_net$node2),]
  
  # get fraction of similarities above 0.5
  
  high_frac <- nrow(jacc_net[jacc_net$jaccard > 0.5,])/(nrow(jacc_net)-nrow(jacc_dist))
  
  print(paste(round(high_frac, digits = 4)*100, 
              "% of Geneset-pairs have Jaccard similarity above 0.5"))
  
  # generate heatmap with matrix
  
  pdf(paste("./GSEA_results/", sign_names[i], 
            "_significant_GeneSet_similarity_matrix.pdf", sep = ""),
      onefile = FALSE)
  
  aheatmap(jacc_dist, scale = "none", labRow = NA, labCol = NA)
  
  dev.off()

  # now I need to use this to find clusters of Gene-Sets
  # perform a hierarchical clustering and plot the result
  
  hclust_res <- hclust(d = as.dist(1-jacc_dist), method="complete")
  
  # open PDF device
  
  pdf(paste("./GSEA_results/", sign_names[i], 
            "_hclust_of_significant_GeneSets.pdf", sep = ""))
  
  # generate plot of hirearchical clustering
  
  plot(hclust_res, labels = FALSE)
  
  # close PDF device
  
  dev.off()
  
  # extract clusters
  
  clusters <- cutree(tree = hclust_res, h = 0.5)
  
  # print out number of clusters obtained for the signature
  
  print(paste(max(clusters), 
              "clusters of enriched Genesets identified for", sign_names[i]))
  
  # add clusters to table of significant results
  
  sig_gsea_res_2 <- merge(sig_gsea_res, clusters, by.x = "pathway", by.y = 0)
  colnames(sig_gsea_res_2)[ncol(sig_gsea_res_2)] <- "cluster"
  
  # now, for each cluster do the following:
  # - select the result with the lowest p-value
  # - collapse the names of the rest of the pathways in the cluster
  # - into an onlye string
  
  # generate empty list in which to store the results
  
  top_gsea_list <- vector("list", length = max(clusters))
  
  for(j in 1:max(clusters)){
    
    # extract enrichment results for corresponding cluster of significantly
    # enriched Gene-Sets
    
    clust_gsea_all <- sig_gsea_res_2[sig_gsea_res_2$cluster == j,]
    
    # if there is only 1 GeneSet in cluster, no need to select
    
    if(nrow(clust_gsea_all) == 1){
      
      # add top GeneSet and the rest of GeneSets
      
      clust_gsea_all$top_GeneSets <- clust_gsea_all$GeneSet
      clust_gsea_all$other_GeneSets_in_cluster <- NA
      clust_gsea_all$GeneSet_N_in_cluster <- nrow(clust_gsea_all)
      
      # add that table to corresponding slot in list
      
      top_gsea_list[[j]] <- clust_gsea_all
      
    } else if (nrow(clust_gsea_all) > 1){
      
      # select top GeneSet
      
      clust_gsea_top <- clust_gsea_all[clust_gsea_all$pval == min(clust_gsea_all$pval),]
      
      # if there are more than 1 GeneSet with the lowest p-value, select one row
      
      if(nrow(clust_gsea_top) >1) {
        
        clust_gsea_top_2 <- clust_gsea_top[1,]
        clust_gsea_top_2$top_GeneSets <- paste(clust_gsea_top$GeneSet, collapse = " ; ")
        
        # get the list of the rest of GeneSets
        
        other_gs <- setdiff(clust_gsea_all$GeneSet, clust_gsea_top$GeneSet)
        
        clust_gsea_top_2$other_GeneSets_in_cluster <- paste(other_gs,
                                                  collapse = " ; ")
        clust_gsea_top_2$GeneSet_N_in_cluster <- nrow(clust_gsea_all)
        
        # store table in list
        
        top_gsea_list[[j]] <- clust_gsea_top_2
     
      } else if (nrow(clust_gsea_top) == 1){
        
        clust_gsea_top$top_GeneSets <- clust_gsea_top$GeneSet
        
        # get the list of the rest of GeneSets
        
        other_gs <- setdiff(clust_gsea_all$GeneSet, clust_gsea_top$GeneSet)
        
        clust_gsea_top$other_GeneSets_in_cluster <- paste(other_gs,
                                                            collapse = " ; ")
        clust_gsea_top$GeneSet_N_in_cluster <- nrow(clust_gsea_all)
        
        # store table in list
        
        top_gsea_list[[j]] <- clust_gsea_top
         
        }
      }
    }
     
   # put the results across clusters together
  
  clumped_gsea_res <- as.data.frame(do.call("rbind", top_gsea_list))
  
  # select the columns that will be included in the final output
  
  clumped_gsea_res_2 <- clumped_gsea_res[,c("top_GeneSets", "cluster", "GeneSet_N_in_cluster", "size",
                                            "pval", "padj", "ES", "NES", 
                                            "other_GeneSets_in_cluster")]
  
  # rename columns and reorder the table according to p-value
  
  colnames(clumped_gsea_res_2) <- c("Top GeneSets", "Cluster", "GeneSets in cluster",
                                    "Genes in GeneSet", "p_value", "FDR", 
                                    "Enrichment Score", "Normalized Enrichment Score",
                                    "Other GeneSets in cluster")
  clumped_gsea_res_2 <- clumped_gsea_res_2[order(clumped_gsea_res_2$p_value),]
  
  # write table out as an csv file
  
  write.csv(clumped_gsea_res_2, paste("./GSEA_results/", sign_names[i],
                  "_GSEA_enrichment_clumped_significant_results.csv", sep = ""),
            row.names = F, quote = F)
  
  # we also want to visualize the results for the top 10 genesets (or, rather
  # the top GeneSet from the top 10 GeneSet clusters)
  
  # select top 10 results and calculate -log10p-value
  
  clumped_top_10 <- clumped_gsea_res_2[1:10,]
  clumped_top_10$logPval <- -log(clumped_top_10$p_value, 10)
  colnames(clumped_top_10) <- gsub(" ", "", colnames(clumped_top_10))
  
  # insert line break for those cases where more than 1 top GeneSet was identified
  
  clumped_top_10$TopGeneSets <- gsub(" ; ", "\n", clumped_top_10$TopGeneSets)
  
  # open PDF device
  
  pdf(paste("./GSEA_results/top_10_", 
            sign_names[i], "_enriched_GeneSets.pdf", sep = ""),
      width = 12)
  # generate ggplot object
  
  gsea_g <- ggplot(data = clumped_top_10, aes(x = logPval, y = TopGeneSets)) +
    geom_point(aes(colour = NormalizedEnrichmentScore, size = logPval)) + 
    theme_bw() +
    xlab("Enrichment Significance (-log10Palue)") + 
    theme(axis.title.y = element_blank()) + labs(colour = "NES", size = "-log10Pvalue") +
    scale_y_discrete(limits = rev(clumped_top_10$TopGeneSets)) +
    scale_colour_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0) + ggtitle(gsub("_"," ", sign_names[i]))
  print(gsea_g)
  
  # close PDF device
  
  dev.off()
}


# save object in working environment into an .RData file

save.image("./GSEA_results/GSEA_analysis.RData")

