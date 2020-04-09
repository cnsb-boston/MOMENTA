

###########################
#
## Run MOMENTA Enrichments using fgsea
#
##########################


##########################
# 
#  Variables

working_dir <- file.path("some", "path")

rnk_file <- file.path("some", "file")

species <- "Human (9606)"
  # One of "Human (9606)", "Mouse (10090)", "Fruit Fly (7227)", "Yeast (559292)", "E. coli (511145)"

# 
##########################


##########################
# 
#  Load required packages
# 

packages <- c("fgsea") 

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(setdiff(packages, rownames(installed.packages())))  
}  

library(fgsea)

##########################


##########################
# 
# 
# 

# Set working directory
setwd(working_dir)

# Download genesets file
git_url <- "https://raw.githubusercontent.com/cnsb-boston/MOMENTA/master"
if( dir.exists("genesets") == FALSE ) { dir.create("genesets") }
if( file.exists(file.path("genesets","MatchedGeneSets.txt"))==FALSE ){
  download.file(url=file.path(git_url, "genesets", "MatchedGeneSets.txt"),
            destfile= file.path("genesets","MatchedGeneSets.txt") )
}
geneset_lookup <-  read.delim(file.path("genesets","MatchedGeneSets.txt"))
gmt_files <- as.vector(geneset_lookup[which(geneset_lookup[,"Species"]==species),"GeneSet"])
gmt_names <- as.vector(geneset_lookup[which(geneset_lookup[,"Species"]==species),"Name"])

# Set up output dirs
gsea_working_path <- paste("MOMENTA_", Sys.Date(), sep="")
gsea_images_path <- file.path(gsea_working_path, "Images")
gsea_absval_path <- file.path(gsea_working_path, "AbsoluteValueRank")
if( dir.exists(gsea_working_path) == FALSE ) { dir.create(gsea_working_path) }
if( dir.exists(gsea_images_path) == FALSE ) { dir.create(gsea_images_path) }
if( dir.exists(gsea_absval_path) == FALSE ) { dir.create(gsea_absval_path) }

# Run enrichment for each GMT file
for(gmt_index in 1:length(gmt_files)){
  
  dest_gmt_file <-file.path("genesets", gmt_files[gmt_index])
  if( file.exists(dest_gmt_file)==FALSE ){
    download.file(url=file.path(git_url, "genesets", gmt_files[gmt_index]),
                  destfile= dest_gmt_file )
  }

  analysis_name<- paste(gsub(".rnk","",rnk_file), "_", gmt_names[gmt_index], sep="") ;

  for (i in 1:length(analysis_names)){
    
    ranked_features <- read.delim(rnk_file)
    ranked_vector <- ranked_features[,"rank"]
    names(ranked_vector) <- ranked_features[,"GeneName"]
    
    pathways_gmt <- gmtPathways(dest_gmt_file)
    suppressWarnings({
      fgsea_results <- fgsea(pathways=pathways_gmt, 
                             stats=ranked_vector,
                             minSize=15, 
                             maxSize=500, 
                             nperm=10000)
    })
    
    fgsea_out <- fgsea_results[,c(1,1,2,3)]
    colnames(fgsea_out) <- c("Term", "Description", "p.Val", "FDR")
    fgsea_out[,"Phenotype"] <- sign(fgsea_results[,"NES"])
    fgsea_out[,"Genes"] <- apply(fgsea_results[,"leadingEdge"], 1, function(x) paste(as.character(unlist(x["leadingEdge"])), collapse=","))
    fgsea_out[,"Genes"] <- apply(fgsea_out, 1, function(x) gsub(" ", "", x["Genes"]))
    fgsea_out[,"NES"] <- (fgsea_results[,"NES"])
    fgsea_out[,"ES"] <- (fgsea_results[,"ES"])
    fgsea_out[,"Gene_Hits"] <- apply(fgsea_results, 1, function(x) length(unlist(x["leadingEdge"])) )
    fgsea_out[,"Gene_Total"] <- fgsea_results[,"size"]
    output_filename <- file.path(gsea_working_path, paste("fgsea_", analysis_name, ".txt", sep=""))
    write.table(fgsea_out, file=output_filename, sep="\t", row.names=FALSE, quote=FALSE);
    
    topPathways <- fgsea_results[order(pval, decreasing=F), pathway]
    output_filename <- file.path(gsea_images_path, paste("fgsea_", analysis_name, "_allPathways.pdf", sep=""))
    pdf(output_filename, width=4, height=3)
    for(pathway_index in 1:length(topPathways)){ try({
      print(plotEnrichment(pathways_gmt[[topPathways[pathway_index] ]], ranked_vector) +
              labs(title=topPathways[pathway_index]) + theme(plot.title=element_text(size=6)) )
    }, silent=TRUE) }
    dev.off()
    
    ranked_features <- read.delim(rnk_file)
    ranked_vector <- abs(ranked_features[,"rank"])
    names(ranked_vector) <- ranked_features[,"GeneName"]
    pathways_gmt <- gmtPathways(dest_gmt_file)
    suppressWarnings({
      fgsea_results <- fgsea(pathways=pathways_gmt, 
                             stats=ranked_vector,
                             minSize=15, 
                             maxSize=500, 
                             nperm=10000)
    })
    
    fgsea_out <- fgsea_results[,c(1,1,2,3)]
    colnames(fgsea_out) <- c("Term", "Description", "p.Val", "FDR")
    fgsea_out[,"Phenotype"] <- sign(fgsea_results[,"NES"])
    fgsea_out[,"Genes"] <- apply(fgsea_results[,"leadingEdge"], 1, function(x) paste(as.character(unlist(x["leadingEdge"])), collapse=","))
    fgsea_out[,"Genes"] <- apply(fgsea_out, 1, function(x) gsub(" ", "", x["Genes"]))
    fgsea_out[,"NES"] <- (fgsea_results[,"NES"])
    fgsea_out[,"ES"] <- (fgsea_results[,"ES"])
    fgsea_out[,"Gene_Hits"] <- apply(fgsea_results, 1, function(x) length(unlist(x["leadingEdge"])) )
    fgsea_out[,"Gene_Total"] <- fgsea_results[,"size"]
    output_filename <- file.path(gsea_absval_path, paste("fgsea_", analysis_name, "_AbsVal.txt", sep=""))
    write.table(fgsea_out, file=output_filename, sep="\t", row.names=FALSE, quote=FALSE);
    
  }
}

#
##########################