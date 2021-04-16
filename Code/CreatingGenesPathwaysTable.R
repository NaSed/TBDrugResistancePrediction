rm(list=ls())

homePath <- '/Users/nafiseh/Google Drive/DrugResistance/Data/KEGG/'

require(KEGGREST)

genes_pathways <- keggLink("pathway", "mtu")

# There is a difference between these and the online KEGG pathways
# For example here Rv11270c is included in the pathway mtu00010 while on the webpage it is not included. So, we extract name of genes in another way
# genes <- names(genes_pathways) 
# genes <- gsub('mtu:','', genes)

pathways <- unname(genes_pathways)

pathways <- unique(gsub('path:', '', pathways))
genes_pathways <- DataFrame()
for (i in pathways){
  tmp <- keggGet(i)
  all <- tmp[[1]]$GENE
  if (!isEmpty(all)){
  Locus <- all[seq(1,length(all),by=2)]
  name <- sub(";.*", "", tmp[[1]]$GENE[seq(2,length(all), by=2)])
  ind <- which(!grepl(";", tmp[[1]]$GENE[seq(2,length(all), by=2)], fixed = TRUE))
  name[ind] <- Locus[ind]
  genes_pathways <- rbind(genes_pathways,cbind(Name = name, Locus=Locus, Pathway = rep(i, length(all))))
  }
}


write.csv(genes_pathways, file = paste(homePath, 'genes_pathways.csv', sep=''), row.names = FALSE)



