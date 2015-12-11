library(biomaRt)
library(Gviz)
Hs=useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org") # use this one when biomart.org is down
Hs = useDataset("hsapiens_gene_ensembl", Hs)


validation = read.csv('~/Computing/cancer/SMC_RNA/TCGA_LAML_data/TCGA_AML_fusion_verification_matrix.csv')
head(validation)
allHCNGs = unlist(strsplit(as.character(validation$TCGA.ID),split = "-"))
fusionPairs = matrix(allHCNGs,nrow = length(allHCNGs)/2,ncol = 2, byrow = TRUE)
head(fusionPairs)
fusionPairsENSG = matrix(NA,nrow = nrow(fusionPairs),ncol = 10)


x = getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "chromosome_name", "start_position", "end_position"),filters="hgnc_symbol",values=toupper(fusionPairs[,1]), mart=Hs)    
fusionPairsENSG[,1] = toupper(fusionPairs[,1])
fusionPairsENSG = as.data.frame(fusionPairsENSG)
fusionPairsENSG[,2:5] = x[match(fusionPairsENSG[,1],x$hgnc_symbol),2:5]


x = getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "chromosome_name", "start_position", "end_position"),filters="hgnc_symbol",values=toupper(fusionPairs[,2]), mart=Hs)    
fusionPairsENSG[,6] = toupper(fusionPairs[,2])
fusionPairsENSG[,7:10] = x[match(fusionPairsENSG[,6],x$hgnc_symbol),2:5]


completePairs = na.omit(fusionPairsENSG)


rm(x)


idKey = read.delim('~/Computing/cancer/SMC_RNA/TCGA_LAML_data/TCGA_LAML_RNAseq_data_summary.tsv')
idKey$mainID = substr(idKey$barcode,start = 9, stop = 12)


idKey$mainID[which(as.character(idKey$analysis_id) == "05911238-a724-4748-896d-33f6a5f8438d")]
idKey$analysis_id[which(as.character(idKey$mainID) %in% substr(colnames(validation),start = 2,stop=5)[5:84])]





head(validation)

shortValidation = validation[which(toupper(as.character(validation$TCGA.ID)) %in% paste(completePairs[,1],completePairs[,6],sep="-")),]

for (i in 1:nrow(shortValidation)){
  sampsToCheck = colnames(shortValidation)[which(shortValidation[i,] %in% c(1,2))]
  idsToCheck = idKey$analysis_id[which(as.character(idKey$mainID) %in% substr(sampsToCheck,start = 2,stop=5))]
  print(as.character(idsToCheck))
  posInCompletePairs = match(toupper(as.character(shortValidation$TCGA.ID)[i]), paste(completePairs[,1],completePairs[,6],sep="-"))
  
  for (j in 1:idsToCheck) {
    bgFile = paste("", idsToCheck[j])
    dTrack = DataTrack(range = bgFile, genome = "hg38",type = "l", chromosome = as.character(completePairs[posInCompletePairs,3]), start = as.character(completePairs[posInCompletePairs,4]), end = as.character(completePairs[posInCompletePairs,5]), name = "bedGraph")
    plotTracks(dTrack, main = shortValidation$TCGA.ID[i])
  }
}

#bgFile = "~/Computing/cancer/SMC_RNA/TCGA_LAML_data/05911238-a724-4748-896d-33f6a5f8438dAligned.sortedByCoord.out.bedgraph"
#dTrack2 <- DataTrack(range = bgFile, genome = "hg38",type = "l", chromosome = "19", name = "bedGraph")
#class(dTrack2)
#plotTracks(dTrack2)
