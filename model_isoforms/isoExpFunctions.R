
# Noise is modeled by a gamma distribution. 
# Noise is either added or subtracted from original value accoring to binomial variable.
addNoise=function(inModel){
  dirChange = rbinom(n = nrow(inModel), size = 1, prob = 0.5)
  fracChange = rgamma(n = nrow(inModel),shape = 2,scale = 1)
  #hist(fracChange)
  altModel = inModel
  posChanges = which(dirChange == 1)
  altModel$TPM[posChanges] = inModel$TPM[posChanges] + (fracChange[posChanges]*inModel$TPM[posChanges])
  negChanges = which(dirChange == 0)
  altModel$TPM[negChanges] = inModel$TPM[negChanges] - (fracChange[negChanges]*inModel$TPM[negChanges])
  altModel$TPM[which(altModel$TPM < 0)] = 0
  preSum = sum(altModel$TPM)
  altModel$TPM = altModel$TPM/preSum * 1e6
  #plot(log(altModel$TPM), log(model$TPM))
  return(altModel)
}


# create diploid set
makeDiploid=function(inModel){
  # Assumes input data frame has columns "transcript_id" and "gene_id" and "TPM"
  diploid = rbind(inModel[,1:8], inModel[,1:8])
  diploid$transcript_id = as.character(diploid$transcript_id)
  diploid$gene_id = as.character(diploid$gene_id)
  #  head(diploid)
  diploid$transcript_id[1:nrow(inModel)] = paste(inModel$transcript_id, "hap1", sep = "-")
  diploid$transcript_id[(nrow(inModel)+1):nrow(diploid)] = paste(inModel$transcript_id, "hap2", sep = "-")
  diploid$gene_id[1:nrow(inModel)] = paste(inModel$gene_id, "hap1", sep = "-")
  diploid$gene_id[(nrow(inModel)+1):nrow(diploid)] = paste(inModel$gene_id, "hap2", sep = "-")
  #  head(diploid)
  #  tail(diploid)
  #  diploid[(nrow(inModel)-5):(nrow(inModel)+5),]
  print(paste('sum diploid TPM:', sum(diploid$TPM)))
  return(diploid)
}
