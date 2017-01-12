#! /usr/bin/env Rscript
# Manipulates learned isoform TPM values to adjust for dipoid genome simulation.
# KKD for Sage Bionetworks
# Mar 9, 2016

library(argparse)
require(data.table)

parser = ArgumentParser(description='Generate diploid TPM files for read simulation.')
parser$add_argument('--TPM', dest="inpath", type="character", required=TRUE, help='Isoform TPM file from RSEM.')
parser$add_argument('--gtf', type="character", required=TRUE, help='Gene models from which to simulate reads.')
parser$add_argument('--targetDepth', default = 20, type="integer", help='Total million reads to simulate [default %(default)s].')
parser$add_argument('--wd', default = getwd(), type="character", help='Directory for output files [default %(default)s].')
#parser$add_argument('--codeDir', default = getwd(), type="character", help='Directory where R code functions file is located [default %(default)s].')



args = parser$parse_args()
source(Sys.which("isoExpFunctions.R"))
inpath = args$inpath
gtf = read.delim(args$gtf,header = FALSE, comment.char="#")
targetDepth = args$targetDepth
setwd(args$wd)

outName = paste(basename(inpath), 'modDiploid', targetDepth, sep = "_")
outNameFus = paste(basename(inpath), 'modDiploidFusionOnly', targetDepth, sep = "_")

# Read in RSEM model
model = read.delim(inpath)

# look at data
print(paste("sum model TPM", sum(as.numeric(model$TPM))))
sum(model$length)
sum(model$effective_length)
sum(model$expected_count)
sum(model$FPKM)

hist(log(model$TPM), col = "honeydew1", main = "TPM")
hist(log(model$length), col = "honeydew1", main = "all transcripts - length")

print(paste("Fraction not observed", length(which(model$TPM == 0))/ nrow(model)))
print(paste("Fraction not observed", length(which(model$expected_count == 0))/ nrow(model)))
 
# Convert to data.table and remove factors
dtModel <- data.table(model)
dtModel[, gene_id := as.character(gene_id)]
dtModel[, transcript_id := as.character(transcript_id)]

# Calculate total number of transcripts per gene and total tpm per gene
dtModel[, total_tx := .N, by = gene_id]
dtModel[, total_tpm := sum(TPM), by = gene_id] 

# Index the transcripts per gene
dtModel[, index_tx := as.numeric(ave(gene_id, gene_id, FUN=function(x){ sample(length(x)) } )) ]

# Randomly select the number of transcripts to be expressed per gene
#dtModel[, sampled_tx := ceiling(total_tx * runif(1)), by= gene_id]
dtModel[, sampled_tx := rnbinom(1,size=1,prob=.5) + 1, by= gene_id]
dtModel[, sampled_tx := pmin(total_tx,sampled_tx)]

# Randomly select which transcripts will be expressed
dtModel[, tx_exp := runif(1)*(index_tx <= sampled_tx), by = transcript_id]

# Distribute the total gene TPM across the expressed transcripts
dtModel[, new_tpm := tx_exp / sum(tx_exp) * total_tpm, by = gene_id]

# Convert back to data.frame
dtModel <- dtModel[,c("transcript_id","gene_id","length", "effective_length", "expected_count", "new_tpm","FPKM","IsoPct")]
colnames(dtModel) = colnames(model)
model = as.data.frame(dtModel)

# Look at data
hist(log(model$TPM), col = "honeydew1", main = "TPM")

# Modify TPM values to introduce noise to original model.
altModel = addNoise(model)
print(paste('correlation:', cor(log(altModel$TPM), log(model$TPM), method = "spearman")))
print(paste('sum altModel TPM:', sum(altModel$TPM)))


# Get values for the fusion genes
makeFusionMatrix=function(inGTF, inModel){
  head(inGTF)
  geneOnlyGtf = inGTF[inGTF$V3 == "gene",]
  fusions = data.frame(geneOnlyGtf$V1)
  fusions[,2] = geneOnlyGtf$V1
  fusions[,3] = geneOnlyGtf$V5
  fusions[,4] = geneOnlyGtf$V5
  fusions[,5:8] = 0
  colnames(fusions) = colnames(inModel)
  return(fusions)
}
fusions = makeFusionMatrix(inGTF = gtf, inModel = altModel)


# Remove donors from main matrix of isoform values.
fusPartners = t(sapply(as.list(fusions$transcript_id),function(x) {unlist(strsplit(as.character(x),split = '-'))}))
colnames(fusPartners) = c("donor", "acceptor")
toRemove = match(fusPartners[,1], altModel$transcript_id)
altModel = altModel[-toRemove,]



# For sampling fusion expression values, use only the TPM distribution greater than median. 
medianNonZero = median(log(altModel$TPM[altModel$TPM > 0]))
greaterThanMedian = altModel$TPM[log(altModel$TPM) > medianNonZero]
names(greaterThanMedian) = altModel$transcript_id[log(altModel$TPM) > medianNonZero]
fusionExp = data.frame(sampled = sample(greaterThanMedian, size = nrow(fusions), replace = FALSE))
hist(log(fusionExp$sampled), col = "honeydew")


# Adjust fusion TPM and other TPM to sum to 1e6
fractionalAdjustment = 1+(1-((sum(fusionExp$sampled) + sum(altModel$TPM)) / 1e6))
fusionExp$adj = fusionExp$sampled * fractionalAdjustment
altModel$TPM = altModel$TPM * fractionalAdjustment
sum(altModel$TPM) + sum(fusionExp$adj)


# Assign sampled expression values to fusions and originator transcripts.
# Split summed exp value between fusion and donor allele according to splitVal distribution.
splitVal = rnorm(n = nrow(fusPartners),mean = 0.5,sd = 0.03) # selects allelic distribution
fusionExp$fusionAllele = fusionExp$adj*splitVal
fusionExp$donorAllele = fusionExp$adj*(1-splitVal)


# plot TPM vs length of fusions
fusions[,6] = fusionExp$fusionAllele
plot(fusions[,3],log(fusions[,6]), xlab = "fusion length", ylab = "log TPM of fusion")
colnames(fusions) = colnames(model)



# Adjust TPM values between pairs in diploid set.
# *** extract this to a function that can also be used above
diploid = makeDiploid(altModel)
splitVal = rnorm(n = nrow(altModel),mean = 0.5,sd = 0.03)
hist(splitVal)
makeZero = which(altModel$TPM == 0)
splitVal[makeZero] = 0
originalTPM = diploid$TPM[1:length(splitVal)]
diploid$TPM[1:length(splitVal)] = originalTPM*splitVal
diploid$TPM[(length(splitVal)+1):nrow(diploid)] = originalTPM*(1-splitVal)
sum(diploid$TPM)


# Add back in the originator genes
originators = model[toRemove,]
originators$TPM = fusionExp$donorAllele
sum(diploid$TPM, originators$TPM, fusionExp$fusionAllele)
diploidAug = rbind(diploid, originators)


# Plot
hist(log(diploidAug$TPM), col = "honeydew1", main = "TPM")
hist(log(altModel$TPM), col = "honeydew1", main = "TPM")
hist(log(altModel$length), col = "honeydew1", main = "all transcripts - length")
print(paste("Fraction not observed", length(which(diploidAug$TPM == 0))/ nrow(diploidAug)))



# calculate fraction reads for fusions
fusionReads = sum(fusions$TPM)*targetDepth
otherReads = sum(as.numeric(diploidAug$TPM))*targetDepth
print(paste('Number of fusion reads to simulate:', round(fusionReads), sep = ' '))
print(paste('Number of other reads to simulate:', round(otherReads), sep = ' '))
sum(fusionReads, otherReads)/1e6
sum(fusions$TPM, as.numeric(diploidAug$TPM))/1e6

diploidAug$TPM = formatC(x = diploidAug$TPM,digits = 4)

write.table(diploidAug, file = paste(args$wd, outName, sep = "/"), append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE, )
write.table(fusions, file = paste(args$wd, outNameFus, sep = "/"), append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE, )
