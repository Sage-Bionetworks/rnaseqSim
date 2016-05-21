#! /usr/bin/env Rscript
# Manipulates learned isoform TPM values to adjust for dipoid genome simulation.
# KKD for Sage Bionetworks
# Mar 9, 2016

library(argparse)

parser = ArgumentParser(description='Generate diploid TPM files for read simulation.')
parser$add_argument('--TPM', dest="inpath", type="character", required=TRUE, help='Isoform TPM file from RSEM.')
parser$add_argument('--gtf', type="character", required=TRUE, help='Gene models from which to simulate reads.')
parser$add_argument('--expThresh', default = 2, type="double", help='Minimum log TPM value for fusion gene [default %(default)s].')
parser$add_argument('--targetDepth', default = 20, type="integer", help='Total million reads to simulate [default %(default)s].')
parser$add_argument('--wd', default = getwd(), type="character", help='Directory for output files [default %(default)s].')


args = parser$parse_args()

inpath = args$inpath
gtf = read.delim(args$gtf,header = FALSE)
expThreshold = args$expThresh
targetDepth = args$targetDepth
setwd(args$wd)

outName = paste(basename(inpath), 'modDiploid', expThreshold, targetDepth, sep = "_")
outNameFus = paste(basename(inpath), 'modDiploidFusionOnly', expThreshold, targetDepth, sep = "_")

# Read in RSEM model
model = read.delim(inpath)

# look at data
print(paste("sum model TPM", sum(as.numeric(model$TPM))))
sum(model$length)
sum(model$effective_length)
sum(model$expected_count)
sum(model$FPKM)

#hist(log(model$TPM), col = "honeydew1", main = "TPM")
#hist(log(model$length), col = "honeydew1", main = "all transcripts - length")

print(paste("Fraction not observed", length(which(model$TPM == 0))/ nrow(model)))
print(paste("Fraction not observed", length(which(model$expected_count == 0))/ nrow(model)))
#hist(log(model$length[which(model$TPM > 0)]), col = "honeydew1", main = "observed transcripts - length")


# Modify TPM values to introduce noise to original model.
# Noise is modeled by a gamma distribution.
dirChange = rbinom(n = nrow(model), size = 1, prob = 0.5)
fracChange = rgamma(n = nrow(model),shape = 2,scale = 1)
#hist(fracChange)
altModel = model
posChanges = which(dirChange == 1)
altModel$TPM[posChanges] = model$TPM[posChanges] + (fracChange[posChanges]*model$TPM[posChanges])
negChanges = which(dirChange == 0)
altModel$TPM[negChanges] = model$TPM[negChanges] - (fracChange[negChanges]*model$TPM[negChanges])
altModel$TPM[which(altModel$TPM < 0)] = 0
preSum = sum(altModel$TPM)
altModel$TPM = altModel$TPM/preSum * 1e6
#plot(log(altModel$TPM), log(model$TPM))
print(paste('correlation:', cor(log(altModel$TPM), log(model$TPM), method = "spearman")))
print(paste('sum altModel TPM:', sum(altModel$TPM)))

# create diploid set
temp = altModel
diploid = rbind(temp[,1:8], temp[,1:8])
diploid$transcript_id = as.character(diploid$transcript_id)
diploid$gene_id = as.character(diploid$gene_id)
head(diploid)
diploid$transcript_id[1:nrow(model)] = paste(temp$transcript_id, "hap1", sep = "-")
diploid$transcript_id[(nrow(model)+1):nrow(diploid)] = paste(temp$transcript_id, "hap2", sep = "-")
diploid$gene_id[1:nrow(model)] = paste(temp$gene_id, "hap1", sep = "-")
diploid$gene_id[(nrow(model)+1):nrow(diploid)] = paste(temp$gene_id, "hap2", sep = "-")
head(diploid)
tail(diploid)
diploid[(nrow(model)-5):(nrow(model)+5),]
sum(diploid$TPM)


# Adjust TPM values between pairs in diploid set.
# Assume expression is allocated between alleles according to a value drawn from a normal distribution centered on 0.5.
splitVal = rnorm(n = nrow(model),mean = 0.5,sd = 0.03)
hist(splitVal)
makeZero = which(model$TPM == 0)
splitVal[makeZero] = 0
originalTPM = diploid$TPM[1:length(splitVal)]
diploid$TPM[1:length(splitVal)] = originalTPM*splitVal
diploid$TPM[(length(splitVal)+1):nrow(diploid)] = originalTPM*(1-splitVal)
sum(diploid$TPM)

#hist(log(diploid$TPM), col = "honeydew1", main = "TPM")
#hist(log(model$TPM), col = "honeydew1", main = "TPM")
#hist(log(model$length), col = "honeydew1", main = "all transcripts - length")

print(paste("Fraction not observed", length(which(diploid$TPM == 0))/ nrow(diploid)))
#hist(log(diploid$length[which(diploid$TPM > 0)]), col = "honeydew1", main = "observed transcripts - length")


head(diploid)
diploid$TPM = formatC(x = diploid$TPM,digits = 5)
head(diploid)




# get values for the fusion genes
head(gtf)
geneOnlyGtf = gtf[gtf$V3 == "gene",]
fusions = data.frame(geneOnlyGtf$V1)
fusions[,2] = geneOnlyGtf$V1
fusions[,3] = geneOnlyGtf$V5
fusions[,4] = geneOnlyGtf$V5
fusions[,5:8] = 0
head(fusions)


# use only the TPM distribution greater than expThreshold (in log) for sampling
greaterThanLogZero = model$TPM[log(model$TPM) > expThreshold]
fusionExp = sample(greaterThanLogZero, size = nrow(geneOnlyGtf), replace = FALSE)
hist(log(fusionExp), col = "honeydew")


fusions[,6] = fusionExp
head(fusions)
toZero = match(fusionExp, model$TPM)
diploid$TPM[toZero] = 0
diploid$TPM[(nrow(model)+toZero)] = 0


# plot TPM vs length of fusions
plot(fusions[,3],log(fusions[,6]), xlab = "fusion length", ylab = "log TPM of fusion")
colnames(fusions) = colnames(diploid)


# calculate fraction reads for fusions
fusionReads = sum(fusions$TPM)*targetDepth
otherReads = sum(as.numeric(diploid$TPM))*targetDepth
print(paste('Number of fusion reads to simulate:', round(fusionReads), sep = ' '))
print(paste('Number of other reads to simulate:', round(otherReads), sep = ' '))
sum(fusionReads, otherReads)/1e6
sum(fusions$TPM, as.numeric(diploid$TPM))/1e6



write.table(diploid, file = paste(args$wd, outName, sep = "/"), append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE, )
write.table(fusions, file = paste(args$wd, outNameFus, sep = "/"), append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE, )