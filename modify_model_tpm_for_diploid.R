#! /usr/bin/env Rscript
# Manipulates learned isoform TPM values to adjust for dipoid genome simulation.
# KKD for Sage Bionetworks
# Mar 9, 2016

library(argparse)

parser = ArgumentParser(description='Generate diploid TPM files for read simulation.')
parser$add_argument('--model', dest="inpath", type="character", required=TRUE, help='RSEM model file.')
parser$add_argument('--gtf', type="character", required=TRUE, help='Gene models from which to simulate reads.')
parser$add_argument('--expThresh', default = 2, type="double", help='Minimum log TPM value for fusion gene [default %(default)s].')
parser$add_argument('--targetDepth', default = 20, type="integer", help='Total million reads to simulate [default %(default)s].')
parser$add_argument('--wd', default = getwd(), type="character", help='Directory for output files [default %(default)s].')


args = parser$parse_args()

setwd(get(args$wd))
inpath = get(args$accumulate)
gtf = read.delim(get(args$gtf),header = FALSE)
expThreshold = get(args$expThresh)
targetDepth = get(args$targetDepth)

outpath = paste(inpath, 'modDiploid', expThreshold, targetDepth, sep = "_")
outpathFus = paste(inpath, 'modDiploidFusionOnly', expThreshold, targetDepth, sep = "_")

# Read in RSEM model
model = read.delim(inpath)
head(model)

# look at data
sum(model$TPM)
sum(model$length)
sum(model$effective_length)
sum(model$expected_count)
sum(model$FPKM)

hist(log(model$TPM), col = "honeydew1", main = "TPM")
hist(log(model$length), col = "honeydew1", main = "all transcripts - length")

print(paste("Fraction not observed", length(which(model$TPM == 0))/ nrow(model)))
print(paste("Fraction not observed", length(which(model$expected_count == 0))/ nrow(model)))
hist(log(model$length[which(model$TPM > 0)]), col = "honeydew1", main = "observed transcripts - length")


# create diploid set
temp = model
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

hist(log(diploid$TPM), col = "honeydew1", main = "TPM")
hist(log(model$TPM), col = "honeydew1", main = "TPM")
hist(log(model$length), col = "honeydew1", main = "all transcripts - length")

print(paste("Fraction not observed", length(which(diploid$TPM == 0))/ nrow(diploid)))
hist(log(diploid$length[which(diploid$TPM > 0)]), col = "honeydew1", main = "observed transcripts - length")


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
print(paste('Number of fusion reads to simulate:', fusionReads, sep = ' '))
print(paste('Number of other reads to simulate:', otherReads, sep = ' '))
sum(fusionReads, otherReads)/1e6
sum(fusions$TPM, as.numeric(diploid$TPM))/1e6


write.table(diploid, file = outpath, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE, )
write.table(fusions, file = outpathFus, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE, )
