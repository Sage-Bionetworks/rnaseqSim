# Manipulates learned isoform TPM values to adjust for dipoid genome simulation.
# KKD for Sage Bionetworks
# Mar 9, 2016

setwd('~/Computing/cancer/SMC_RNA/OICR_samples/')
expThreshold = 1 # in log TPM scale
targetDepth = 100 # in million reads

# star
inpath = "~/Computing/cancer/SMC_RNA/OICR_samples/star-ref/CPCG_0340.isoforms.results"
outpath = paste(inpath, 'modDiploid', expThreshold, targetDepth, sep = "_")
outpathFus = paste(inpath, 'modDiploidFusionOnly', expThreshold, targetDepth, sep = "_")

# drop suffix numbers from transcript and gene ids -- not necessary for star refs
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
gtf = read.delim("~/Computing/cancer/SMC_RNA/sim1/unfiltered_sim1_filtered.gtf",header = FALSE)
head(gtf)
fusions = data.frame(gtf$V1)
fusions[,2] = gtf$V1
fusions[,3] = gtf$V5
fusions[,4] = gtf$V5
fusions[,5:8] = 0
head(fusions)


# use only the TPM distribution greater than 0 (in log) for sampling
greaterThanLogZero = model$TPM[log(model$TPM) > expThreshold]
fusionExp = sample(greaterThanLogZero, size = nrow(gtf), replace = FALSE)
hist(log(fusionExp), col = "honeydew")


fusions[,6] = fusionExp
head(fusions)
toZero = match(fusionExp, model$TPM)
diploid$TPM[toZero] = 0
diploid$TPM[(nrow(model)+toZero)] = 0


# plot TPM vs length of fusions
plot(fusions[,3],log(fusions[,6]), xlab = "fusion length", ylab = "log TPM of fusion")
colnames(fusions) = colnames(diploid)


# alter TPM of fusions file based on target read depth
targetTPM = fusions$TPM * targetDepth # remember that targetDepth is in million reads
fusions$TPM = targetTPM
fusions

write.table(diploid, file = outpath, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE, )
write.table(fusions, file = outpathFus, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE, )
