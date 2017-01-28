#! /usr/bin/env Rscript
# KKD for Sage Bionetworks
# Oct. 19, 2016

library(synapseClient)
synapseLogin()
library(argparse)
library(yaml)

parser = ArgumentParser(description='Model file for RSEM FASTQ simulation based on input parameters.')
parser$add_argument('--param', dest="param", type="character", required=TRUE, help='Input parameter file (YAML-format).')
parser$add_argument('--inModel', dest="synid", type="character", required=TRUE, help='Synapse ID of input model.')
parser$add_argument('--out', dest="out", type="character", required = TRUE, help='Name of output model file OR Synapse ID of model file to replace.')


args = parser$parse_args()

#####################
# Functions to modify model parameters
#####################

# Generates a new distribution with a gamma shape. Minimum cannot be less than read length.
makeNewInsertSizeDist=function(targetInsertSize,upper=510,lower=100){
  newdist = rgamma(n = 5000,shape = 5,scale = 0.5)
  md = median(newdist)
  
  multFactor = (targetInsertSize+30)/md 
  distSpread = seq(lower,upper)
  probs = dgamma(x=distSpread/multFactor,shape = 5,scale = 0.5)
  probsScaled = probs/sum(probs)
  plot(x=distSpread,y=probsScaled)
  total = sum(probsScaled)
  
  output = list(spread=paste(lower,upper,length(distSpread),sep = " "), probs=paste(probsScaled, collapse = " "))
  return(output)
}


# Generates a distribution for 3' coverage bias
makeBiasedDist=function(numBins=20){
  newdist = hist(rbeta(1000, 1.7,1.1),breaks = (numBins-1))
  probs = newdist$density/numBins
  output = paste(probs, collapse = " ")
  return(output)
}

params = yaml.load_file(args$param)
modelEnt = synGet(args$synid)
model = read.delim(getFileLocation(modelEnt),header = FALSE,as.is = TRUE)
newModel = data.frame()

#model_type # 0, single-end, no quality score; 1, single-end, quality score; 2, paired-end, no quality score; 3, paired-end, quality score
newModel[1,1] = '3'

# forward probability
if(params$stranded == TRUE) newModel[2,1] = '1' else  newModel[2,1] = '0.5' 

# Set insert size distribution
res = makeNewInsertSizeDist(params$insertSize,lower = params$readLength)
newModel[3,1] = res$spread
newModel[4,1] = res$probs


# Read length distribution
span = 1
newModel[5,1] = paste(params$readLength-span,params$readLength,span,sep = ' ')
newModel[6,1] = '1'

# Set read coverage bias across transcripts
if(params$coverageBias == 'uniform'){ 
  newModel[7,1] = '0' 
  temp = data.frame("V1"=model[10:nrow(model),1])
  finalModel = rbind(newModel,temp)
} else {
  newModel[7,1] = '1'
  newModel[8,1] = '20'
  newModel[9,1] = makeBiasedDist()
  temp = data.frame("V1"=model[10:nrow(model),1])
  finalModel = rbind(newModel,temp)
  
}

annotations = list(coverage_bias=params$coverageBias,insert_size=params$insertSize,original_sample=modelEnt@fileHandle$fileName)

if (substr(args$out,start = 1,stop = 3) == "syn") {
  finalModelEnt = synGet(args$out)
  write.table(finalModel,file = finalModelEnt@filePath, append = FALSE,quote = FALSE,row.names = FALSE,col.names = FALSE)
  synSetAnnotations(finalModelEnt) = annotations
  finalModelEnt = synStore(finalModelEnt)
} else {
  write.table(finalModel,file = args$out,append = FALSE,quote = FALSE,row.names = FALSE,col.names = FALSE)
  finalModelEnt = File(path = args$out, synapseStore = TRUE, parentId="syn7506232")
  synSetAnnotations(finalModelEnt) = annotations
  finalModelEnt = synStore(finalModelEnt)
}

