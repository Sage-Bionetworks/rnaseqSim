#! /usr/bin/env Rscript
# KKD for Sage Bionetworks
# Oct. 19, 2016

library(synapseClient)
synapseLogin()


# Generates a new distribution with a hard minimum (100) and gamma shape
makeNewInsertSizeDist=function(targetInsertSize,upper=510,lower=100){
  newdist = rgamma(n = 1000,shape = 5,scale = 0.5)
#  hist(newdist, breaks = 40)
  md = median(newdist)
  
  multFactor = (targetInsertSize+30)/md
  distMin = lower
  distMax = upper
  distSpread = seq(distMin,distMax)
  probs = dgamma(x=distSpread/multFactor,shape = 5,scale = 0.5)
  probsScaled = probs/sum(probs)
  plot(x=distSpread,y=probsScaled)
  total = sum(probsScaled)
  
  write.table(probs, file = "insertString.txt", quote = FALSE, sep = " ",row.names = FALSE,col.names = FALSE,eol = " ")
}


# With peak at 150
# newdist = rgamma(n = 1000,shape = 5,scale = 0.5)
# hist(newdist, breaks = 40)
# md = median(newdist)
# 
# multFactor = 180/md
# distMin = 100
# distMax = 470
# distSpread = seq(distMin,distMax)
# probs = dgamma(x=distSpread/multFactor,shape = 5,scale = 0.5)
# probsScaled = probs/sum(probs)
# plot(x=distSpread,y=probsScaled)
# total = sum(probsScaled)
