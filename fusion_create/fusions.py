import random
#import time
import gffutils
from gffutils import biopython_integration

MAX_TOSS_NUM = 100


def getTranscript(db, gene):
    trans = []
    result =  db.children(gene,featuretype = "transcript")
    for item in result:
        if item.source == 'protein_coding':
            trans.append(item)
    if (len(trans) == 0):
    	return None,None
    else:
    	tId = random.randint(0,len(trans)-1)
    return trans[tId].strand,trans[tId]['transcript_id'][0]

    
def getJunctionAtExonBoundary(db, tranId, strand, isDonor):
# TODO: figure out why this is sometimes returning empty lists    
    exons = []
    result =  db.children(tranId, featuretype='exon', order_by='start')
    for item in result:
        exons.append(item)
    if len(exons)-2 < 0:
        return False,999,999
    if (strand == '+' and isDonor) or (strand == '-' and (not isDonor)):    #Kristen: please check all the genes has '+' or '-' not '.'
        eId = random.randint(0,len(exons)-2)
        fusExons = exons[0:(eId+1)] 
        return True,exons[eId].end,fusExons
    elif (strand == '-' and isDonor) or (strand == '+' and (not isDonor)):
        eId = random.randint(1,len(exons)-1)
        fusExons = exons[eId:]
        return True,exons[eId].start,fusExons     #gffutils is 1-based
    else:
        return False,999,999


def isStay(pStay):
    # the function has pStay probability return True
    p = random.random()
    if p < pStay:
        return True
    else:
        return False


def getRandomFusions(db, names, num=5, pStay=0.0):
    # db is the database from module.py 
    # names is a vector of the ENSG gene ids (protein coding genes only) from module.py
    # num: number of fusions to simulate
    # pStay: the probability of staying in the same gene pair to generate another fusion isoform. Set to 0.0 to get only one isoform. 
    #random.seed(time.time())    #for final code, add parameters for seed
    
    res = list() # the list to store the dictionaries for fusions
                 # donorTranId, acceptorTranId donorJunction, acceptorJunction

    if len(names) <2:
        print("Not enough protein coding genes.")
        exit(1)

    total=0
    tossed=0
    while total < num:
        # Select genes.
        dId = random.randint(0,len(names)-1)
        aId = random.randint(0,len(names)-1)

        # Discard the result if the genes selected are the same.
        if dId==aId:
           tossed = tossed + 1
           if tossed > MAX_TOSS_NUM:
              print("Tossed > "+str(MAX_TOSS_NUM)+" times in generating a pair of genes.")
              exit(1)
           continue
            
        dGene = db[names[dId]]
        aGene = db[names[aId]]
        
        # Decide whether to keep the same transcript pair for the next fusion event.
        tossed2 = 0
        keepSame = True
        while keepSame is True:
            keepSame = isStay(pStay)
            # Choose transcripts
            dStrand,dTran = getTranscript(db, dGene)
            aStrand,aTran = getTranscript(db, aGene)
            if (dTran is None) or (aTran is None): continue
            # Choose junctions
            dIsSucess,dJunction,dExons = getJunctionAtExonBoundary(db, dTran, dStrand, True)
            aIsSucess,aJunction,aExons = getJunctionAtExonBoundary(db, aTran, aStrand, False)             
            if dIsSucess and aIsSucess:
                dExonSF = list()
                aExonSF = list()
                for exon in dExons:
                   dExonSF.append(biopython_integration.to_seqfeature(exon))
                for exon in aExons:
                   aExonSF.append(biopython_integration.to_seqfeature(exon))
                if (len(dExonSF) > 0) and (len(aExonSF) > 0):
                	# adjust junction positions to be 0-based
                   res.append({'donorExons':dExonSF,'acceptorExons':aExonSF,'dJunction':dJunction-1,'aJunction':aJunction-1})
                
                   total = total + 1
            else:            
               tossed2 = tossed2  + 1
               if tossed2 > MAX_TOSS_NUM:
                  print("Tossed > "+str(MAX_TOSS_NUM)+" times in generating fusion junctions.")
                  exit(1)
    return(res)
     
