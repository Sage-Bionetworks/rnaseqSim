import random
import time
import gffutils
from gffutils import biopython_integration

MAX_TOSS_NUM = 100


def getTranscript(db, gene):
    trans = []
    for item in db.children(gene,featuretype = "transcript"):
        if item.source == 'protein_coding':
            trans.append(item)
    if (len(trans) == 0):
    	return None,None
    else:
    	tId = random.randint(0,len(trans)-1)
    return trans[tId].strand,trans[tId]['transcript_id'][0]
    
    
def getJunctionAtExonBoundary(db, tranId, strand, isDonor):
    exons = []
    for item in db.children(tranId, featuretype='exon', order_by='start'):
        exons.append(item)
    if len(exons)-2 < 0:
        return False,999,999
    if (strand == '+' and isDonor) or (strand == '-' and (not isDonor)):    #Kristen: please check all the genes has '+' or '-' not '.'
        eId = random.randint(0,len(exons)-2)
        return True,exons[eId].end,eId
    else:
        eId = random.randint(1,len(exons)-1)
        return True,exons[eId].start,eId     #gffutils is 1-based


def isStay(pStay):
    # the function has pStay probability return True
    p = random.random()
    if p < pStay:
        return True
    else:
        return False


def getRandomFusions(db, names, num=20, pStay=0.0):
    # db is the database from module.py 
    # names is a vector of the ENSG gene ids (protein coding genes only) from module.py
    # num: number of fusions to simulate
    # pStay: the probability of staying in the same gene pair to generate another fusion isoform. Set to 0.0 to get only one isoform. 
    random.seed(time.time())    #for final code, add parameters for seed
    
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
            dIsSucess,dJunction,dEid = getJunctionAtExonBoundary(db, dTran, dStrand, True)
            aIsSucess,aJunction,aEid = getJunctionAtExonBoundary(db, aTran, aStrand, False)             
            if dIsSucess and aIsSucess:
                donor = db[dTran]
                donor['donorJunction'] = dJunction
                donor['junctionExonNum'] = dEid
                acceptor = db[aTran]
                acceptor['acceptorJunction'] = aJunction
                acceptor['junctionExonNum'] = aEid
                res.append({'donor':donor, 'acceptor':acceptor })
                #print dTran,aTran,dJunction,aJunction    #for test purpose
                
                total = total + 1
            else:            
               tossed2 = tossed2  + 1
               if tossed2 > MAX_TOSS_NUM:
                  print("Tossed > "+str(MAX_TOSS_NUM)+" times in generating fusion junctions.")
                  exit(1)
    return(res)
     

def convertToSeqObj(fusionDict):
   """Converts results into Seq objects."""
   
   donorSeq = biopython_integration.to_seqfeature(fusionDict['donor'])
   acceptorSeq = biopython_integration.to_seqfeature(fusionDict['acceptor'])
   return donorSeq,acceptorSeq
   