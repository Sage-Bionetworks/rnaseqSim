import sys
import os
import gffutils
import random
import time

MAX_TOSS_NUM = 100


def getTranscript(db, gene):
    trans = []
    transcripts = db.children(gene,featuretype = "transcript")
    for item in transcripts:
        trans.append(item)
    tId = random.randint(0,len(trans)-1)
    return trans[tId].strand,trans[tId]['transcript_id'][0]
    
def getJunctionAtExonBoundary(db, gene, tranId, strand, isDonor):
    exons = []
    for item in db.children(gene, featuretype='exon', order_by='start'):
        if item['transcript_id'][0] == tranId:
            exons.append(item)
    if len(exons)-2 < 0:
        return False,0
    if (strand == '+' and isDonor) or (strand == '-' and (not isDonor)):    #Kristen: please check all the genes has '+' or '-' not '.'
        eId = random.randint(0,len(exons)-2)
        return True,exons[eId].end
    else:
        eId = random.randint(1,len(exons)-1)
        return True,exons[eId].start     #Kristen: please check it is 1 based or 0 based

def isStay(pStay):
    p = random.random()
    if p < pStay:
        return True
    else:
        return False

def run_fusions():
    
    random.seed(time.time())    #for real code, add parameters for seed
    
    db = gffutils.FeatureDB('tmp2.sqlite3')    #for real code, use def run_fusions(db,names,num,pStay):
    names = []
    names.append('ENSG00000272647')
    names.append('ENSG00000106261')
    names.append('ENSG00000146826')
    names.append('ENSG00000146834')
    names.append('ENSG00000087077')
    num=10
    pStay=0.0

    #real code

    res = list()

    if len(names) <2:
        print "Not enough protein coding genes."
        exit(1)

    total=0
    tossed=0
    while True:
        tossed = tossed + 1
        if tossed > MAX_TOSS_NUM:
            print "Tossed > "+str(MAX_TOSS_NUM)+" times in generating a pair of genes."
            exit(1)
        #select 
        dId = random.randint(0,len(names)-1)
        aId = random.randint(0,len(names)-1)
        print "dId,aId=",dId,aId
        if dId==aId:
            continue
            
        #choose transcripts
        dGene = db[names[dId]]
        aGene = db[names[aId]]
        tossed2 = 0
        while True:
            tossed2 = tossed2  + 1
            if tossed2 > MAX_TOSS_NUM:
                print "Tossed > "+str(MAX_TOSS_NUM)+" times in generating fusion junctions."
                exit(1)
            dStrand,dTran = getTranscript(db, dGene)
            aStrand,aTran = getTranscript(db, aGene)
            dIsSucess,dJunction = getJunctionAtExonBoundary(db, dGene, dTran, dStrand, True)
            aIsSucess,aJunction = getJunctionAtExonBoundary(db, aGene, aTran, aStrand, False)             
            if dIsSucess and aIsSucess:
                res.append({'donorTranId':dTran, 'acceptorTranId':aTran, 'donorJunction':dJunction, 'acceptorJunction':aJunction })
                print dTran,aTran,dJunction,aJunction    #for test purpose
                if isStay(pStay):
                    continue
                else:
                    total = total + 1
                    print total
                    if total >= num:
                        return res
                    else:
                        break
            else:
                continue 
     

if __name__ == '__main__':
    run_fusions()
