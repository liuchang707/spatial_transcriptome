#!/usr/env/python

import re
import sys
import os 
dirr="/storage/peiweikeLab/guochenyu/Liuchang/project2/02.STtools/liver/02.all_liver/02.binsize300_all_C2"
origibarcode=open(dirr+'/collapsedBarcodes.csv')
origigenes=open(dirr+'/collapsedGenes.csv')
normalmatrix=open(dirr+'/collapsedMatrix.mtx')

norb=open('collapsedBarcodes.normal.csv','w')
norgenes=open('collapsedGenes.normal.csv','w')
norm=open('collapsedMatrix.normal.mtx','w')
tdb=open('collapsedBarcodes.td.csv','w')
tdgenes=open('collapsedGenes.td.csv','w')
tdm=open('collapsedMatrix.td.mtx','w')

#### get barcode
normalbar={}
tdbar={}
nbc=0
tbc=0
for line in origibarcode:
    l=line
    line=line.strip().split(',')
    line1=line[0].strip('"')
    if re.search('Colla',line[1]): ##!!! collapse Collase 
        #print ("y")
        if re.search('Collase_tile_1_(.*)_0_0_',line[1]):
            aa=re.search('Collase_tile_1_(.*)_0_0_',line[1])
            #print(aa.group(1))
            if int(aa.group(1))>=2102 and  int(aa.group(1))<=2107:
                nbc+=1
                new='"'+str(nbc)+'"'+','+line[1]+'\n'
                normalbar[line1]=[nbc,new]
            if int(aa.group(1))>=2116 and  int(aa.group(1))<=2119:
                tbc+=1
                new='"'+str(tbc)+'"'+','+line[1]+'\n'
                #print('yes')
                tdbar[line1]=[tbc,new]
#### get genes
genes={}
gc=0
count=0
gsort={}
for line in origigenes:
    count+=1
    if count>1:
        l=line
        line=line.strip().split(',')
        line1=line[0].strip('"')
        line2=line[1].strip('"')
        if not re.search('^Gm(\d+)',line2):
            gc+=1
            genes[line1]=l
            gsort[line1]=gc


#### read matrix and first change barcode count
norg={}
norG=''
norpos={}
norPos=''
nortotal=0
norout=[]
tdg={}
tdG=''
tdpos={}
tdPos=''
tdtotal=0
tdout=[]
c=0

for line in normalmatrix:
    c+=1
    if c>=3:
        line=line.strip().split()
        if line[0] in genes:
            if line[1] in normalbar:
                norout.append(line[0]+' '+str(normalbar[line[1]][0])+' '+line[2]+'\n')
                if line[0] not in norg:
                    norg[line[0]]=1
                if line[1] not in norpos:
                    norPos+=normalbar[line[1]][1]
                    norpos[line[1]]=1
                    print(line[1])
                nortotal+=1
            if line[1] in tdbar:
                tdout.append(line[0]+' '+str(tdbar[line[1]][0])+' '+line[2]+'\n')
                if line[0] not in tdg:
                    tdg[line[0]]=1
                if line[1] not in tdpos:
                    tdpos[line[1]]=1
                    tdPos+=tdbar[line[1]][1]
                tdtotal+=1

#####get final genes 
norgenes.write("\"\",\"x\"\n")
norcount=0
nordic={}
tdgenes.write("\"\",\"x\"\n")
tdcount=0
tddic={}
for i in sorted(gsort.items(),key=lambda item:item[1]):
    if i[0] in norg:
        norcount+=1
        norgenes.write('\"'+str(norcount)+'\",'+genes[i[0]].split(',')[1])
        nordic[i[0]]=norcount
    if i[0] in tdg:
        tdcount+=1
        tdgenes.write('\"'+str(tdcount)+'\",'+genes[i[0]].split(',')[1])
        tddic[i[0]]=tdcount

#####get final matrix
norm.write('%%MatrixMarket matrix coordinate integer general\n'+str(len(norg))+' '+str(len(norpos))+' '+str(nortotal)+'\n')
for key in norout:
    k=key.strip().split()
    newk=str(nordic[k[0]])+' '+k[1]+' '+k[2]+'\n'
    norm.write(newk)

tdm.write('%%MatrixMarket matrix coordinate integer general\n'+str(len(tdg))+' '+str(len(tdpos))+' '+str(tdtotal)+'\n')
for key in tdout:
    k=key.strip().split()
    newk=str(tddic[k[0]])+' '+k[1]+' '+k[2]+'\n'
    tdm.write(newk)



norb.write("\"\",\"x\"\n"+norPos)

tdb.write("\"\",\"x\"\n"+tdPos)

norm.close()
tdm.close()
norb.close()
tdb.close()
norgenes.close()
tdgenes.close()
origibarcode.close()
origigenes.close()
normalmatrix.close()

