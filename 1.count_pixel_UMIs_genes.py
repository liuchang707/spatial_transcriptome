#!/usr/env/python

import sys
import os
import gzip

matrix=gzip.open(sys.argv[1],"rt")
coord=gzip.open(sys.argv[2],"rt")
barcodes=gzip.open(sys.argv[3],"rt")
features=gzip.open(sys.argv[4],"rt")

out1=open(sys.argv[5],'w')
out2=open(sys.argv[6],'w')
dic={}
### read matrix 
lcount=0
for line in matrix:
    lcount+=1
    if lcount>3:
        line=line.strip().split(' ')
        if line[1] not in dic:
            dic.update({line[1]:{line[0]:int(line[2])}})
        else:
            dic[line[1]].update({line[0]:int(line[2])})
### read coord
coor={}
for line in coord:
    line=line.strip().split()
    a=line[1]+'_'+line[2]+'_'+line[3]+'_'+line[4]
    coor[line[0]]=a

### read barcodes
bar=[]
for line in barcodes:
    line=line.strip().split()
    bar.append(line[0])
#print(bar[0])

### read features
feat=[]
for line in features:
    line=line.strip().split()
    feat.append(line[0])

out1.write('index\tspatialcoordinates\tUMIs\n')
out2.write('index\tspatialcoordinates\tgenes\n')

for key1 in dic:
    umi=0
    gene=0
    g=''
    for key2 in dic[key1]:
        umi+=dic[key1][key2]
        gene+=1
        g+=feat[int(key2)-1]+';'
    index=int(key1)-1
    #print(bar[index])
    if bar[index] in coor:
        out1.write(bar[index]+'\t'+coor[bar[index]]+'\t'+str(umi)+'\n')
        out2.write(bar[index]+'\t'+coor[bar[index]]+'\t'+str(gene)+'\t'+g+'\n')
    else:
        print(bar[index])

matrix.close()
coord.close()
barcodes.close()
features.close()
out1.close()
out2.close()



