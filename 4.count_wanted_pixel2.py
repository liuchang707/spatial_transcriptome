#!/usr/env/python

import re 
import sys
import os

### open UMI gene count file; wanted region file
umiinf=open(sys.argv[1])
geneinf=open(sys.argv[2])
coorinf=open(sys.argv[3])
### open outfile
umiout=open(sys.argv[4],'w')
geneout=open(sys.argv[5],'w')
pdout=open(sys.argv[6],'w')
umiout.write('sample\tspatialcoordinates\tUMIs\n')
geneout.write('sample\tspatialcoordinates\tgenes\n')
### make coor dic
coor={}
for line in coorinf:
    if not re.search('lane',line):
        line=line.strip().split()
        coor.update({line[0]:[int(line[1]),int(line[2]),int(line[3]),int(line[4])]})

### count 30pixel2 umi 
umi={}
pixel={}
for line in umiinf:
    if not re.search('index',line):
        line=line.strip().split()
        a=line[1].split('_')
        lt=a[0]+'_'+a[1]
        if lt in coor:
            if int(a[2])>=coor[lt][0] and int(a[2])<=coor[lt][1] and int(a[3])>=coor[lt][2] and int(a[3])<=coor[lt][3]:
                if int(line[2])<=100:
                    #umiout.write('liver\t'+line[1]+'\t'+line[2]+'\n')
                    aa=(int(a[2])-coor[lt][0])//30 ### x ,30 pixel= 1 weimi
                    bb=(int(a[3])-coor[lt][2])//30 ### y, 30....
                    cc=str(aa)+'_'+str(bb)
                    if lt not in umi:
                        umi.update({lt:{cc:int(line[2])}})
                        pixel.update({lt:{cc:1}})
                    elif cc not in umi[lt]:
                        umi[lt].update({cc:int(line[2])})
                        pixel[lt].update({cc:1})
                    else:
                        umi[lt][cc]+=int(line[2])
                        pixel[lt][cc]+=1


for key1 in umi:  
    for key2 in umi[key1]:
        #if pixel[key1][key2]>=800: ## filter last x y pixel2 
        pdout.write('liver\t'+key1+'\t'+key2+'\t'+str(pixel[key1][key2])+'\n')
        umiout.write('liver\t'+key1+'\t'+key2+'\t'+str(umi[key1][key2])+'\n')

### 
gene={}
for line in geneinf:
    if not re.search('index',line):
        line=line.strip().split()
        a=line[1].split('_')
        lt=a[0]+'_'+a[1]
        if lt in coor:
            if int(a[2])>=coor[lt][0] and int(a[2])<=coor[lt][1] and int(a[3])>=coor[lt][2] and int(a[3])<=coor[lt][3]:
                if int(line[2])<=100:
                    #geneout.write('liver\t'+line[1]+'\t'+line[2]+'\n')
                    aa=(int(a[2])-coor[lt][0])//30
                    bb=(int(a[3])-coor[lt][2])//30
                    cc=str(aa)+'_'+str(bb)
                    if lt not in gene:
                        gene.update({lt:{cc:int(line[2])}})
                    elif cc not in gene[lt]:
                        gene[lt].update({cc:int(line[2])})
                    else:
                        gene[lt][cc]+=int(line[2])

for key1 in gene:
    for key2 in gene[key1]:
        geneout.write('liver\t'+key1+'\t'+key2+'\t'+str(gene[key1][key2])+'\n')



umiinf.close()
geneinf.close()
coorinf.close()
umiout.close()
geneout.close()
pdout.close()
