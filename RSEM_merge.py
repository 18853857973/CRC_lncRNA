#!/usr/bin/python3


#input two paramater
#col column
#output file name
import os
import sys
mydict = {}
allfile=os.listdir(os.getcwd())
myfile=[]
col=sys.argv[1]
outputfile=sys.argv[2]

fprint=open(outputfile+'.txt','w')

for file in allfile:
    if file.endswith('lncRNA.genes.results'):
        myfile.append(file)
        print (file)
        for line in open(file, 'r'):
            key,value = line.strip().split('\t')[0],line.strip().split('\t')[int(col)-1]
            if key in mydict:
                mydict[key] = mydict[key] + '\t' + value
            else:
                mydict[key] = value

for mf in myfile:
    fprint.write(mf.strip('.genes.results')+'\t')
fprint.write('\r')

for key,value in mydict.items():
    fprint.write(key + '\t' + value +'\r')
fprint.close()

print ('you chose the '+col+' column')

