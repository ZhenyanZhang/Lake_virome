import os
import sys
from Bio import SeqIO
import re


raw_fasta = sys.argv[1]
eggnog_file = sys.argv[2]
checkv_dir = sys.argv[3]
outfile = sys.argv[4]

sample = os.path.basename(raw_fasta).split(".")[0]

d={}
d2={}
if os.path.exists(eggnog_file):
    f=open(eggnog_file)
    lines=f.readlines()
    for line in lines:
        if not line.strip():
            continue
        if line.startswith("#"):
            continue
        data = line.strip().split('\t')
        query_id = "_".join(data[0].split('_')[:-1])
        ogs = data[4]
        description = data[7].lower()
        d2[query_id]="no"
        if query_id not in d:
            d[query_id] = []
        if "Viruses" in ogs:
            flag="yes"
        elif description=='-':
            flag = 'yes'
        elif 'hypothetical' in description:
            flag='yes'
        elif 'phage' in description or "virus" in description:
            flag='yes'
        elif 'unknow' in description:
            flag='yes'
        else:
            flag='no'
        for key_word in "capsid,phage,terminase,base plate,baseplate,prohead,virion,virus,viral,tape measure,tapemeasure neck,tail,head,bacteriophage,prophage,portal,DNA packaging,T4,p22,holin".split(","):
            if key_word.lower() in description:
                flag='yes'
                d2[query_id]="yes"
                break
        d[query_id].append(flag)
    f.close()
else:
    print(eggnog_file)


file1 = os.path.join(checkv_dir,"proviruses.fna")
file2 = os.path.join(checkv_dir,"viruses.fna")
file3 = os.path.join(checkv_dir,"quality_summary.tsv")

keep_dict = {}
for k,v in d.items():
    #contain key gene words
    if k in d2:
        keep_dict[k]=1
        continue
    #70% protein virus or unknown
    all_flag = len(v)
    yes_flag = v.count('yes')
    name = '_'.join(k.split("_")[:-1])
    if yes_flag/all_flag>=0.7:
        keep_dict[k]=1

tmp = {}
index=1
w=open(outfile,'w+')
#proviruses
if os.path.exists(file1):
    for rec in SeqIO.parse(file1,'fasta'):
        name = '_'.join(rec.id.split("_")[:-1])
        if name not in keep_dict:
            continue
        new_name = sample+'_'+str(index)
        w.write(">"+new_name+' '+name+"\n"+str(rec.seq)+"\n")
        tmp[name]=1
        index+=1
if os.path.exists(file2):
    for rec in SeqIO.parse(file2,'fasta'):
        name = rec.id
        if name not in keep_dict:
            continue
        new_name = sample+'_'+str(index)
        w.write(">"+new_name+' '+name+"\n"+str(rec.seq)+"\n")
        tmp[name]=1
        index+=1
for rec in SeqIO.parse(raw_fasta,'fasta'):
    name = rec.id
    if name in keep_dict and name not in tmp:
        new_name = sample+'_'+str(index)
        w.write(">"+new_name+' '+name+name+"\n"+str(rec.seq)+"\n")
        index+=1
w.close()

