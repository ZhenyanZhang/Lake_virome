import os
import sys
from Bio import SeqIO

protein_file = sys.argv[1]
db_file = sys.argv[2]
db_ann = sys.argv[3]
outdir = sys.argv[4]
if not os.path.exists(outdir):
	os.makedirs(outdir)

cmd = f"diamond blastp --query {protein_file} --db {db_file} --outfmt 6 --out {outdir}/blast.out --threads 40 --max-target-seqs 1 --max-hsps 1"
print(cmd)

protein_dict = {}
for rec in SeqIO.parse(protein_file,'fasta'):
	name = '_'.join(rec.id.split("_")[:-1])
	if name not in protein_dict:
		protein_dict[name] = 0
	protein_dict[name] += 1

tax_dict = {}
f=open(db_ann)
lines=f.readlines()
for line in lines:
	data = line.strip().split('\t')
	tax_dict[data[0]] = ';'.join(data[1].split(";")[:5]) + '\t' + data[2]
f.close()

d={}
dd={}

f=open(outdir+'/blast.out')
lines=f.readlines()
for line in lines:
	data = line.strip().split('\t')
	name = '_'.join(data[0].split("_")[:-1])
	bitscore = float(data[-1])
	if bitscore>=0:
		name1=data[1].split("|")[0]
		if name1 in tax_dict:
			if name not in dd:
				dd[name] = {}
			tax = tax_dict[name1]
			if tax not in dd[name]:
				dd[name][tax] = 0
			dd[name][tax] += 1
f.close()

w=open(outdir+"/phage.tax.xls",'w+')
for k,vdict in dd.items():
	for k1,v1 in sorted(vdict.items(),key=lambda x:x[1],reverse=True):
		all_protein_n = protein_dict[k]
		print(k,v1,all_protein_n,k1)
		if v1/all_protein_n>=0.5:
			d[k] = 1
			w.write(k+'\t'+k1+'\n')
		break

w.close()

