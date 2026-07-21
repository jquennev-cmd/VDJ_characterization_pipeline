# module load Biopython
# module load ncbi-db

from Bio import SeqIO
import subprocess
import sys

# takes original nr files and exports fasta files based on keyword searches

db_fi = '/private_stores/mirror/ncbi-blastdb/20230124/db/nr'
species='scrofa'
keywords = ['immunoglobulin']
data = {}
i=0
for record in SeqIO.parse(db_fi, 'fasta'):
    if species in record.description.lower() and any(x in record.description.lower() for x in keywords):
        data[record.id] = (record.description, str(record.seq))
        c+=1
    i+=1 
    if i%100_000 == 0:
        print(i, c)
with open('nr_Sscrofa_Ig.fasta', 'w') as f:
    for i in data.keys():
        f.write('>'+data[i][0]+'\n')
        f.write(data[i][1])
        f.write('\n')
