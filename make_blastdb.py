# module load Biopython
# module load ncbi-db

from Bio import SeqIO
import subprocess
import sys

# takes original nr files and exports fasta files based on keyword searches

db_fi = '/private_stores/mirror/ncbi-blastdb/20230124/db/nr'

data = {}
i=0 ; c=0
for record in SeqIO.parse(db_fi, 'fasta'):
    if 'scrofa' in record.description and 'immunoglobulin' in record.description:
        data[record.id] = (record.description, str(record.seq))
        c+=1
    if c == 1739:
        break
    i+=1 
    if i%100_000 == 0:
        print(i, c)
with open('nr_Sscrofa_Ig.fasta', 'w') as f:
    for i in data.keys():
        f.write('>'+data[i][0]+'\n')
        f.write(data[i][1])
        f.write('\n')
