import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIXML
# module load BLAST+


import itertools
import random
import sys
import subprocess

# define custom blast database file location
blast_db_fi='nr_Sscrofa_Ig.fasta'

# processes trinity output rather than inchworm
di = 'spades_'

fis = sys.argv[1:-1]
out_fi_name = sys.argv[-1]
print('processing:', fis)
print('outputting to:', out_fi_name)
intended_product_size = 361 # update as needed
size_threshold = 0.4 # update as needed
min_prod_size = intended_product_size*(1-size_threshold)
max_prod_size = intended_product_size*(1+size_threshold)
seqs = {}
for fi in fis:
    path = di+fi+'/transcripts.fasta'
    try:
        for rec in SeqIO.parse(path, 'fasta'):
            # name = fi.split('.')[0][8:] + '_' + rec.id.replace('TRINITY_', '')
            name = fi + '_' + rec.id
            if min_prod_size <= len(rec.seq) <= max_prod_size:
                if not 'N' in str(rec.seq): # some contigs have poly Ns for some reason
                    seqs[name] = str(rec.seq)
    except FileNotFoundError:
        print('error, could not find:', path)#; input()
        continue
# sys.exit()
        
rev_d = {}
for k in seqs.keys():
        try:
            rev_d[seqs[k]] = rev_d[k] + [k]
            print('foo')
        except KeyError:
            rev_d[seqs[k]] = [k]
t = {}
for k in rev_d.keys():
    t['_'.join(rev_d[k])] = k
seqs = t
# sys.exit()

print('total number of contigs:', len(seqs.keys()))
# BLAST analysis


def allFrameTranslate(seq):
    # expects string to be translated in all 3 frames. Does not do rev-comp
    translations = []
    frames = [0,1,2]
    for x in frames:
        #print(x)
        translations.append(str(Seq(seq[x:((len(seq)-x)//3)*3+x]).translate()))# trims off miniumum excess of sequence to supress Bio warning
    return translations

# test code
#print(allFrameTranslate('TGCCACCATGAAGAAGACATAGCGTTTTTGCTCGCCTTGATGTTTGTTTTTTCAATAGC'))

def BLASTcall(query, db=blast_db_fi):
    # calls blastp and outputs xml file
    # blastp -query test.fasta -db nr_Sscrofa.fasta -out blastP2.xml -outfmt 5
    # print('blasting', query)
    with open('temp.fasta', 'w') as f:
        f.write('>temp\n'+query+'\n')
    blast_cmd = ['blastp', '-query', 'temp.fasta', '-db', db, '-out', 'blastP.xml', '-max_target_seqs', '1800', '-outfmt', '5']
    subprocess.call(blast_cmd)
    return

def BLASTread(fi):
    Ig_detected = 0
    with open('blastP.xml', 'r') as handle:
        recs = NCBIXML.parse(handle)
        for rec in recs:
            for alignment in rec.alignments:
                if 'immuno' in alignment.title:
                    print('Ig detected') 
                    if alignment.hsps[0].align_length/alignment.length >= 0.75:
                        print('passes alignment length threshold')
                        Ig_detected +=1 
                    else:
                        print('does not pass alignment threshold')
    return Ig_detected
# with open('temp', 'w') as f:
#     f.write('\n'.join(seqs.keys()))
# sys.exit()

data = []
for x, k in enumerate(seqs.keys()):
    # k= '77_0dpc_S17_NODE_11_length_1764_cov_33.171629_g9_i0'
    if x % 50 == 0:
        print('BLAST processed', x, 'contigs')
    translations = []
    translations+=allFrameTranslate(seqs[k])
    translations+=allFrameTranslate(str(Seq(seqs[k]).reverse_complement()))
    translate_data = []
    for i, prot in enumerate(translations):
        # print(row.id)
        if '*' in prot:
            translate_data.append([prot] + [-99.0, i])
        else:
            #BLAST against Sscrofa 
            #print(prot)
            BLASTcall(prot)
            Igs = BLASTread('temp.fasta')
            translate_data.append([prot] + [Igs, i])
        #     if k == '77_0dpc_S17_NODE_11_length_1764_cov_33.171629_g9_i0':
        #         print(Igs)
        # print('frame', i) ;input()
    # if k =='77_0dpc_S17_NODE_11_length_1764_cov_33.171629_g9_i0':
    #     sys.exit()
            # sys.exit()
    # print (translate_data)
    #print(translate_data)
    best_translation = max(translate_data, key=lambda sublist: sublist[1])
    data.append([k, seqs[k]]+best_translation)
    # if len(data) % 50 == 0:
    #     print(len(data));break
print(data)
df = pd.DataFrame(data, columns=['id', 'seq', 'best_translation', 'score', 'frame'])
df.to_csv('SPAdes_contigs_completeDataTable_'+out_fi_name+'.csv')
# BlAST analysis complete

print(df.shape)
df = df[df['score'] > 0]
print(df.shape)
new_seq = []
for i, row in df.iterrows():
    if row['frame'] > 2:
        new_seq.append(str(Seq(row['seq']).reverse_complement()))
    else:
        new_seq.append(row['seq'])
df['seq'] = new_seq
df = df.loc[df['score'] > 0] # need to remove contigs which don't pass filtering

def sliding_window(s, n):
    out = []
    for i in range(len(s)-n+1):
        out.append(s[i:i+n])
    return set(out)
#x = sliding_window('gaggtgaagctggtggagtctggaggaggcctggtgcagcctggggggtctctgagactctcctgtgtcggctctggattcaccttcagtagtacctacattaactgggtccgccaggctccagggaaggggctggagtggctggcagctattagtactagtggtggtagcacctactacgcagactctgtgaagggccgattcaccatctccagagacgactcccagaacacggcctatctgcaaatgaacagcctgagaaccgaagacacggcccgctattactgtgcaacagggaattgctatagctatggtgctagttgctatagtgacgccatggctacttagattcgtggggccagggcatcctggtcaccgtctcctcag', 20)
# for e in x:
#     print(e)

def gen_kmers(seqs, n):
    print('generating all sample kmers')
    kmers = []
    for s in seqs:
        k = sliding_window(s, n)
        kmers = list(set(kmers) | k)
    print(len(kmers), 'total kmers')
    print('returning all sample kmers')
    return kmers
#gen_kmers(['gaggtgaagctggtggagtctggaggaggcctg','aaggtgaagctggtggagtctggaggaggcctg'], 20)
sample_kmers = gen_kmers(df.loc[df['score'] > 0, 'seq'].to_list(), 20)
#print(sample_kmers)

def gen_spacer_seq(kmers, l):
    # kmer set subtraction process takes too much RAM
    # bases = ['a','t','g','c']
    # all_kmers = set([''.join(p) for p in itertools.product(bases, repeat=20)])
    # allowed_kmers = all_kmers - kmers
    
    # trying random seq generation method
    valid_kmers = [] # kmers not found in sample set
    while len(valid_kmers) < 30:
        prev_len = len(valid_kmers)
        x = ''.join(random.choices('atcg', k=20))
        #print(x)
        if not x in kmers or not x in valid_kmers:
            valid_kmers.append(x)
        if len(valid_kmers) > prev_len:
            print(len(valid_kmers))
    switch = 0 ; c = 0
    while switch == 0:
        c+=1
        syn_seq = 'tttttt'+''.join(random.choices(valid_kmers, k=10))+'tttttt'
        for s in valid_kmers:
            if s in syn_seq:
                print('re-generate syn_seq', c)
                continue
        switch = 1
    print(syn_seq)
    return syn_seq
    

#gen_spacer_seq(['catcctggtcaccgtctcct', 'ccccagctggctggattact'], 200)
syn_seq = gen_spacer_seq(sample_kmers, 200)

genome_seq = ''+syn_seq
gff_data = []
for i,c in enumerate(df['seq']):
    # print(c)
    l0 = len(genome_seq)+1
    genome_seq += c
    l1 = len(genome_seq)+1
    genome_seq += syn_seq
    #gff_data.append(['42', 'processed_transcript', 'transcript', str(l0), str(l1), '.', '+', '.'])
    gff_data.append(['42', '.', 'mRNA', str(l0), str(l1), '.', '+', '.', 'ID=mrna'+str(i+1)+';Name='+df.iloc[i]['id'].replace(';', '_')])

with open('SPAdes_contigSynGen_'+out_fi_name+'.fa', 'w') as f:
    f.write('>42\n')
    f.write(genome_seq+'\n')
with open('SPAdes_contigSynGen_'+out_fi_name+'.gff', 'w') as f:
    f.write('##gff-version 3\n')
    for line in gff_data:
        f.write('\t'.join(line))
        f.write('\n')

# print(df)
