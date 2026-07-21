import sys
from Bio import SeqIO

# expected call: python process_bams.py gff [samples(cov files)]
'''
WORKFLOW:
starting with:
    gff of genome
    nt level coverage
goal:
    what contigs have "good" coverage
        contiguous coverage across most of the contig
        it's ok if coverage is partial or 0 at edges
output:
    table of contigs with:
        position
        nt seq
        prot seq
        minimum coverage
        maximum coverage
        number of coverage segments
        coverage segments
we'll do better quant & DE using more specialized tools. This is just to find what contigs we should focus on
known good contigs:
    78_0dpc_S20_NODE_411_length_308_cov_39775.175097_g399_i0
'''
# constants definitions
max_allowed_coverage_segments = 1
percent_cov_threshold = 0.8

####### FUNCTION DEFINITIONS ###############
def read_in_gff(fi):
    # reads in gff file and exports dict of {feature_name:(chr, start, end)}
    # also builds a ID:mrna name dictionary for later use
    gff_dict = {}
    id2names_dict = {}
    with open(fi, 'r') as f:
        for i, line in enumerate(f):
            if not i == 0:
                line = line.strip().split('\t')
                mrna_id = float(line[8].split(';')[0].split('a')[1]) # just retains the mRNA ID numbner, optimized for JQ synthetic genome naming scheme. Casting to float for memory optimization
                # this should still work with merged synthetic genomes, since mrnaID are now [genome#].#
                name = line[8].split('Name=')[1].split(';')[0] # just get full mrna name section. second split if for future-proofing
                chromosome = line[0]
                start = int(line[3])
                end = int(line[4])
                gff_dict[mrna_id] = (chromosome, start, end) 
                id2names_dict[mrna_id] = name
                #print(line);sys.exit()
    return gff_dict, id2names_dict

def read_in_coverage(fi): 
    # reads in .coverage file from samtools depth export. 0 coverage positions not reported
    # exports dict of {(chr, pos):depth}
    c = {}
    print('loading coverage of:', fi)
    with open(fi, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip().split('\t')
            # print(i, line)
            c[(line[0], int(line[1]))] = int(line[2])
    print('loading complete')
    return c

def process_feat_coverage(l):
    # accepts list of ints reffering to feature coverage
    # returns:
    # minimum coverage
    # maximum coverage
    # number of coverage segments
    # coverage segments
    # percentage covered
    min_cov = min(l)
    max_cov = max(l)
    segs = []
    start = None
    for i, x in enumerate(l):
        if not x == 0 and start == None:
            start = i
        elif x == 0 and type(start) == int:
            segs.append((start, i-1))
            start = None
    if type(start) == int:
        segs.append((start, i))
    p_cov = sum([1 if x>0 else 0 for x in l])/len(l) #; print(p_cov) ; sys.exit()
    return min_cov, max_cov, len(segs), segs, p_cov
# print(process_feat_coverage([0,1,1,0,5,5, 0])) ;sys.exit()

def frame_translate(nt):
    # takes as input nt sequence and returns translation with least stops
    t = []
    for x in [0,1,2]:
        i = nt[x:]
        i+='N'*(3-(len(i)%3))
        t.append(i.translate())
    # print(t)
    c = 9999999
    best = ''
    for e in t:
        if e.count('*') < c:
            best = e
            c = e.count('*')
    return best

def fetch_seqs(feature_id, fasta):
    # assuming genome fasta file is in same dir
    print(fasta)
    for rec in fasta.keys():
        if rec == gff[feature_id][0]: # same chromosome
            nt_seq = fasta[rec][gff[feature_id][1]:gff[feature_id][2]-1]
            prot_seq = frame_translate(nt_seq)
           # print(nt_seq, prot_seq, gff[feature_id][1], gff[feature_id][2]-1);sys.exit()
    return str(nt_seq), str(prot_seq)

####### FUNCTION DEFINITIONS END ###############

# get gff file path
try:
    gff_fi = sys.argv[1]
except IndexError:
    print('defaulting gff file to: contigSynGen_spades_synGen.gff')
    gff_fi='contigSynGen_spades_synGen.gff'
#print('gff file is', gff_fi)
    
# load in genome sequence for later use
fasta = {}
for rec in SeqIO.parse(gff_fi.replace('.gff', '.fa'), 'fasta'):
    fasta[rec.id] = rec.seq
    
# collect .coverage files
if len(sys.argv) > 2:
    sample_fis = sys.argv[2:]
else:
    sample_fis = ['test.coverage']
    print('defaulting to these sample files', sample_fis)

print(gff_fi, sample_fis)    

# import gff file
gff, id2name = read_in_gff(gff_fi)

# import coverage data
print(sample_fis)
coverages = {}
for sample_fi in sample_fis:
    coverages[sample_fi] = read_in_coverage(sample_fi)
#print(coverages)

# begin coverage calculations/feature
# dict sctruct: {feature_name:(chr, start, end)}

data = []
header = ['mrnaID', 'name', 'nt_seq', 'best_prot_seq']
for sample in coverages.keys():
    s = sample.split('/')[-1].replace('.coverage', '').replace('.sorted', '')
    header += [s+'_min_cov', s+'_max_cov', s+'_num_segments', s+'_segments', s+'_percent_cov']
for f, feat in enumerate(gff.keys()):
    pos_range = range(gff[feat][1], gff[feat][2])
    chrom = gff[feat][0]
    nt_seq, prot_seq = fetch_seqs(feat, fasta)
    # num_segs_across_samples = []
    # percent_covs = []
    test_vals = []
    for i, sample in enumerate(coverages.keys()):
        feat_coverage = []
        for x in pos_range:
            try:
                feat_coverage.append(coverages[sample][(chrom, x)])
            except KeyError:
                feat_coverage.append(0)
        # print(feat_coverage)
        min_cov, max_cov, num_segments, segments, percent_cov = process_feat_coverage(feat_coverage)
        # if feat == 349:
        #     print(min_cov, max_cov, num_segments, segments, percent_cov) 
        #     print(pos_range)
        #     print(feat_coverage, len(feat_coverage))
        #     print(id2name[feat])
        #     print(nt_seq, len(nt_seq))
            # sys.exit()
        # num_segs_across_samples.append(num_segments)
        # percent_covs.append(percent_cov)
        test_vals.append((num_segments, percent_cov))
        if i == 0:
            line = [feat, id2name[feat], nt_seq, prot_seq, min_cov, max_cov, num_segments, segments, percent_cov]
        else:
            line += [min_cov, max_cov, num_segments, segments, percent_cov]
        # print(line)#;sys.exit()
    # print(line, num_segs_across_samples, percent_covs)
    # if min(num_segs_across_samples) <= max_allowed_coverage_segments and max(percent_covs) >= percent_cov_threshold:
    if  any([(x <= max_allowed_coverage_segments and y >= percent_cov_threshold) for x, y in test_vals]):
        print('appending', line)
        data.append(line)
    # if f % 1_000 == 0:
    # if len(data) == 100:
    #     print(feat)
    #     break
        #sys.exit()

        
# print('exiting before file rewrite');sys.exit()
print(len(data))
with open('contig_coverage.csv', 'w') as f:
    f.write(','.join(header)+'\n')
    for line in data:
        t = []
        for e in line:
            if type(e) == list:
                tt = ''
                for x in e:
                    tt += str(x[0])+'-'+str(x[1])+'; '
                t.append(tt[:-2])
            else:
                t.append(e)
        f.write(','.join([str(x) for x in t]))
        f.write('\n')