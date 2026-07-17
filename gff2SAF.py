import sys


'''
saf format example:
GeneID	Chr	Start	End	Strand
497097	chr1	3204563	3207049	-
497097	chr1	3411783	3411982	-
497097	chr1	3660633	3661579	-
gff sample line:
42      .       mRNA    213     1368    .       +       .       ID=mrna1;Name=71_0dpi_S19_NODE_9_length_1155_cov_59527.853261_g8_i0
'''


header = ['GeneID', 'Chr', 'Start', 'End', 'Strand']

# gff = 'contigSynGen_spades_synGen_offset1_largeIndex.gff'
gff = sys.argv[1]+'.gff'
data = []
with open(gff, 'r') as f:
    for i, line in enumerate(f):
        if not i ==0:
            line = line.strip().split('\t')
            data.append([line[-1].split(';')[0].split('mrna')[1], line[0], line[3], line[4], '+'])
            
with open(gff[:-4]+'.saf', 'w') as f:
    f.write('\t'.join(header)+'\n')
    for line in data:
        f.write('\t'.join(line)+'\n')