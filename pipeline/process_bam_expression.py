import pandas as pd
import sys



sample_dir=sys.argv[1]
genome_header=sys.argv[2]
samples=sys.argv[3:]

# test variables
# samples = ['67_0dpi_S4_newSamples+originalSamples.sorted_SPAdes_contigSynGen_mergedSamples.readCounts', '68_0dpi_S8_newSamples+originalSamples.sorted_SPAdes_contigSynGen_mergedSamples.readCounts', '68_14dpi_S7_newSamples+originalSamples.sorted_SPAdes_contigSynGen_mergedSamples.readCounts', '68_21dpi_S27_newSamples+originalSamples.sorted_SPAdes_contigSynGen_mergedSamples.readCounts', '68_7dpi_S28_newSamples+originalSamples.sorted_SPAdes_contigSynGen_mergedSamples.readCounts', '77_0dpc_S17_newSamples+originalSamples.sorted_SPAdes_contigSynGen_mergedSamples.readCounts', '77_7dpc_S35_newSamples+originalSamples.sorted_SPAdes_contigSynGen_mergedSamples.readCounts']
# genome_header='SPAdes_contigSynGen_'

df = pd.read_csv('contig_coverage.csv')

# add read counts
d_readCounts = {}
for s in samples:
    rcs = {}
    prev_id = 0
    with open(sample_dir+'/'+s, 'r') as f:
        for i, line in enumerate(f):
            if i > 1:
                line = line.strip().split('\t')
                mrnaID = float(line[0])
                # if not mrnaID - prev_id == 1:
                #     print('records are not numeric');break
                # prev_id = mrnaID
                rcs[mrnaID] = int(line[-1])
        d_readCounts[s.replace('.readCounts','')] = rcs

# construct [[sample counts], [], ...] for each sample matching mrnaID
# toAppend = [[] for x in ]
header = [x for x in d_readCounts.keys()]
toAppend = []
for s in header:
    t=[]
    for i, row in df.iterrows():
        try:
            t.append(d_readCounts[s][row.iloc[0]])
        except KeyError:
            print('shouldn\'t be here');sys.exit()
    toAppend.append(t)

df = pd.concat([df, pd.DataFrame(dict(zip([x+'.readCounts' for x in header], toAppend)))], axis=1)

# df.to_csv('temp.csv');sys.exit()
# convert readCounts based on coverage
'''
test expression from previous script
 # if min(num_segs_across_samples) <= max_allowed_coverage_segments and max(percent_covs) >= percent_cov_threshold:
    if  any([(x <= max_allowed_coverage_segments and y >= percent_cov_threshold) for x, y in test_vals])
'''

header2 = [x+'.covAdjustedCounts' for x in header]
toAppend =[[] for x in header]
for i, row in df.iterrows():
    for x, sample in enumerate(header):
        s = sample.replace('.sorted_SPAdes_contigSynGen_mergedSamples', '')
        print(s+'_num_segments')
        print(s+'.sorted.bam.coverage_percent_cov')
        t = row[s+'_num_segments'] < 2 and row[s+'_percent_cov'] >= 0.8
        # print(t);sys.exit()
        if t:
            toAppend[x].append(row[sample+'.readCounts'])
        else:
            toAppend[x].append(0)
df = pd.concat([df, pd.DataFrame(dict(zip(header2, toAppend)))], axis=1)


# normalize readCounts
# normalization is needed because there is significant (~5x) total readCount differences between samples
# using CPM 
# normalize using original readCounts, not coveraged adjusted ones:
# logic: using coveraged adjusted counts will lower the denominator in calculation, inflating expression values for samples with not a lot of well covered contigs.
# repeat coverage adjustment after calculation
# actually, don't have to repeat, just use totalCounts from original calculation

toAppend = []
toAppend2 = []
for sample in header:
    totalCounts = sum(df[sample+'.readCounts'])
    t = []
    t2 = []
    for i, row in df.iterrows():
        t.append(round(row[sample+'.readCounts']*1000000.0/totalCounts, 3))
        t2.append(round(row[sample+'.covAdjustedCounts']*1000000.0/totalCounts, 3))
        # print(          row[sample+'.covAdjustedCounts'])
        # print(t2); sys.exit()
    toAppend.append(t)
    toAppend2.append(t2)
df = pd.concat([df, pd.DataFrame(dict(zip([x+'.CPM' for x in header], toAppend)))], axis=1)
df = pd.concat([df, pd.DataFrame(dict(zip([x+'.covAdjustedCPM' for x in header], toAppend2)))], axis=1)


df.to_csv('USDA_pigImmunity_kmerContig_masterTable.csv', index = False)

t = [x.replace('_'+genome_header,'').replace('.sorted.bam','') for x in df.columns]
df.columns=t
todrop = [i for i, x in enumerate(t) if 'coverage' in x]
todrop
df2 = df.drop(df.columns[todrop], axis=1)
df2.columns
df2.to_csv('USDA_pigImmunity_kmerContig_expTable.csv', index = False)
