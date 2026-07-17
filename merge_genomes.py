import sys



# expected call: python merge_genomes.py [GENOMES] output_name
# expects gff files to be same as genome name


genomes = sys.argv[1:-1]
if len(genomes) == 0:
    print('error in arg parsing, genomes not found')
    sys.exit()
elif len(genomes) == 1:
    print('only 1 genome specified, need more to merge')
    sys.exit()
print('genomes:', genomes)
output_fi_name = sys.argv[-1]
print(output_fi_name)

# write merged fasta
with open(output_fi_name+'.fa', 'w') as f:
    for fi in genomes:
        with open(fi, 'r') as g:
            for line in g:
                f.write(line)

# write merged annoation file
# issue: merged genome needs better geneID than 1,2,3 etc
# they get duplicated when genomes are merged
# soln: use periods denoting genome position in assembly
with open(output_fi_name+'.gff', 'w') as f:
    for x, fi in enumerate(genomes):
        with open(fi.replace('.fa', '.gff'), 'r') as g:
            for i, line in enumerate(g):
                if x > 0 and i == 0:
                    continue # skips header line on all subsequent annotation files
                else:
                    f.write(line.replace('mrna', 'mrna'+str(x)+'.'))

