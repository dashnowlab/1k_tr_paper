import gzip
import sys

if sys.argv[1].endswith('.gz'):
    fh = gzip.open(sys.argv[1], 'rt')
else:
    fh = open(sys.argv[1], 'r')

out = ''
if sys.argv[1].endswith('.bed.gz'):
    out = sys.argv[1][:-7] + '.atarva.bed'
elif sys.argv[1].endswith('.bed'):
    out = sys.argv[1][:-4] + '.atarva.bed'
else:
    out = sys.argv[1] + '.atarva.bed'

out_fh = open(out, 'w')

for line in fh:
    chrom, start, end, motif, _ = line.strip().split('\t')
    print(chrom, start, end, motif, len(motif), f'{chrom}-{start}-{end}-{motif}', sep='\t', file=out_fh)

fh.close()
out_fh.close()
