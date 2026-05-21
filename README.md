# MORPHA
MOlecular effect Range of germline Pathogenic variation in global protein HAplotypes

Example run (needs to be in this order):
1. ./genename2uniqproteinhaplo.sh MEFV
2. ./genename2pathogenicvariants.sh MEFV
3. ./genename2pathohaplo.sh MEFV


### notes
UCSC API: 0-based
https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=16;start=3243404;end=3243407
-> CAT

ClinVar VCF: 1-based
POS 3243405 -> CAT

