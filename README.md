# MORPHA
MOlecular effect Range of germline Pathogenic variation in global protein HAplotypes

Example run (needs to be in this order):
1. ./genename2uniqproteinhaplo.sh MEFV
   - Using GRCh38 and 1000G data, reconstruct and group all 1KG haplotypes.
2. ./genename2pathogenicvariants.sh MEFV
   - Obtain all known pathogenic variants
3. ./genename2pathohaplo.sh MEFV
   - Spike pathogenic variants into 1KG haplotypes and group them.
4. ./genename2foldedhaplo MEFV
   - Fold all unique haplotype proteins, 1KG and spiked.
   - Expects [esmfold-server](https://github.com/joerivandervelde/esmfold-server) to run on your network, default address http://192.168.2.9:8000, adjust accordingly.

General script dependencies:
- samtools 1.23.1 (also tested with 1.13, seems OK)
- bcftools 1.23.1
- CrossMap 0.7.0
- gffread 0.12.9
- tabix 1.17 (also tested with 1.13, seems OK)


### notes
UCSC API: 0-based
https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=16;start=3243404;end=3243407
-> CAT

ClinVar VCF: 1-based
POS 3243405 -> CAT

## Quantify protein stability example
```
for i in {1..57}; do
   echo "--- Processing Variant $i ---"
   /Applications/FoldX/foldx5_Mac_0/foldx_20261231 --command=RepairPDB --pdb=MEFV_protein_variant_${i}.pdb
   /Applications/FoldX/foldx5_Mac_0/foldx_20261231 --command=Stability --pdb=MEFV_protein_variant_${i}_Repair.pdb
done
```

## Ubuntu WSL install and run

### install the following
sudo apt install samtools # old but might be OK, else compile from src like bcftools below
sudo apt install gffread
sudo apt install python3-pip
sudo apt install tabix
pip3 install CrossMap

### install latest BCFtools from source because default Ubuntu install is too old
sudo apt-get update
sudo apt-get install -y build-essential bzip2 libbz2-dev xz-utils liblzma-dev zlib1g-dev libcurl4-gnutls-dev libssl-dev
wget https://github.com/samtools/bcftools/releases/download/1.23.1/bcftools-1.23.1.tar.bz2
tar -xjvf bcftools-1.23.1.tar.bz2
cd bcftools-1.23.1
./configure --prefix=$HOME/.local
make
make install
export PATH="$HOME/.local/bin:$PATH"

# before running, fix line endings if checkout was done under Windows
sed -i 's/\r//' genename2uniqproteinhaplo.sh
sed -i 's/\r//' genename2pathogenicvariants.sh
sed -i 's/\r//' genename2pathohaplo.sh
sed -i 's/\r//' genename2foldedhaplo.sh
