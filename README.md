# MORPHA
MOlecular effect Range of germline Pathogenic variation in global protein HAplotypes

Example run (in this order):
1. `./genename2uniqproteinhaplo.sh <gene> <RefSeq protein ID> <front signal peptide cutoff, optional, default = 0>`
   - Using GRCh38 and 1000G data, reconstruct and group all 1KG haplotypes.
   - If front signal peptide cutoff is provided, the first N residues will be removed from each output protein sequence.
2. `./genename2pathogenicvariants.sh <gene> <RefSeq protein ID> `
   - Obtain pathogenic variants from ClinVar and prepare for merging with 1000G haplotypes.
3. `./genename2pathohaplo.sh <gene> <localized vcf.gz> <front signal peptide cutoff, optional, default = 0>`
   - Spike pathogenic variants into 1KG haplotypes and group them.
   - A custom variant vcf.gz can be supplied here, skipping step 2.
     - Made with `bgzip -c abc.vcf > abc.vcf.gz` followed by `tabix -p vcf abc.vcf.gz`.
     - BEWARE: the VCF `POS` field must be calculated as: `(hg38 location) - (GFF transcript start location) + 1`.
     - Also, `ID` field may NOT include `:`.
4. `./genename2foldedhaplo.sh <gene>`
   - Fold all unique haplotype proteins, 1KG and spiked.
   - Expects [esmfold-server](https://github.com/joerivandervelde/esmfold-server) to run on your network, default address http://192.168.2.9:8000, adjust accordingly.
5. Merge haplotype data with IGSR sample metadata
   - Run `merge_haplo_with_igsrmd.R` (input gene inside script)
   - Inputs: `<gene>_protein_haplotypes.fa`, `<gene>_protein_path_haplotypes.fa` and `viz/igsr_samples.tsv`
   - Output: `<gene>_protein_haplotypes_matched_to_igsr.tsv`
5. Quantify protein structural properties of choice
   - Using FoldX, Rosetta, GROMACS, SPURS, etc etc
6. Merge your quantifications with sample haplotype data
   - Use `getQuant-FoldX.R` for instance
   - Inputs: `<gene>_protein_haplotypes_matched_to_igsr.tsv` and your protein quantifications (e.g. `*.fxout` files)
   - Output: `<gene>_protein_haplotypes_quant.tsv`
6. Visualize results
   - Use `world-plot-1KG-subpop.R`
   - Inputs: `<gene>_protein_haplotypes_quant.tsv` and `viz/igsr_samples.tsv`
   - Outputs: world plot according to parameters


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
### FoldX
```
for i in {9..57}; do
   /Applications/FoldX/foldx5_Mac_0/foldx_20261231 --command=RepairPDB --pdb=MEFV_protein_variant_${i}.pdb
   /Applications/FoldX/foldx5_Mac_0/foldx_20261231 --command=Stability --pdb=MEFV_protein_variant_${i}_Repair.pdb
done
for i in {1..275}; do
   /Applications/FoldX/foldx5_Mac_0/foldx_20261231 --command=RepairPDB --pdb=MEFV_protein_path_variant_${i}.pdb
   /Applications/FoldX/foldx5_Mac_0/foldx_20261231 --command=Stability --pdb=MEFV_protein_path_variant_${i}_Repair.pdb
done
```
### Rosetta
```
for i in {1..7}; do
   /Applications/rosetta.binary.m1.release-408/main/source/bin/relax.static.macosclangrelease \
      -in:file:s TTR_protein_path_variant_${i}_Repair.pdb \
      -relax:fast \
      -nstruct 10 \
      -out:file:scorefile TTR_protein_path_variant_${i}_Repair.sc \
      -out:path:pdb /tmp/
   /Applications/rosetta.binary.m1.release-408/main/source/bin/relax.static.macosclangrelease \
      -in:file:s TTR_protein_variant_${i}_Repair.pdb \
      -relax:fast \
      -nstruct 10 \
      -out:file:scorefile TTR_protein_variant_${i}_Repair.sc \
      -out:path:pdb /tmp/
done
```

## Ubuntu WSL install and run

### install the following
```
sudo apt install samtools # old but might be OK, else compile from src like bcftools below
sudo apt install gffread
sudo apt install python3-pip
sudo apt install tabix
pip3 install CrossMap
```

### install latest BCFtools from source because default Ubuntu install is too old
```
sudo apt-get update
sudo apt-get install -y build-essential bzip2 libbz2-dev xz-utils liblzma-dev zlib1g-dev libcurl4-gnutls-dev libssl-dev
wget https://github.com/samtools/bcftools/releases/download/1.23.1/bcftools-1.23.1.tar.bz2
tar -xjvf bcftools-1.23.1.tar.bz2
cd bcftools-1.23.1
./configure --prefix=$HOME/.local
make
make install
export PATH="$HOME/.local/bin:$PATH"
```

# before running, fix line endings if checkout was done under Windows
```
sed -i 's/\r//' genename2uniqproteinhaplo.sh
sed -i 's/\r//' genename2pathogenicvariants.sh
sed -i 's/\r//' genename2pathohaplo.sh
sed -i 's/\r//' genename2foldedhaplo.sh
```

# MEFV pathogenic variants from INFEVERS (instead of ClinVar)
- Download TSV from https://infevers.umai-montpellier.fr/web/download.php?n=1
- Filter lines if needed (i.e. select only "Pathogenic", "Benign", etc)
- Grab cDNA and format variants like this: NM_000243.3:c.724A>G
- Run variants through VEP at https://www.ensembl.org/Homo_sapiens/Tools/VEP/
- Download result as VCF
- Get genome start coordinate from MEFV MANE GFF, line containing 'transcript', transcript_start = 3242027
- Calculate and replace by local start coordinates: genomic_start - transcript_start + 1
- Replace chromosome '16' by 'chr16'
- Remove 'NM_000243.3:' from ID field
- Compress: `bgzip -c MEFV-infevers-pathogenic-local.vcf > MEFV-infevers-pathogenic-local.vcf.gz`
- Index: `tabix -p vcf MEFV-infevers-pathogenic-local.vcf.gz`
- Run haplotype builder with these variants:  `./genename2pathohaplo.sh MEFV /path/to/MEFV-infevers-pathogenic-local.vcf.gz`

# Trio of genes
Analysis of 3 hallmark pathogenic mutations in known disease genes TTR-V50M, GBA1-N409S and MLH1-L622H.
- Research question: When does genetic background modulate pathogenicity?
- Criteria: pathogenic missense mutation in relatively small protein, variable/reduced penetrance, ancestral/founder haplotypes, haplotypes proposed to explain phenotypic/penetrance differences, protein-stability/folding mechanism.
Notes:
- TTR has a 20-residue signal peptide to be removed (https://www.uniprot.org/uniprotkb/P02766/entry)
- GBA1 has a 39-residue signal peptide to be removed (https://www.uniprot.org/uniprotkb/P04062/entry)
- MLH1 does not have a signal peptide
```
./genename2uniqproteinhaplo.sh TTR 20
./genename2uniqproteinhaplo.sh GBA1 39
./genename2uniqproteinhaplo.sh MLH1

./genename2pathohaplo.sh TTR results/TTR/TTR-V50M.vcf.gz 20
./genename2pathohaplo.sh GBA1 results/GBA1/GBA1-N409S.vcf.gz 39
./genename2pathohaplo.sh MLH1 results/MLH1/MLH1-L622H.vcf.gz

./genename2foldedhaplo.sh TTR
./genename2foldedhaplo.sh GBA1
./genename2foldedhaplo.sh MLH1
```
Control sample: TP53-R175H
- Pathogenic missense mutation in small protein but no established R175H-specific ancestral-haplotype mechanism explaining variation in penetrance.
```
./genename2uniqproteinhaplo.sh TP53
./genename2pathohaplo.sh TP53 results/TP53/TP53-R175H.vcf.gz
./genename2foldedhaplo.sh TP53
```
Then run R scripts merge_haplo_with_igsrmd.R, getQuant-FoldX.R, etc.

Birgit suggestion: CHEK2 c.1100delC cancer risk

NOTE that chromosome X needs manual download and rename:
```
cd data
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi
mv 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
mv 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
```
