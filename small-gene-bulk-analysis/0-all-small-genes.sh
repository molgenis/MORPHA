wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.refseq_protein.faa.gz

zcat < MANE.GRCh38.v1.5.refseq_protein.faa.gz | awk '
    /^>/ {
        if (id != "") print id "," length(seq)
        id=substr($1,2)
        seq=""
        next
    }
    {
        seq=seq $0
    }
    END {
        if (id != "") print id "," length(seq)
    }
' > protein_lengths.csv

wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz

gawk -F'\t' '
    BEGIN {
        OFS=","
        print "gene,NP_identifier,NM_identifier,AA_length"
    }

    NR==FNR {
        split($0, a, ",")
        aa[a[1]] = a[2]
        next
    }

    /^#/ { next }

    {
        prot = $7
        if (prot in aa && aa[prot] <= 100)
            print $4, $7, $6, aa[prot]
    }
' protein_lengths.csv <(zcat < MANE.GRCh38.v1.5.summary.txt.gz) \
> short_MANE_proteins.csv




