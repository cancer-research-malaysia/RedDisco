note: the refMrna.fa was a RefSeq-based hg38 transcriptome fasta downloaded directly from `wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz`

I ran `awk '/^>/ {print $1 "." $2; next} {print}' refMrna.fa > hg38_refseq_Mrna_clean.fa` to reformat the header lines for each transcript to match the ones used in SnpEff
