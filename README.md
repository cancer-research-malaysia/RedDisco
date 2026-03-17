# REdDisco – RNA Editing Discovery @ CRMY

This repository documents the analysis files and scripts to explore RNA editing events.

## Tools Tested

We tested 3 different tools; Reditools v1, Reditools v3, and SPRINT, for RNA editing event discovery. 

### Reditools v1

1. Rewrote `rediportal2recoding.py` to be compatible with the hg38 version of REDIportal's `TABLE_1.txt`. (`rediportal2recoding-v2.py`)
	- This needs to be incorporated manually into the dockerized version


_________________________________________________________________________________________________________

## Output Files

The REdDisco pipeline currently produced a specific output directory, which follows the following format (`SAMPLE-NAME_REP`):

```
└── SAMPLE-NAME_R1
    ├── DnaRna_SAMPLE-NAME_R1
    │   └── parameters.txt
    ├── final-annotated-SnpEff
    │   ├── SAMPLE-NAME_R1_gene_level_editing_subset_protein_coding_gt2_sites.tsv
    │   ├── SAMPLE-NAME_R1_gene_level_editing.tsv
    │   ├── SAMPLE-NAME_R1-AG-Subs-Only-Sites-Freq-10pct.ann.vcf
    │   ├── SAMPLE-NAME_R1-AG-Subs-Only-Sites-Freq-10pct.vcf
    │   ├── SAMPLE-NAME_R1.known_sites_freq10pct_genic.tsv
    │   └── SAMPLE-NAME_R1.known_sites_freq10pct.tsv
    └── postprocessed
        ├── final_output_files
        │   ├── SAMPLE-NAME_R1-AG-Subs-Only-Sites-Freq-10pct.tsv
        │   ├── SAMPLE-NAME_R1-allEditing.tsv
        │   ├── SAMPLE-NAME_R1-knownEditing-labeled.tsv
        │   └── SAMPLE-NAME_R1-novelEditing.tsv
        ├── SAMPLE-NAME_R1-ALU--novelEditing.txt
        ├── SAMPLE-NAME_R1-first
        │   └── DnaRna_743461156
        │       ├── badreads.txt
        │       ├── outPosReads_743461156
        │       ├── outReads_743461156
        │       ├── parameters.txt
        │       └── reads.psl
        ├── SAMPLE-NAME_R1-firstalu
        │   └── DnaRna_149009986
        │       └── parameters.txt
        ├── SAMPLE-NAME_R1-knownEditing
        ├── SAMPLE-NAME_R1-NONALU-NONREP--novelEditing.txt
        ├── SAMPLE-NAME_R1-second
        │   └── DnaRna_636381558
        │       └── parameters.txt
        └── summary_stats_files
            ├── SAMPLE-NAME_R1-editingSummary.txt
            └── SAMPLE-NAME_R1-FinalDiscoverySummary.txt
```

1. The main output directories within a sample replicate (e.g. `SAMPLE-NAME_R1`) output directory are `DnaRna_REDITOOLS-ID`, `post-processed`, and `final-annotated-SnpEff`.

2. The `post-processed` folder contains various plaintext files containing either all editing sites that pass the stringent filtering, novel editing sites, and AG-edited sites (with editing frequency above or equal to 10%).

3. The `final-annotated-SnpEff` folder contains various files that contain only the AG-edited sites. The `tsv` files are tab-delimited files that can be opened in Excel or text editor of your choice. 

4. In the `final-annotated-SnpEff` folder, the `.known_sites_freq10pct.tsv` files contain only the bona fide known AG-edited sites passing the editing frequency filter, while the `_genic.tsv` version is the subset of sites that are NOT intergenic.

5. Finally, the `_gene_level_editing.tsv` is the file that contains the aggregated data of sites at gene level. The `_gene_level_editing_subset_protein_coding_gt2_sites.tsv` file is the subset of the former, containing only the genes that have more than two AG-edited sites AND are themselves protein coding genes. 


