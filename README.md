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

The main output directories within a sample replicate (e.g. `SAMPLE-NAME_R1`) output directory are `DnaRna_REDITOOLS-ID`, `post-processed`, and `final-annotated-SnpEff`.


