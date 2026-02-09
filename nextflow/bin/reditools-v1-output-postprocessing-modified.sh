#!/bin/bash

# REDItools DNA-RNA Variant Filtering Pipeline (Modified Filtering Logic)
# Optimized for: Shallow WGS, WES, and High-Depth WGS.
# Strategy: Treat WGS data as a tool for excluding false positives (it supports novel site discovery) rather than a barrier for entry (removing known RE sites due to lack of coverage).

set -e  
set -u  
set -o pipefail 

# [Arguments and Configuration - Remains same as your original]
if [ "$#" -ne 10 ]; then
    echo "Usage: $0 <RNA_BAM> <EDITING_DIR> <GENOME_FA> <SPLICE_SITES> <NOCHR_FILE> <RMSK_GTF> <SNP_GTF> <REDIPORTAL_GTF> <THREADS> <SAMPLE_ID>"
    exit 1
fi

RNA_BAM=$1 # Aligned.sortedByCoord.out.bam
EDITING_DIR=$2 # this has to be the output folder of REDItools DNA-RNA run prefixed with DnaRna_*
GENOME_FA=$3 #"inputs/refs/GRCh38.primary_assembly.genome.fa"
SPLICE_SITES=$4 #inputs/refs/gencode.v49.primary_assembly.annotation.splicesites.txt
EXClUDE_CONTIGS=$5 #inputs/refs/excluded_contigs.txt
RMSK_GTF=$6 #inputs/refs/rmsk.sorted.gtf.gz
SNP_GTF=$7 #inputs/refs/snp151.sorted.gtf.gz
REDIPORTAL_GTF=$8 #inputs/refs/REDIPortal_hg_38_atlas.gtf
THREADS=$9
SAMPLE_ID=${10} 

# ============================================================================
# MAIN PIPELINE
# ============================================================================

echo "======================================================================="
echo "REDItools DNA-RNA Variant Filtering Pipeline"
echo "========================================================================"
echo "Starting: $(date)"
echo ""

# Step 1: Navigate to output folder
echo "[Step 01] Checking output RedItools directory..."
# print the current directory name
echo "Current directory: $(pwd)"
echo "Sample ID: ${SAMPLE_ID}"

# Check for outTable file
if ls "${EDITING_DIR}"/outTable_* 1> /dev/null 2>&1; then
    echo "Found outTable file(s):"
    ls "${EDITING_DIR}"/outTable_*
    # count the number of outTable files
    NUM_FILES=$(ls "${EDITING_DIR}"/outTable_* | wc -l)
    if [ "${NUM_FILES}" -gt 1 ]; then
        echo "Multiple outTable files found (${NUM_FILES}). This should not happen. Exiting."
        exit 1
    fi
    # copy the outTable file to the current directory with a standardized name
    cp "${EDITING_DIR}"/outTable_* "outTable_${SAMPLE_ID}"
    echo "Copied outTable file to outTable_${SAMPLE_ID}"
else
    echo "No outTable file found. Exiting."
    exit 1
fi

##############################################################################
# Step 2: We drop WGS coverage requirements. We only keep sites where RNA data shows substitution but DNA data is invariant.
# This ensures that sites with 0 WGS coverage (common in shallow WGS/WES) survive the first pass filtering.

echo "[Step 02] Keeping sites that are variants in RNA data AND invariant in DNA data..."
awk -v FS="\t" '{if ($8!="-" && $13=="-") print}' "outTable_${SAMPLE_ID}" > "outTable_${SAMPLE_ID}.filt.out" && \
###############################################################################
# Step 3: Annotate with RepeatMasker
echo "[Step 03] Annotating with RepeatMasker..." && \
AnnotateTable.py -a "${RMSK_GTF}" -n rmsk -i "outTable_${SAMPLE_ID}.filt.out" -o "outTable_${SAMPLE_ID}.filt.out.rmsk" -u && \
################################################################################
# Step 4: Annotate with dbSNP
echo "[Step 04] Annotating with dbSNP..." && \
AnnotateTable.py -a "${SNP_GTF}" -n snp151 -i "outTable_${SAMPLE_ID}.filt.out.rmsk" -o "outTable_${SAMPLE_ID}.filt.out.rmsk.snp" -u
################################################################################
# Step 5: REDIportal Annotation
echo "[Step 05] Annotating known RNA editing sites from REDIportal..."

if [ ! -f "outTable_${SAMPLE_ID}.filt.out.rmsk.snp" ]; then
    echo "Annotated file outTable_${SAMPLE_ID}.filt.out.rmsk.snp not found. Exiting."
    exit 1
fi

AnnotateTable.py -a "${REDIPORTAL_GTF}" -n ed -k R -c 1 -i "outTable_${SAMPLE_ID}.filt.out.rmsk.snp" -o "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed" -u

# Use this annotated file for subsetting
INPUT_FILE="outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed"

################################################################################
# Step 6: Create Position Sets
echo "[Step 06] Selecting RNA-seq quality sets (Set1: ≥5 RNAseq reads, 1 mismatch, Set2: ≥10 RNAseq reads, 3 mismatches, freq ≥0.1)..."
selectPositions.py -i "$INPUT_FILE" -c 5 -v 1 -f 0.0 -o "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1" && \
# count lines
NUM_SET1=$(wc -l < "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1") && \
echo "Number of positions selected in set 1: ${NUM_SET1}"

selectPositions.py -i "$INPUT_FILE" -c 10 -v 3 -f 0.1 -o "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2" && \
# count lines
NUM_SET2=$(wc -l < "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2")
echo "Number of positions selected in set 2: ${NUM_SET2}"
###############################################################################
# Step 7: Select ALU sites; where $1!="chrM" excludes mitochondrial reads, substr($16,1,3)=="Alu" selects positions in Alu elements, $17=="-" removes known SNPs, $8!="-" requires only variant positions...
echo "[Step 07] Selecting ALU sites with modified filtering logic..."
awk -v FS="\t" '{
    is_known = ($19=="ed");
    # Pass if: (Site is in REDIportal) OR (Site is Novel AND has 10+ WGS reads)
    pass_wgs = (is_known || $10>=10);
    if ($1!="chrM" && substr($16,1,3)=="Alu" && $17=="-" && $8!="-" && pass_wgs) print
}' "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1" > "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1.ALU" && \
# count lines
NUM_ALU=$(wc -l < "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1.ALU")
echo "Number of ALU sites selected: ${NUM_ALU}"

################################################################################
# Step 8: Select REP NON ALU sites; where $1!="chrM" excludes mitochondrial reads, substr($16,1,3)!="Alu" selects positions in non-Alu elements, $15!="-" excludes non-annotated sites, $15!="Simple_repeat" excludes positions in simple repeats, $15!="Low_complexity" excludes sites in low-complexity regions, $17=="-" removes known SNPs, $8!="-" requires only variant positions, $10>=10 selects sites covered by ≥10 WGS reads, $14<=0.05 considers WGS homozygous positions (with the minor allele frequency <0.05) and $9>=0.1 includes only sites with an editing level ≥0.1.....
echo "[Step 08] Selecting REP NON ALU sites with modified filtering logic..."
awk -v FS="\t" '{
    is_known = ($19=="ed");
    # Pass if: (Known) OR (Novel AND 10+ WGS reads AND low DNA variant freq)
    pass_wgs = (is_known || ($10>=10 && $14<=0.05));
    if ($1!="chrM" && substr($16,1,3)!="Alu" && $15!="-" && $15!="Simple_repeat" && $15!="Low_complexity" && $17=="-" && $8!="-" && $9>=0.1 && pass_wgs) print
}' "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2" > "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONALU" && \
# count lines
NUM_NONALU=$(wc -l < "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONALU")
echo "Number of REP NON ALU sites selected: ${NUM_NONALU}"

################################################################################
# Step 9: Select NON REP sites; where $1!="chrM" excludes mitochondrial reads, substr($16,1,3)!="Alu" selects positions in non-Alu elements, $15=="-" includes sites in non-repetitive elements, $17=="-" removes known SNPs, $8!="-" requires only variant positions, $10>=10 selects sites covered by ≥10 WGS reads, $14<=0.05 considers WGS homozygous positions (with the minor allele frequency <0.05) and $9>=0.1 includes only sites with an editing level ≥0.1....
echo "[Step 09] Selecting NON REP sites with modified filtering logic..."
awk -v FS="\t" '{=
    is_known = ($19=="ed");
    pass_wgs = (is_known || ($10>=10 && $14<=0.05));
    if ($1!="chrM" && substr($16,1,3)!="Alu" && $15=="-" && $17=="-" && $8!="-" && $9>=0.1 && pass_wgs) print
}' "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2" > "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONREP" && \
# count lines
NUM_NONREP=$(wc -l < "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONREP")
echo "Number of NON REP sites selected: ${NUM_NONREP}"

################################################################################
# Step 10: Extract known editing events
echo "[Step 10] Extracting known editing events..."

# first test if the three annotated files exist
if [ -f "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1.ALU" ] && [ -f "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONALU" ] && [ -f "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONREP" ]; then
    echo "All three annotated files found. Proceeding to extract known editing events."
else
    echo "One or more annotated files are missing. Exiting."
    exit 1
fi

cp "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set1.ALU" "${SAMPLE_ID}-ALU"
cp "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONALU" "${SAMPLE_ID}-NONALU"
cp "outTable_${SAMPLE_ID}.filt.out.rmsk.snp.ed.set2.NONREP" "${SAMPLE_ID}-NONREP"
cat "${SAMPLE_ID}-ALU" "${SAMPLE_ID}-NONALU" "${SAMPLE_ID}-NONREP" > "${SAMPLE_ID}-ALU-NONALU-NONREP" && \
awk -v FS="\t" '{if ($19=="ed") print}' "${SAMPLE_ID}-ALU-NONALU-NONREP" > "${SAMPLE_ID}-knownEditing"

################################################################################
# Step 11: The previous step retains known editing events in separate files for ALU, REP NON ALU and NON REP sites. The next steps process novel editing candidates separately for ALU and for REP NON ALU + NON REP sites.
echo "[Step 11-A] Processing novel editing candidates separately for ALU and for REP NON ALU + NON REP sites..."
echo "First, extract novel REP NON ALU and NON REP sites..."
cat "${SAMPLE_ID}-NONALU" "${SAMPLE_ID}-NONREP" > "${SAMPLE_ID}-NONALU-NONREP--novel" && \
awk -v FS="\t" '{if ($19!="ed") print}' "${SAMPLE_ID}-NONALU-NONREP--novel" > "${SAMPLE_ID}--pos.txt"
# convert to GFF
echo "Now convert NON ALU and NON REP novel sites to GFF..."
TableToGFF.py -i "${SAMPLE_ID}--pos.txt" -s -t -o "${SAMPLE_ID}--pos.gff"

# Step 11-B: Extract novel ALU sites
echo "[Step 11-B] Extracting novel ALU sites..."
awk -v FS="\t" '{if ($19!="ed") print}' "${SAMPLE_ID}-ALU" > "${SAMPLE_ID}--pos-ALU.txt"
# convert to GFF
echo "Now convert ALU novel sites to GFF..."
TableToGFF.py -i "${SAMPLE_ID}--pos-ALU.txt" -s -t -o "${SAMPLE_ID}--pos-ALU.gff"
###############################################################################
# Step 12: Run REDItoolDnaRna.py on ALU sites with stringent criteria to recover potential editing events
echo "[Step 12] Run REDItoolDnaRna.py on ALU sites (stringent)..."

# test if gff.gz file exists first
if [ ! -f "${SAMPLE_ID}--pos-ALU.sorted.gff.gz" ];
then
    echo "GFF file ${SAMPLE_ID}--pos-ALU.sorted.gff.gz not found. Exiting."
    exit 1
fi

echo "GFF file found. Proceeding with REDItoolDnaRna.py..."

REDItoolDnaRna.py -s 2 -g 2 -S -t ${THREADS} \
    -i "${RNA_BAM}" \
    -f "${GENOME_FA}" \
    -c 5,5 -q 30,30 -m 255,255 -O 5,5 \
    -p -u -a 11-6 -l -v 1 -n 0.0 -e \
    -T "${SAMPLE_ID}--pos-ALU.sorted.gff.gz" \
    -w "${SPLICE_SITES}" \
    -k "${EXClUDE_CONTIGS}" \
    -R -o "${SAMPLE_ID}-firstalu"
###############################################################################
# Step 13: Run REDItoolDnaRna.py on REP NON ALU and NON REP sites using stringent criteria to recover RNA reads harboring reference mismatches
echo "[Step 13] Run REDItoolDnaRna.py on REP NON ALU and NON REP sites (less stringent)..."
# test if gff.gz file exists first
if [ ! -f "${SAMPLE_ID}--pos.sorted.gff.gz" ];
then
    echo "GFF file ${SAMPLE_ID}--pos.sorted.gff.gz not found. Exiting."
    exit 1
fi

echo "GFF file found. Proceeding with REDItoolDnaRna.py..."
REDItoolDnaRna.py -s 2 -g 2 -S -t ${THREADS} \
    -i "${RNA_BAM}" \
    -f "${GENOME_FA}" \
    -c 10,10 -q 30,30 -m 255,255 -O 5,5 \
    -p -u -a 11-6 -l -v 3 -n 0.1 -e \
    -T "${SAMPLE_ID}--pos.sorted.gff.gz" \
    -w "${SPLICE_SITES}" \
    -k "${EXClUDE_CONTIGS}" \
    --reads -R --addP -o "${SAMPLE_ID}-first"
###############################################################################
# Step 14: Run pblat on RNAseq reads harboring reference mismatches and select for multimapping reads
echo "[Step 14] Running pblat to identify multi-mapping reads with reference mismatches..."
FIRST_DNARNA=$(find "${SAMPLE_ID}-first" -type d -name "DnaRna_*" | head -n1)

pblat -t=dna -q=rna -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 "${GENOME_FA}" "${FIRST_DNARNA}/outReads_$(basename ${FIRST_DNARNA} | sed 's/DnaRna_//')" "${FIRST_DNARNA}/reads.psl"

# first check reads.psl file exists or not
if [ ! -f "${FIRST_DNARNA}/reads.psl" ]; then
    echo "reads.psl file not found in ${FIRST_DNARNA}. Exiting."
    exit 1
fi

readPsl.py "${FIRST_DNARNA}/reads.psl" "${FIRST_DNARNA}/badreads.txt"

###############################################################################
# Step 15: Extract and deduplicate reads
echo "[Step 15] Extracting and deduplicating reads..."
sort -k1,1 -k2,2n -k3,3n "${FIRST_DNARNA}/outPosReads_$(basename ${FIRST_DNARNA} | sed 's/DnaRna_//')" | mergeBed > bed && \
samtools view -@ ${THREADS} -L bed -h -b "${RNA_BAM}" > "${SAMPLE_ID}_bed.bam" && \
samtools sort -@ ${THREADS} -n "${SAMPLE_ID}_bed.bam" -o "${SAMPLE_ID}_bed_ns.bam" && \
samtools fixmate -@ ${THREADS} -m "${SAMPLE_ID}_bed_ns.bam" "${SAMPLE_ID}_bed_ns_fx.bam" && \
samtools sort -@ ${THREADS} "${SAMPLE_ID}_bed_ns_fx.bam" -o "${SAMPLE_ID}_bed_ns_fx_st.bam" && \
samtools markdup -r -@ ${THREADS} "${SAMPLE_ID}_bed_ns_fx_st.bam" "${SAMPLE_ID}_bed_ns_fx_st_dedup.bam" && \
samtools index "${SAMPLE_ID}_bed_ns_fx_st_dedup.bam"
###############################################################################
# Step 16: Re-run REDItoolDnaRna.py with deduplicated reads

# check input file
if [ ! -f "${SAMPLE_ID}_bed_ns_fx_st_dedup.bam" ];
then
    echo "Deduplicated BAM file ${SAMPLE_ID}_bed_ns_fx_st_dedup.bam not found. Exiting."
    exit 1
fi

echo "[Step 16] Re-run REDItoolDnaRna.py on REP NON ALU and NON REP sites using stringent criteria,
deduplicated reads and mis-mapping info..."

REDItoolDnaRna.py -s 2 -g 2 -S -t ${THREADS} -i "${SAMPLE_ID}_bed_ns_fx_st_dedup.bam" -f "${GENOME_FA}" -c 10,10 -q 30,30 -m 255,255 -O 5,5 -p -u -a 11-6 -l -v 3 -n 0.1 -e -T "${SAMPLE_ID}--pos.sorted.gff.gz" -w "${SPLICE_SITES}" -k "${EXClUDE_CONTIGS}" -b "${FIRST_DNARNA}/badreads.txt" -R --rmIndels -o "${SAMPLE_ID}-second"

###############################################################################
# COMPLETION MESSAGE

echo ""
echo "========================================================================"
echo "Pipeline completed successfully!"
echo "========================================================================"
echo "Finished: $(date)"
echo ""


echo "Pipeline finished successfully: $(date)"
