#!/bin/bash

# REDItools DNA-RNA variant filtering pipeline
# Based on Nature Protocols procedure

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS
# ============================================================================

# Input directories
RMSK="inputs/refs/rmsk.sorted.gtf.gz"
SNP_DIR="inputs/refs/snp151.sorted.gtf.gz"

ALIGNMENT_DIR="../../../Alignment"
GENOME_DIR="../../../genome_hg19"
GENCODE_DIR="../../../Gencode_annotation"
PBLAT_DIR="../../../pblat"

# Input files
RNA_BAM="${ALIGNMENT_DIR}/SRR1258218_Aligned.sortedByCoord.out.bam"
GENOME_FA="${GENOME_DIR}/GRCh37.primary_assembly.genome.fa"
SPLICE_SITES="${GENCODE_DIR}/gencode.v30lift37.splicesites.txt"
NOCHR_FILE="${GENOME_DIR}/nochr"

# REDItoolDnaRna.py output folder (modify as needed)
EDITING_DIR="editing_chr4_1_191154276"
DNARNA_DIR="DnaRna_39649917"
SAMPLE_ID="39649917"
CHR="chr4"

# Annotation files
RMSK_GTF="${RMSK_DIR}/rmsk.sorted.gtf.gz"
SNP_GTF="${SNP_DIR}/snp151.sorted.gtf.gz"
REDIPORTAL_GTF="${REDIPORTAL_DIR}/atlas.gtf.gz"

# Processing parameters
THREADS=4

# ============================================================================
# MAIN PIPELINE
# ============================================================================

echo "========================================================================"
echo "REDItools DNA-RNA Variant Filtering Pipeline"
echo "========================================================================"
echo "Starting: $(date)"
echo ""

# Step 1: Navigate to output folder
echo "[Step 9-10] Navigating to REDItoolDnaRna.py output folder..."
cd "${EDITING_DIR}/${DNARNA_DIR}"
##############################################################################
# Step 2: Filter invariant positions and low coverage sites
echo "[Step 11] Filtering invariant positions and sites with <10 WGS reads..."
q='"'
awk -v FS="\t" '{if ($8!="-" && $10>=10 && $13=="-") print}' "outTable_${SAMPLE_ID}" > "outTable_${SAMPLE_ID}.out"
###############################################################################
# Step 3: Annotate with RepeatMasker
echo "[Step 12] Annotating with RepeatMasker..."
python "${REDITOOLS_DIR}/accessory/AnnotateTable.py" \
    -a "${RMSK_GTF}" \
    -n rmsk \
    -i "outTable_${SAMPLE_ID}_${CHR}.out" \
    -o "outTable_${SAMPLE_ID}_${CHR}.out.rmsk" \
    -u
################################################################################
# Step 4: Annotate with dbSNP
echo "[Step 13] Annotating with dbSNP..."
python "${REDITOOLS_DIR}/accessory/AnnotateTable.py" \
    -a "${SNP_GTF}" \
    -n snp151 \
    -i "outTable_${SAMPLE_ID}_${CHR}.out.rmsk" \
    -o "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp" \
    -u
################################################################################
# Step 5: Create first set of positions (≥5 RNA reads, 1 mismatch)
echo "[Step 14] Selecting positions: ≥5 RNAseq reads, 1 mismatch..."
python "${REDITOOLS_DIR}/accessory/selectPositions.py" \
    -i "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp" \
    -c 5 -v 1 -f 0.0 \
    -o "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.sel1"
###############################################################################
# Step 6: Create second set of positions (≥10 RNA reads, 3 mismatches, freq≥0.1)
echo "[Step 15] Selecting positions: ≥10 RNAseq reads, 3 mismatches, freq≥0.1..."
python "${REDITOOLS_DIR}/accessory/selectPositions.py" \
    -i "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp" \
    -c 10 -v 3 -f 0.1 \
    -o "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.sel2"
###############################################################################
# Step 7: Select ALU sites
echo "[Step 16] Selecting ALU sites..."
awk -v FS="\t" '{if ($1!="chrM" && substr($16,1,3)=="Alu" && $17=="-" && $8!="-" && $10>=10 && $13=="-") print}' \
    "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.sel1" > "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.alu"
################################################################################
# Step 8: Select REP NON ALU sites
echo "[Step 17] Selecting REP NON ALU sites..."
awk -v FS="\t" '{if ($1!="chrM" && substr($16,1,3)!="Alu" && $15!="-" && $15!="Simple_repeat" && $15!="Low_complexity" && $17=="-" && $8!="-" && $10>=10 && $14<=0.05 && $9>=0.1) print}' \
    "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.sel2" > "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.nonalu"
################################################################################
# Step 9: Select NON REP sites
echo "[Step 09] Selecting NON REP sites..."
awk -v FS="\t" '{if ($1!="chrM" && substr($16,1,3)!="Alu" && $15=="-" && $17=="-" && $8!="-" && $10>=10 && $14<=0.05 && $9>=0.1) print}' \
    "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.sel2" > "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.nonrep"

# Step 19: Annotate with REDIportal
echo "[Step 19] Annotating ALU sites with REDIportal..."
python "${REDITOOLS_DIR}/accessory/AnnotateTable.py" \
    -a "${REDIPORTAL_GTF}" \
    -n ed -k R -c 1 \
    -i "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.alu" \
    -o "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.alu.ed" \
    -u

echo "[Step 19] Annotating REP NON ALU sites with REDIportal..."
python "${REDITOOLS_DIR}/accessory/AnnotateTable.py" \
    -a "${REDIPORTAL_GTF}" \
    -n ed -k R -c 1 \
    -i "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.nonalu" \
    -o "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.nonalu.ed" \
    -u

echo "[Step 19] Annotating NON REP sites with REDIportal..."
python "${REDITOOLS_DIR}/accessory/AnnotateTable.py" \
    -a "${REDIPORTAL_GTF}" \
    -n ed -k R -c 1 \
    -i "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.nonrep" \
    -o "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.nonrep.ed" \
    -u

# Step 20: Extract known editing events
echo "[Step 20] Extracting known editing events..."
mv "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.alu.ed" alu
mv "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.nonalu.ed" nonalu
mv "outTable_${SAMPLE_ID}_${CHR}.out.rmsk.snp.nonrep.ed" nonrep

cat alu nonalu nonrep > alu-nonalu-nonrep
awk -v FS="\t" '{if ($19=="ed") print}' alu-nonalu-nonrep > knownEditing

# Step 21: Convert REP NON ALU and NON REP novel sites to GFF
echo "[Step 21] Converting REP NON ALU and NON REP novel sites to GFF..."
cat nonalu nonrep > nonalu-nonrep
awk -v FS="\t" '{if ($19!="ed") print}' nonalu-nonrep > pos.txt

python "${REDITOOLS_DIR}/accessory/TableToGFF.py" \
    -i pos.txt -s -t -o pos.gff

# Step 22: Convert ALU novel sites to GFF
echo "[Step 22] Converting ALU novel sites to GFF..."
awk -v FS="\t" '{if ($19!="ed") print}' alu > posalu.txt

python "${REDITOOLS_DIR}/accessory/TableToGFF.py" \
    -i posalu.txt -s -t -o posalu.gff

# Step 23: Run REDItoolDnaRna.py on ALU sites with stringent criteria
echo "[Step 23] Running REDItoolDnaRna.py on ALU sites (stringent)..."
python "${REDITOOLS_DIR}/main/REDItoolDnaRna.py" \
    -s 2 -g 2 -S -t ${THREADS} \
    -i "${RNA_BAM}" \
    -f "${GENOME_FA}" \
    -c 5,5 -q 30,30 -m 255,255 -O 5,5 \
    -p -u -a 11-6 -l -v 1 -n 0.0 -e \
    -T posalu.sorted.gff.gz \
    -w "${SPLICE_SITES}" \
    -k "${NOCHR_FILE}" \
    -R -o firstalu

# Step 24: Run REDItoolDnaRna.py on REP NON ALU and NON REP sites
echo "[Step 24] Running REDItoolDnaRna.py on REP NON ALU and NON REP sites..."
python "${REDITOOLS_DIR}/main/REDItoolDnaRna.py" \
    -s 2 -g 2 -S -t ${THREADS} \
    -i "${RNA_BAM}" \
    -f "${GENOME_FA}" \
    -c 10,10 -q 30,30 -m 255,255 -O 5,5 \
    -p -u -a 11-6 -l -v 3 -n 0.1 -e \
    -T pos.sorted.gff.gz \
    -w "${SPLICE_SITES}" \
    -k "${NOCHR_FILE}" \
    --reads -R --addP -o first

# Step 25: Run pblat and identify multi-mapping reads
echo "[Step 25] Running pblat to identify multi-mapping reads..."
FIRST_DNARNA=$(find first -type d -name "DnaRna_*" | head -n1)
FIRST_ID=$(basename ${FIRST_DNARNA} | sed 's/DnaRna_//')

"${PBLAT_DIR}/pblat" \
    -t=dna -q=rna -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 \
    "${GENOME_FA}" \
    "${FIRST_DNARNA}/outReads_${FIRST_ID}" \
    reads.psl

python "${REDITOOLS_DIR}/accessory/readPsl.py" reads.psl badreads.txt

# Step 26: Extract and deduplicate reads
echo "[Step 26] Extracting and deduplicating reads..."
sort -k1,1 -k2,2n -k3,3n "${FIRST_DNARNA}/outPosReads_${FIRST_ID}" | \
    mergeBed > bed

samtools view -@ ${THREADS} -L bed -h -b "${RNA_BAM}" > SRR1258218_bed.bam
samtools sort -@ ${THREADS} -n SRR1258218_bed.bam -o SRR1258218_bed_ns.bam
samtools fixmate -@ ${THREADS} -m SRR1258218_bed_ns.bam SRR1258218_bed_ns_fx.bam
samtools sort -@ ${THREADS} SRR1258218_bed_ns_fx.bam -o SRR1258218_bed_ns_fx_st.bam
samtools markdup -r -@ ${THREADS} SRR1258218_bed_ns_fx_st.bam SRR1258218_bed_dedup.bam
samtools index SRR1258218_bed_dedup.bam

# Step 27: Re-run REDItoolDnaRna.py with deduplicated reads
echo "[Step 27] Re-running REDItoolDnaRna.py with deduplicated reads..."
python "${REDITOOLS_DIR}/main/REDItoolDnaRna.py" \
    -s 2 -g 2 -S -t ${THREADS} \
    -i SRR1258218_bed_dedup.bam \
    -f "${GENOME_FA}" \
    -c 10,10 -q 30,30 -m 255,255 -O 5,5 \
    -p -u -a 11-6 -l -v 3 -n 0.1 -e \
    -T pos.sorted.gff.gz \
    -w "${SPLICE_SITES}" \
    -R -k "${NOCHR_FILE}" \
    -b badreads.txt --rmIndels -o second

# Step 28: Collect filtered editing candidates
echo "[Step 28] Collecting filtered editing candidates..."
python "${REDITOOLS_DIR}/NPscripts/collect_editing_candidates.py"
sort -k1,1 -k2,2n editing.txt > editing_sorted.txt

# Step 29: Generate statistics
echo "[Step 29] Generating editing statistics..."
python "${REDITOOLS_DIR}/NPscripts/get_Statistics.py"

echo ""
echo "========================================================================"
echo "Pipeline completed successfully!"
echo "========================================================================"
echo "Finished: $(date)"
echo ""
echo "Output files:"
echo "  - editing_sorted.txt: Final RNA editing candidates for ${CHR}"
echo "  - editingStats.txt: Editing statistics and A-to-I enrichment"
echo "  - knownEditing: Known editing events from REDIportal"
echo ""
echo "Next steps:"
echo "  1. Review editing_sorted.txt for candidate sites"
echo "  2. Check editingStats.txt for A-to-I enrichment"
echo "  3. Validate candidates using independent methods"
echo "========================================================================"
