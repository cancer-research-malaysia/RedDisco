// Reditools v1 wrangling
process WRANGLE_FINAL_OUTPUTS_REDITOOLS_V1 {
    tag "$sampleId"
    label 'finalWranglingReditools'
    container "${params.container__reditools}"
    publishDir "${params.outputDir}/reditools-v1/${sampleId}/final", mode: 'copy'

    input:
    tuple val(sampleId), path(known), path(pos), path(pos_alu)
    tuple val(_sampleId), path(firstalu_dir), path(second_dir)

    output:
    path("${sampleId}-collected-editing-candidates--sorted.txt"), emit: final_table
    path("${sampleId}-collected-editing-candidates--sorted--editingStats.txt")
    path("${sampleId}-collected-editing-candidates--sorted--AG-sub-only.txt")
    
    script:
    """
    # Step 1: Collect filtered ALU, REP NON ALU and NON REP sites
    echo "[Step 1] Collecting filtered editing candidates..."

    if collect_editing_candidates-v2.py --prefix "${sampleId}" --known "${sampleId}-knownEditing" --pos "${sampleId}--pos.txt" --posalu "${sampleId}--pos-ALU.txt" --output "${sampleId}-collected-editing-candidates.txt" && \
    sort -k1,1 -k2,2n "${sampleId}-collected-editing-candidates.txt" > "${sampleId}-collected-editing-candidates--sorted.txt"; then
      echo "Collected editing candidates successfully."
    else
      echo "Error collecting editing candidates. Exiting."
      exit 1
    fi
    # Step 2: Generate statistics
    echo "[Step 2] Generating editing statistics..."
    if get_statistics-v2.py --input "${sampleId}-collected-editing-candidates--sorted.txt" --output "${sampleId}-collected-editing-candidates--sorted--editingStats.txt" && \
    filter_A-to-I_sites_only.py --input "${sampleId}-collected-editing-candidates--sorted.txt" --output "${sampleId}-collected-editing-candidates--sorted--AG-sub-only.txt"; then
      echo "Generated statistics successfully."
    else
      echo "Error generating statistics. Exiting."
      exit 1
    fi


    """
}
