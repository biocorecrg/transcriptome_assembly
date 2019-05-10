#!/usr/bin/env nextflow


/* 
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '0.1'

params.help            = false
params.resume          = false

log.info """

╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╔╦╗┬─┐┌─┐┌┐┌┌─┐┌─┐┬─┐┬┌─┐┌┬┐┌─┐┌┬┐┌─┐  ╔═╗┌─┐┌─┐┌─┐┌┬┐┌┐ ┬ ┬ ┬
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦   ║ ├┬┘├─┤│││└─┐│  ├┬┘│├─┘ │ │ ││││├┤   ╠═╣└─┐└─┐├┤ │││├┴┐│ └┬┘
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝   ╩ ┴└─┴ ┴┘└┘└─┘└─┘┴└─┴┴   ┴ └─┘┴ ┴└─┘  ╩ ╩└─┘└─┘└─┘┴ ┴└─┘┴─┘┴ 
                                                                                       
====================================================
BIOCORE@CRG Transcriptome Quantification - N F  ~  version ${version}
====================================================
pairs                               : ${params.pairs}
transcripts                         : ${params.transcripts}
transmap                            : ${params.transmap}
output                              : ${params.output}
email                               : ${params.email}
"""

if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

util_scripts_image_path = "/usr/local/bin/trinityrnaseq/util/"
support_scripts_image_path = "${util_scripts_image_path}/support_scripts"

/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.pairs )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.pairs}" }  
    .set { read_pairs }

transcripts = file(params.transcripts)
transmap = file(params.transmap)
outputMatrix = ${params.output}/"Expression"
outputEvaluation = ${params.output}/"Evaluation"

/*
*  Prepare the reference for alignment and abundance estimation
*/
process prepReference {
    label "big_mem_cpus"

    input:
    file(transcripts)
    file(transmap)

    output:
    file ("${transcripts}.salmon_quasi.idx") into indexed_transcripts

    script:
    """
    ${util_scripts_image_path}/align_and_estimate_abundance.pl --transcripts ${transcripts} --est_method salmon --gene_trans_map ${transmap} --prep_reference
    """
}

/*
*  Align the reads
*/
process alignReadsToTranscritps {
    tag("$pair_id")
    label "big_cpus"

    input:
    file(transcripts)
    file(indexed_transcripts)
    set pair_id, file(pairs) from read_pairs
    
    output: 
    file("${pair_id}") into abund_estimates
         
    script:
    """
    ${util_scripts_image_path}/align_and_estimate_abundance.pl --thread_count ${task.cpus} --transcripts ${transcripts} --seqType fq --left ${pairs[0]} --right ${pairs[1]} --est_method salmon --trinity_mode --output_dir ${pair_id}
    """
}

/*
*  Build the matrix
*/
process BuildMatrices {
    publishDir outputMatrix, mode: 'copy'

    input:
    file(transcripts)
    file(transmap)
    file(estimates) from abund_estimates.collect()
    
    output: 
    file("salmon.gene.TPM.not_cross_norm") into tpm_matrix_for_evaluation
         
    script:
    def estimates_params = estimates.collect{ "$it/quant.sf" }.join(' ') 

    """
    ${util_scripts_image_path}/abundance_estimates_to_matrix.pl \
            --est_method salmon \
            --gene_trans_map ${transmap} \
            --name_sample_by_basedir \
            ${estimates_params}
    """
}

process evaluateTranscription {
    publishDir outputEvaluation, mode: 'copy'

    input:
    file(tpm_matrix_for_evaluation)

    output: 
    file("${tpm_matrix_for_evaluation}.counts_by_min_TPM*") 
         
    script:
    """
    ${util_scripts_image_path}/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
		${tpm_matrix_for_evaluation} | tee ${tpm_matrix_for_evaluation}.counts_by_min_TPM
	
	plot_TPM.r  ${tpm_matrix_for_evaluation}.counts_by_min_TPM
    """
}



workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

