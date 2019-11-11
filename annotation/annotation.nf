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
BIOCORE@CRG Transcriptome Annotation - N F  ~  version ${version}
====================================================
peptide sequences                   : ${params.peps}
cds sequences                       : ${params.cdss}
annotation in gff3                  : ${params.gff3}
transcripts                         : ${params.transcripts}
email                               : ${params.email}
genetic code                        : ${params.geneticode}
output (output folder)              : ${params.output}
diamondDB (uniprot or uniRef90)     : ${params.diamondDB}
pfamDB (pfam database path)         : ${params.pfamDB}
minProtSize (minimum protein sized) : ${params.minProtSize}
batch_diam                          : ${params.batch_diam}
batch_pfam                          : ${params.batch_pfam}
"""

if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

minContigSize   = (params.minProtSize*3)
outputfolder    = "${params.output}"
outputAnnotation= "${outputfolder}/Annotation"

util_scripts_image_path = "/usr/local/src/TransDecoder-TransDecoder-v5.5.0/util/"

/*
 * Creates the `peptides` channel 
 */
 
Channel
    .fromPath( params.peps )                                             
    .ifEmpty { error "Cannot find any peptide file: ${params.peps}" }  
    .into { peptides_for_diamond; peptides_for_pfam; peptides_for_prediction }


Channel
    .fromPath( params.gff3 )                                             
    .ifEmpty { error "Cannot find any gff3 annotation file: ${params.gff3}" }  
    .set { gff3_for_prediction }

Channel
    .fromPath( params.cdss )                                             
    .ifEmpty { error "Cannot find any cds file: ${params.cdss}" }  
    .set { cdss_for_prediction }

Channel
    .fromPath( params.transcripts )                                             
    .ifEmpty { error "Cannot find any transcript file: ${params.transcripts}" }  
    .set { transcripts_for_prediction }

process diamondSearch {  
    tag "$pep_batches"  
    label("big_cpus")

    input:
    file(pep_batches) from peptides_for_diamond.splitFasta( by: params.batch_diam, file: true )

    output:
    file ("blastp.outfmt6") into blastout
    
    // fixing the --max_target_seqs problem
    script:
    """
    diamond blastp --sensitive -d ${params.diamondDB} -q ${pep_batches} -p ${task.cpus} > diamond.out
    awk 'BEGIN{print \$0;id=\$1}{if (id!=\$1){print \$0; id=\$1} }' diamond.out | awk '{if (\$0!="") print}' > blastp.outfmt6
    rm diamond.out
    """
}

process pfam_search {  
    tag "$pep_batches"      
    when params.pfamDB != ""
    label("big_cpus")
    
    input:
    file(pep_batches) from peptides_for_pfam.splitFasta( by: params.batch_pfam, file: true )

    output:
    file ("pfam.domtblout") into pfamout
    
    script:
    """
    hmmscan --cpu ${task.cpus} --domtblout pfam.domtblout ${params.pfamDB} ${pep_batches}
    """
}

process concatenateBlastRes {    
    publishDir outputAnnotation, mode: 'copy'
    input:
    file("blastres*") from blastout.collect()

    output:
    file ("blastp.all.results") into blastoutall

    script:
    """
    cat blastres* >> blastp.all.results
    """
}

process concatenatePfamRes {    
    publishDir outputAnnotation, mode: 'copy'

    when params.pfamDB != ""

    input:
    file("pfamout*") from pfamout.collect()

    output:
    file ("pfamout.all.results") into pfamoutres

    script:
    """
    cat pfamout* >> pfamout.all.results
    """
}

if (params.pfamDB == "") {
    pfamoutall = file("EMPTY_FILE")
} else {
    pfamoutall = pfamoutres
}
/*
*/
process transcoderPredict {    
    publishDir outputAnnotation, mode: 'copy'
    label("big_mem")

    input:
    file(blastoutall)
    file (pfam) from pfamoutall
    file(transcripts_for_prediction)
    file(peptides_for_prediction)
    file(gff3_for_prediction)
    file(cdss_for_prediction)

    output:
    set file("*.transdecoder.bed"), file("*.transdecoder.cds"), file("*.transdecoder.gff3"), file("*.transdecoder.pep") into finalRes

    script:
    def filter = pfam.name != 'EMPTY_FILE' ? "--retain_pfam_hits ${pfam}" : ''

    """
    mkdir ${transcripts_for_prediction}.transdecoder_dir
    ln -s \$PWD/`basename ${gff3_for_prediction}` ${transcripts_for_prediction}.transdecoder_dir/
    ln -s \$PWD/`basename ${cdss_for_prediction}` ${transcripts_for_prediction}.transdecoder_dir/
    ln -s \$PWD/`basename ${peptides_for_prediction}` ${transcripts_for_prediction}.transdecoder_dir/
    ${util_scripts_image_path}/compute_base_probs.pl ${transcripts_for_prediction} 1 > Trinity.fasta.transdecoder_dir/base_freqs.dat;
    TransDecoder.Predict --no_refine_starts -G ${params.geneticode} -t ${transcripts_for_prediction} ${filter} --retain_blastp_hits ${blastoutall}
    """
}




workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

/*
* send mail

workflow.onComplete {
    def subject = 'Transcriptome assembly execution'
    def recipient = "${params.email}"
    def attachment = "${outputMultiQC}/multiqc_report.html"

    ['mail', '-s', subject, '-a', attachment, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}



*/


