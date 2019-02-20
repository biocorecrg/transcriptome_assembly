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
BIOCORE@CRG Transcriptome Assembly - N F  ~  version ${version}
====================================================
pairs                         : ${params.pairs}
email                         : ${params.email}
minsize (after filtering)     : ${params.minsize}
orientation (RF or FR)        : ${params.orientation}
output (output folder)        : ${params.output}
"""

if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"
if (params.orientation != "RF" && params.orientation != "FR")  exit 1, "Orientation not allowed. Only FR or RF" 


outputfolder    = "${params.output}"
outputQC        = "${outputfolder}/QC"
outputTrimmed   = "${outputfolder}/trimmedReads"
outputMultiQC   = "${outputfolder}/multiQC"
outputMapping   = "${outputfolder}/Alignments"
est_folder      = "${outputfolder}/Estimated_counts"
rep_folder      = "${outputfolder}/Reports"


/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.pairs )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.pairs}" }  
    .into { read_pairs; raw_reads_for_trimming; raw_reads_for_assembly }

Channel
    .fromPath( params.pairs )
    .set{ reads_for_fastqc }                                           


/*
 * Run FastQC on raw data
*/
process QConRawReads {
    publishDir outputQC
    tag { read }

    input:
    file(read) from reads_for_fastqc

    output:
    file("*_fastqc.*") into raw_fastqc_files

    script:
    def qc = new QualityChecker(input:read, cpus:task.cpus)
    qc.fastqc()
}

 /*
 * Trim reads with skewer (optional). Make empty channels for multiQC
 */
process trimReads {
    publishDir outputTrimmed
    tag { pair_id }
    label 'rnaseq_cpus'

    input:
    set pair_id, file(reads) from (raw_reads_for_trimming )

    output:
    set val("pair1"), file("*-trimmed-pair1.fastq.gz") into trimmed_pair1_for_assembly
    set val("pair2"), file("*-trimmed-pair2.fastq.gz") into trimmed_pair2_for_assembly
    file("*trimmed*.fastq.gz") into filtered_read_for_QC
    file("*trimmed.log") into logTrimming_for_QC
     
    script:
    def trimmer = new Trimmer(reads:reads, id:pair_id, min_read_size:params.minsize, cpus:task.cpus)
    trimmer.trimWithSkewer()
}

process fastqcTrim {
    tag { filtered_read }
    publishDir outputQC

 	afterScript 'mv *_fastqc.zip `basename *_fastqc.zip _fastqc.zip`_sel_fastqc.zip'

    input:
    file(filtered_read) from filtered_read_for_QC.flatten()

    output:
    file("*_sel_fastqc.zip") into trimmed_fastqc_files

    script:
    def qc = new QualityChecker(input:filtered_read, cpus:task.cpus)
    qc.fastqc()
}

//trimmed_pair_for_assembly = trimmed_pair1_for_assembly.mix(trimmed_pair2_for_assembly).groupTuple().collect()


process Inchworm {

    input:
    set val(pair1), file(pair1) from trimmed_pair1_for_assembly.groupTuple()
    set val(pair2), file(pair2) from trimmed_pair2_for_assembly.groupTuple()


    script:
    def pair1_list = pair1.join(',')
    def pair2_list = pair2.join(',')
    
    """
    echo ${pair1_list} ${pair2_list}
    """
}


/*
* send mail

workflow.onComplete {
    def subject = 'indropSeq execution'
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



