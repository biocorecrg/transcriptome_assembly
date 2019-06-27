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
BIOCORE@CRG De Novo Transcriptome Assembly - N F  ~  version ${version}
====================================================
pairs                               : ${params.pairs}
email                               : ${params.email}
minsize (after filtering)           : ${params.minsize}
genetic code                        : ${params.geneticode}
strandness                          : ${params.strandness}
output (output folder)              : ${params.output}
minProtSize (minimum protein sized) : ${params.minProtSize}
"""

if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

minContigSize   = (params.minProtSize*3)
outputfolder    = "${params.output}"
outputQC        = "${outputfolder}/QC"
outputTrimmed   = "${outputfolder}/trimmedReads"
outputMultiQC   = "${outputfolder}/multiQC"
outputMapping   = "${outputfolder}/Alignments"
outputAssembly  = "${outputfolder}/Assembly"
outputAnnotation= "${outputfolder}/Annotation"

util_scripts_image_path = "/usr/local/bin/trinityrnaseq/util/"
support_scripts_image_path = "${util_scripts_image_path}/support_scripts"


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

if (params.strandness != "FR" && params.strandness != "RF" && params.strandness != "NO" ) exit 1, "Please specify FR , RF or NO in case the data is in stranded or non stranded"


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


process TrinityStep1 {
    label 'big_mem_cpus'
    
    input:
    set val(pair1), file(pair1) from trimmed_pair1_for_assembly.groupTuple()
    set val(pair2), file(pair2) from trimmed_pair2_for_assembly.groupTuple()

    output:
    file ("trinity_out_dir/read_partitions/*/*") into partitions_groups
    
    script:
    def pair1_list = pair1.join(',')
    def pair2_list = pair2.join(',')
    def strand = ""
    if (params.strandness != "NO") {
    	strand = "--SS_lib_type ${params.strandness}"
    }
    
    """
    Trinity --min_contig_length ${minContigSize} ${strand} \
    --seqType fq --max_memory ${task.memory.giga}G \
    --left ${pair1_list} --right ${pair2_list} --CPU ${task.cpus} --no_distributed_trinity_exec
    """
}


process TrinityStep2 {
    label 'increase_mem'
    tag { partitions_group }
    
    input:
    file(partitions_group) from partitions_groups.flatten()

    output:
    file ("Trinity_sub.fasta") into components
    set val("${partitions_group}"), file ("Trinity_sub.fasta") into components_for_transcoder
    
    script:
    def strand = ""
    if (params.strandness != "FR") {
    	strand = "--SS_lib_type F"
    }    
    """
    ls ${partitions_group} -l | awk -F'/' '{print "mkdir " \$(NF-1)}' | sh;
    ls ${partitions_group} -l | awk -F'/' '{print "mkdir " \$(NF-1)"/"\$(NF)}' | sh;
    OUTFOLDER=`ls ${partitions_group} -l | awk -F"/" '{print \$(NF-1)"/"\$(NF)"/"}'`;
    for i in ${partitions_group}/*.fa; do \
    Trinity --single \$i --min_contig_length ${minContigSize} --output \$OUTFOLDER`basename \$i`.out ${strand} --CPU 1 --max_memory ${task.memory.giga}G --run_as_paired --seqType fa --trinity_complete --full_cleanup --no_distributed_trinity_exec; done;
    find \$OUTFOLDER -name '*inity.fasta' | ${support_scripts_image_path}/partitioned_trinity_aggregator.pl --token_prefix TRINITY_DN --output_prefix Trinity_sub
    """
}


/*
*/
process TransDecoder {
    tag { partitions_group }
    
    input:
    set val(partitions_group), file(components) from components_for_transcoder

    output:
    file ("longest_orfs.pep") into orfs_for_concatenation
    file ("longest_orfs.cds") into cdss_for_concatenation
    file ("longest_orfs.gff3") into gff3_for_concatenation
    
    script:
    """
	TransDecoder.LongOrfs -m ${params.minProtSize} -t ${components} -G ${params.geneticode} -S 
	cp ${components}.transdecoder_dir/longest* .
    """
}

process collectTrinityRes {
    publishDir outputMultiQC, mode: 'copy', pattern: "Trinity.fasta.stat"
    publishDir outputAssembly, mode: 'copy', pattern: "Trinity.fasta*"
    
    input:
    file("Trinity_sub") from components.collect()

    output:
    file ("Trinity.fasta*")
    
    script:
    """
    cat Trinity_sub* >> Trinity.fasta
    ${support_scripts_image_path}/get_Trinity_gene_to_trans_map.pl Trinity.fasta >  Trinity.fasta.gene_trans_map
    ${util_scripts_image_path}/TrinityStats.pl Trinity.fasta > Trinity.fasta.stat
    """
}

process collectTransDecoderRes {
    publishDir outputAssembly, mode: 'copy', pattern: "longest*"
    
    input:
    file("peps_sub") from orfs_for_concatenation.collect()
    file("cds_sub") from cdss_for_concatenation.collect()
    file("gff3_sub") from gff3_for_concatenation.collect()

    output:
    file ("longest_orfs.*")

    script:
    """
    cat peps_sub* >> longest_orfs.pep
    cat cds_sub* >> longest_orfs.cds
    cat gff3_sub* >> longest_orfs.gff3
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


