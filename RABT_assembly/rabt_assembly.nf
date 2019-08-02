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
BIOCORE@CRG RABT Transcriptome Assembly - N F  ~  version ${version}
====================================================
pairs                               : ${params.pairs}
genome                              : ${params.genome}
annotation                          : ${params.annotation}
minsize (after filtering)           : ${params.minsize}
splitsize (for transdecoder)        : ${params.splitsize}
genetic code                        : ${params.geneticode}
output (output folder)              : ${params.output}
minProtSize (minimum protein sized) : ${params.minProtSize}
strandness                          : ${params.strandness}
maxIntron                           : ${params.maxIntron}
email                               : ${params.email}
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
outputIndex     = "${outputfolder}/Index"
outputCounts    = "${outputfolder}/Counts"

util_scripts_image_path = "/usr/local/bin/trinityrnaseq/util/"
support_scripts_image_path = "${util_scripts_image_path}/support_scripts"

genome_file = file(params.genome)

if( !genome_file.exists() ) exit 1, "Missing genome file: ${genome_file}"
if( params.annotation != "") {
	annotation_file = file(params.annotation)
	if( !annotation_file.exists() ) exit 1, "Missing annotation file: ${annotation_file}"
} else {
	annotation_file = ""
   log.info """
   Executing alignment with NO annotation!!
   """
}
if (params.strandness != "FR" && params.strandness != "RF" && params.strandness != "NO" ) exit 1, "Please specify FR , RF or NO in case the data is in stranded or non stranded"

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

    input:
    set pair_id, file(reads) from (raw_reads_for_trimming )

    output:
    set val("pair1"), file("*-trimmed-pair1.fastq.gz") into trimmed_pair1_for_normalization
    set val("pair2"), file("*-trimmed-pair2.fastq.gz") into trimmed_pair2_for_normalization
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

/*
*/
process TrinityNormalization {
    label 'big_time_cpus'
    
    input:
    set val(pair1), file(pair1) from trimmed_pair1_for_normalization.groupTuple()
    set val(pair2), file(pair2) from trimmed_pair2_for_normalization.groupTuple()

    output:
    file "trinity_out_dir/insilico_read_normalization/*.norm.fq" into norm_reads_first, norm_reads_second

    script:
    def pair1_list = pair1.join(',')
    def pair2_list = pair2.join(',')

    """     
    Trinity --seqType fq --max_memory ${task.memory.giga}G --left ${pair1_list} --right ${pair2_list} --CPU ${task.cpus} --no_run_inchworm
    """
}

process buildIndex {
    publishDir outputIndex
    label 'big_time_cpus'
    
    input:
    file genome_file

    output:
    file "STARgenome" into STARgenomeIndex, STARgenomeIndexForCoverage

    script:
    aligner = new NGSaligner(reference_file:genome_file, index:"STARgenome", annotation_file:annotation_file, read_size:params.minsize-1, cpus:task.cpus)
    aligner.doIndexing("STAR")
}


/*
* Mapping with STAR mapper (first pass)
*/
process firstPassMapping {
    label 'big_mem_cpus'
    
        afterScript 'rm STAR_*/*.bam' 

        input:
        file STARgenome from STARgenomeIndex
        file(reads) from norm_reads_first

        output:
        file "STAR_first/firstSJ.out.tab" into first_pass_junctions

        script:
        def aligner = new NGSaligner(id:"first", reads:reads, index:STARgenome, cpus:task.cpus, output:"STAR_first") 
        aligner.doAlignment("STAR")  
        
}

/*
* Mapping with STAR mapper (second pass)
*/

process secondPassMapping {
    label 'big_mem_cpus'
    publishDir outputCounts, pattern: "STAR_norm/*ReadsPerGene.out.tab",  mode: 'copy'
    publishDir outputQC, pattern: "STAR_norm/*Log.final.out", mode: 'copy'

        input:
        file STARgenome from STARgenomeIndex
        file(reads) from norm_reads_second
        file(first_pass_junctions)
    
        output:
        file "STAR_norm/normAligned.sortedByCoord.out.bam" into STARmappedBam_for_qualimap, STARmappedBam_for_assembly
        file "STAR_norm"  into Aln_folders_for_multiqc
        file "STAR_norm/normReadsPerGene.out.tab" 

        script:
        def aligner = new NGSaligner(id:"norm", reads:reads, index:STARgenome, cpus:task.cpus, output:"STAR_norm", extrapars:"--sjdbFileChrStartEnd ${first_pass_junctions}") 
        aligner.doAlignment("STAR")  
        
}

/*
* Trinity assembly step one
*/
process TrinityAssemblyStep1 {
    label 'assembly'
    
    input:
    file(STARmappedBam_for_assembly)

    output:
    file ("trinity_out_dir/*/*/*") into partitions_groups
    
    script:
    def strand = ""
    if (params.strandness != "NO") {
    	strand = "--SS_lib_type ${params.strandness}"
    }    
    """
    Trinity --genome_guided_bam ${STARmappedBam_for_assembly} \
         --genome_guided_max_intron ${params.maxIntron} ${strand} \
         --min_contig_length ${minContigSize} \
         --no_distributed_trinity_exec \
         --max_memory ${task.memory.giga}G --CPU  ${task.cpus} 
    """
}

/*
* Trinity assembly step two
*/
process TrinityStep2 {
    label 'increase_mem'
    tag { partitions_group.target }
    
    input:
    file(partitions_group) from partitions_groups.flatten()

    output:
    file ("*inity.fasta") optional true into out_trinity

    script:
    def strand = ""
    if (params.strandness != "FR") {
    	strand = "--SS_lib_type F"
    }    
    """
    if [ \$(ls ${partitions_group}/*.trinity.reads | wc -l) -gt 0 ];
    then
        for i in ${partitions_group}/*.trinity.reads; do \
        Trinity --single \$i --min_contig_length ${minContigSize} ${strand} --output ./`basename \$i`.out --CPU 1 --max_memory ${task.memory.giga}G  --seqType fa --trinity_complete --full_cleanup --no_distributed_trinity_exec; done;
    fi
    """
}

/*
* collect trinity results
*/
process collectTrinityRes {
    publishDir outputMultiQC, mode: 'copy', pattern: "Trinity.fasta.stat"
    publishDir outputAssembly, mode: 'copy', pattern: "Trinity.fasta*"
    
    input:
    file("Trinity_sub") from out_trinity.collect()

    output:
    file ("Trinity.fasta") into fasta_to_split
    file ("Trinity.fasta*") 
    
    script:
    """
    find * -name 'Trinity_sub*'  | ${support_scripts_image_path}/GG_partitioned_trinity_aggregator.pl TRINITY > Trinity.fasta;
    ${support_scripts_image_path}/get_Trinity_gene_to_trans_map.pl Trinity.fasta >  Trinity.fasta.gene_trans_map
    ${util_scripts_image_path}/TrinityStats.pl Trinity.fasta > Trinity.fasta.stat
    """
}

/*
*/

process splitTrinityFasta {
    tag { fasta_to_split }
    
    input:
    file(fasta_to_split)

    output:
    file ("pieces_*.fa") into components_for_transcoder
    
    script:
    """
    splitTrinityFasta.sh ${params.splitsize}
    """


}

/*
*/

process TransDecoder {
    tag { components }
    
    input:
    file(components) from components_for_transcoder.flatten()

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

/*
*/
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


