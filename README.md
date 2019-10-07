# ![transcriptome_assembly](https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png) De Novo Transcriptome Assembly
Biocore's de novo transcriptome assembly workflow based on Nextflow

[![DOI](https://zenodo.org/badge/171497634.svg)](https://zenodo.org/badge/latestdoi/171497634)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Nextflow version](https://img.shields.io/badge/nextflow-%E2%89%A50.31.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker Build Status](https://img.shields.io/docker/automated/biocorecrg/trinity_assembly.svg)](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/trinity_assembly)

## Installation

*sh INSTALL.sh* 
it will check the presence of Nextflow in your path, the presence of singularity and will download the BioNextflow library and information about the tools used. 

You need either **Singularity** or **Docker** to launch the pipeline.

## Module denovo_assembly

This module allows to perform de novo assembly and to retrieve both predicted transcripts and proteins.

```bash
╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╔╦╗┬─┐┌─┐┌┐┌┌─┐┌─┐┬─┐┬┌─┐┌┬┐┌─┐┌┬┐┌─┐  ╔═╗┌─┐┌─┐┌─┐┌┬┐┌┐ ┬ ┬ ┬
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦   ║ ├┬┘├─┤│││└─┐│  ├┬┘│├─┘ │ │ ││││├┤   ╠═╣└─┐└─┐├┤ │││├┴┐│ └┬┘
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝   ╩ ┴└─┴ ┴┘└┘└─┘└─┘┴└─┴┴   ┴ └─┘┴ ┴└─┘  ╩ ╩└─┘└─┘└─┘┴ ┴└─┘┴─┘┴ 
                                                                                
====================================================
BIOCORE@CRG Transcriptome Assembly - N F  ~  version 0.1
====================================================
pairs                               : ../test_data/*_{1,2}.fq.gz
email                               : YOUREMAIL@YOURDOMAIN
minsize (after filtering)           : 70
genetic code                        : Universal
strandness                          : RF
output (output folder)              : output
minProtSize (minimum protein sized) : 100
```

## Module RABT_assembly

This module allows to perform de reference annotation based transcript (RABT) assembly and to retrieve both predicted transcripts and proteins.

```bash
╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╔╦╗┬─┐┌─┐┌┐┌┌─┐┌─┐┬─┐┬┌─┐┌┬┐┌─┐┌┬┐┌─┐  ╔═╗┌─┐┌─┐┌─┐┌┬┐┌┐ ┬ ┬ ┬
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦   ║ ├┬┘├─┤│││└─┐│  ├┬┘│├─┘ │ │ ││││├┤   ╠═╣└─┐└─┐├┤ │││├┴┐│ └┬┘
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝   ╩ ┴└─┴ ┴┘└┘└─┘└─┘┴└─┴┴   ┴ └─┘┴ ┴└─┘  ╩ ╩└─┘└─┘└─┘┴ ┴└─┘┴─┘┴ 
                                                                                
====================================================
BIOCORE@CRG Transcriptome Assembly - N F  ~  version 0.1
====================================================
pairs                               : ../test_data2/*_{1,2}.fq.gz
genome                              : ../anno/GRCh38.p12.genome.fa.g
z
annotation                          : ../anno/gencode.v30.annotation.gtf
minsize (after filtering)           : 40
genetic code                        : Universal
output (output folder)              : output
minProtSize (minimum protein sized) : 100
strandness                          : RF
maxIntron                           : 10000
email                               : YOUREMAIL@YOURDOMAIN

```

## Module annotation

This module allows to annotate predicted proteins and transcripts from one of the two assembly modules described before.
```bash
╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╔╦╗┬─┐┌─┐┌┐┌┌─┐┌─┐┬─┐┬┌─┐┌┬┐┌─┐┌┬┐┌─┐  ╔═╗┌─┐┌─┐┌─┐┌┬┐┌┐ ┬ ┬ ┬
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦   ║ ├┬┘├─┤│││└─┐│  ├┬┘│├─┘ │ │ ││││├┤   ╠═╣└─┐└─┐├┤ │││├┴┐│ └┬┘
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝   ╩ ┴└─┴ ┴┘└┘└─┘└─┘┴└─┴┴   ┴ └─┘┴ ┴└─┘  ╩ ╩└─┘└─┘└─┘┴ ┴└─┘┴─┘┴ 
                                                                                
====================================================
BIOCORE@CRG Transcriptome Annotation - N F  ~  version 0.1
====================================================
peptide sequences                   : ../assembly/output/Assembly/lon
gest_orfs.pep
cds sequences                       : ../assembly/output/Assembly/lon
gest_orfs.cds
annotation in gff3                  : ../assembly/output/Assembly/longest_orfs.gff3
transcripts                         : ../assembly/output/Assembly/Trinity.fasta
email                               : YOUREMAIL@YOURDOMAIN
genetic code                        : Universal
output (output folder)              : output
diamondDB (uniprot or uniRef90)     : /nfs/db/uniprot/2018_10/knowledgebase/complete/blast/db/uniprot_sprot.fasta
pfamDB (pfam database path)         : /nfs/db/pfam/Pfam31.0/Pfam-A.hmm
minProtSize (minimum protein sized) : 100
batch_diam                          : 5000
batch_pfam                          : 2000

```

## Module quantify

This module allows to quantify predicted genes obtained from one of the two assembly modules described before.

```bash
╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╔╦╗┬─┐┌─┐┌┐┌┌─┐┌─┐┬─┐┬┌─┐┌┬┐┌─┐┌┬┐┌─┐  ╔═╗┌─┐┌─┐┌─┐┌┬┐┌┐ ┬ ┬ ┬
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦   ║ ├┬┘├─┤│││└─┐│  ├┬┘│├─┘ │ │ ││││├┤   ╠═╣└─┐└─┐├┤ │││├┴┐│ └┬┘
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝   ╩ ┴└─┴ ┴┘└┘└─┘└─┘┴└─┴┴   ┴ └─┘┴ ┴└─┘  ╩ ╩└─┘└─┘└─┘┴ ┴└─┘┴─┘┴ 
                                                                                
====================================================
BIOCORE@CRG Transcriptome Quantification - N F  ~  version 0.1
====================================================
pairs                               : ../test_data/*_{1,2}.fq.gz
transcripts                         : ../assembly/output/Assembly/Trinity.fasta
transmap                            : ../assembly/output/Assembly/Trinity.fasta.gene_trans_map
output                              : output
email                               : YOUREMAIL@YOURDOMAIN

```

