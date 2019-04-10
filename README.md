# ![transcriptome_assembly](https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png) De Novo Transcriptome 
Biocore's de novo transcriptome assembly workflow based on Nextflow


[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Nextflow version](https://img.shields.io/badge/nextflow-%E2%89%A50.31.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker Build Status](https://img.shields.io/docker/automated/biocorecrg/trinity_assembly.svg)](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/trinity_assembly)

## Installation
*sh INSTALL.sh* 
it will check the presence of Nextflow in your path, the presence of singularity and will download the BioNextflow library and information about the tools used. 

You need either **Singularity** or **Docker** to launch the pipeline.

## Module Assembly

### Parameters

|parameter name         | value|
|---------------------------------|------------------------|
|pairs         | "PATH_to_reads/*_{1,2}.fq.gz"|
|email         | YOUREMAIL@YOURDOMAIN|
|minsize|50|
|genetic code|Universal|
|output (output folder)|output|
|minProtSize (minimum protein sized)|100|
