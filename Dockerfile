FROM trinityrnaseq/trinityrnaseq:2.8.4 

ARG SKEWER_VERSION=0.2.2
ARG TOOL_MULTIQC_VERSION=1.1
ARG TRANS_DECODER_VERSION=5.5.0
ARG HMMER_VERSION=3.2.1

# Installing Skewer
RUN bash -c 'curl -k -L https://downloads.sourceforge.net/project/skewer/Binaries/skewer-${SKEWER_VERSION}-linux-x86_64 > /usr/local/bin/skewer'
RUN chmod +x /usr/local/bin/skewer

#Adding perl script for improving multiQC report
RUN bash -c 'curl -k -L  https://github.com/CRG-CNAG/make_tool_desc_for_multiqc/archive/v${TOOL_MULTIQC_VERSION}.tar.gz > tool_ver.tar.gz'
RUN tar -zvxf tool_ver.tar.gz
RUN mv make_tool_desc_for_multiqc-${TOOL_MULTIQC_VERSION}/make_tool_desc_for_multiqc.pl /usr/local/bin/ 
RUN chmod +x /usr/local/bin/make_tool_desc_for_multiqc.pl
RUN rm -fr make_tool_desc_for_multiqc-* v${TOOL_MULTIQC_VERSION}.tar.gz 

# Installing TransDecoder
RUN bash -c 'curl -k -L https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v${TRANS_DECODER_VERSION}.tar.gz > transdec.tar.gz'
RUN tar -zvxf transdec.tar.gz; ln -s $PWD/TransDecoder-TransDecoder-v${TRANS_DECODER_VERSION}/TransDecoder.LongOrfs /usr/local/bin/; ln -s $PWD/TransDecoder-TransDecoder-v${TRANS_DECODER_VERSION}/TransDecoder.Predict /usr/local/bin/ 

# installing HMMER
RUN bash -c 'curl -k -L http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz > hmmer.tar.gz'
RUN tar -zvxf hmmer.tar.gz; cd hmmer-${HMMER_VERSION}; ./configure; make; make install; cd easel; make install; cd ../../

#cleaning
RUN rm -fr *.tar.gz* .tar.bz2
