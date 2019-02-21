FROM trinityrnaseq/trinityrnaseq:2.8.4 

ARG SKEWER_VERSION=0.2.2
ARG MULTIQC_VERSION=1.7
ARG TOOL_MULTIQC_VERSION=1.1
ARG TRANS_DECODER_VERSION=5.5.0

# Installing Skewer
RUN bash -c 'curl -k -L https://downloads.sourceforge.net/project/skewer/Binaries/skewer-${SKEWER_VERSION}-linux-x86_64 > /usr/local/bin/skewer'
RUN chmod +x /usr/local/bin/skewer

# Installing MULTIQC
RUN pip install -Iv https://github.com/ewels/MultiQC/archive/v${MULTIQC_VERSION}.tar.gz 

#Adding perl script for improving multiQC report
RUN bash -c 'curl -k -L  https://github.com/CRG-CNAG/make_tool_desc_for_multiqc/archive/v${TOOL_MULTIQC_VERSION}.tar.gz > tool_ver.tar.gz'
RUN tar -zvxf tool_ver.tar.gz
RUN mv make_tool_desc_for_multiqc-${TOOL_MULTIQC_VERSION}/make_tool_desc_for_multiqc.pl /usr/local/bin/ 
RUN chmod +x /usr/local/bin/make_tool_desc_for_multiqc.pl
RUN rm -fr make_tool_desc_for_multiqc-* v${TOOL_MULTIQC_VERSION}.tar.gz 

# Installing TransDecoder

RUN bash -c 'curl -k -L https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v${TRANS_DECODER_VERSION}.tar.gz > transdec.tar.gz'
RUN tar -zvxf transdec.tar.gz; cp TransDecoder-TransDecoder-v${TRANS_DECODER_VERSION}/TransDecoder.LongOrfs /usr/local/bin/; cp TransDecoder-TransDecoder-v${TRANS_DECODER_VERSION}/TransDecoder.Predict /usr/local/bin/ 

#cleaning
RUN yum clean all
