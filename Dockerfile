FROM centos:centos7

RUN yum install -y epel-release

RUN yum -y install \
	git \
	wget \
	java-1.8.0-openjdk \
	java-1.8.0-openjdk-devel \
	R \
	autoconf \
	automake \
	make \
	gcc \
	perl-Data-Dumper \
	zlib-devel \
	bzip2 \
	bzip2-devel \
	xz-devel \
	curl-devel \
	openssl-devel \
	ncurses-devel \
	graphviz


ENV APPS_ROOT /apps
RUN mkdir -p ${APPS_ROOT}

###############################################
#BWA = 'bwa/intel/0.7.17'

ENV BWA_VERSION 0.7.17

ENV BWA_HOME ${APPS_ROOT}/bwa/${BWA_VERSION}
ENV PATH ${BWA_HOME}:${PATH}

RUN git clone --branch v${BWA_VERSION} https://github.com/lh3/bwa.git ${BWA_HOME}
RUN cd ${BWA_HOME} && make && cd

###############################################
#PICARD = 'picard/2.17.11'

ENV PICARD_VERSION 2.17.11

ENV JAVA_HOME /etc/alternatives/jre
ENV PICARD_HOME ${APPS_ROOT}/picard/${PICARD_VERSION}
ENV PICARD_JAR ${PICARD_HOME}/picard-${PICARD_VERSION}.jar

RUN mkdir -p ${PICARD_HOME}
RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar -O ${PICARD_JAR}

###############################################
#GATK = 'gatk/4.1.3.0'

ENV GATK_VERSION 4.1.3.0

ENV GATK_HOME ${APPS_ROOT}/gatk/${GATK_VERSION}

ENV GATK_LOCAL_JAR ${GATK_HOME}/gatk-package-${GATK_VERSION}-local.jar
ENV GATK_SPARK_JAR ${GATK_HOME}/gatk-package-${GATK_VERSION}-spark.jar
ENV GATK_JAR ${GATK_HOME}/gatk-package-${GATK_VERSION}-local.jar
ENV PATH ${GATK_HOME}:${PATH}

RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip \
        && mkdir ${APPS_ROOT}/gatk \
        && unzip gatk-${GATK_VERSION}.zip \
        && mv gatk-${GATK_VERSION} ${APPS_ROOT}/gatk/${GATK_VERSION} \
        && rm gatk-${GATK_VERSION}.zip

###############################################
#R = 'r/intel/3.4.2'
# INSTALLED MOST CURRENT R

###############################################
#HTSLIB 1.9
ENV HTSLIB_VERSION 1.9
ENV HTSLIB_HOME ${APPS_ROOT}/htslib/${HTSLIB_VERSION}

ENV MANPATH $MANPATH:${HTSLIB_HOME}/share/man
ENV PATH ${PATH}:${HTSLIB_HOME}/bin
ENV LD_LIBRARY_PATH ${HTSLIB_HOME}/lib:${LD_LIBRARY_PATH}
ENV PKG_CONFIG_PATH ${HTSLIB_HOME}/lib/pkgconfig
ENV HTSLIB_HOME ${HTSLIB_HOME}
ENV HTSLIB_INC ${HTSLIB_HOME}/include
ENV HTSLIB_LIB ${HTSLIB_HOME}/lib

RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
	&& tar xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
	&& rm htslib-${HTSLIB_VERSION}.tar.bz2 \
	&& cd htslib-${HTSLIB_VERSION} \
	&& autoheader \
	&& autoconf  \
	&& ./configure --prefix=${HTSLIB_HOME} \
	&& make \
	&& make install

###############################################
#SAMTOOLS = 'samtools/intel/1.9'

ENV SAMTOOLS_VERSION 1.9
ENV SAMTOOLS_HOME ${APPS_ROOT}/samtools/${SAMTOOLS_VERSION}

ENV MANPATH ${SAMTOOLS_HOME}/share/man
ENV PATH ${SAMTOOLS_HOME}/bin:${PATH}
ENV LD_LIBRARY_PATH ${SAMTOOLS_HOME}/lib:${LD_LIBRARY_PATH}
ENV SAMTOOLS_HOME ${SAMTOOLS_HOME}
ENV SAMTOOLS_INC ${SAMTOOLS_HOME}/include
ENV SAMTOOLS_LIB ${SAMTOOLS_HOME}/lib

RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& autoheader \
	&& autoconf -Wno-syntax \
	&& ./configure --prefix=${SAMTOOLS_HOME} --with-htslib=${HTSLIB_HOME} \
	&& make \
	&& make install

###############################################
#SNPEFF = 'snpeff/4.3'

ENV SNPEFF_VERSION 4_3i
ENV SNPEFF_HOME ${APPS_ROOT}/snpeff/${SNPEFF_VERSION}

ENV SNPEFF_JAR ${SNPEFF_HOME}/snpEff.jar
ENV SNPSIFT_JAR ${SNPEFF_HOME}/SnpSift.jar

RUN wget -O snpEff_v${SNPEFF_VERSION}_core.zip  https://sourceforge.net/projects/snpeff/files/snpEff_v${SNPEFF_VERSION}_core.zip/download# \
        && mkdir ${APPS_ROOT}/snpeff \
        && unzip snpEff_v${SNPEFF_VERSION}_core.zip \
        && mv snpEff ${APPS_ROOT}/snpeff/${SNPEFF_VERSION}

###############################################
# R Packages Installation
RUN R -e "install.packages('ggplot2', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('gsalib', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('reshape', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('gplots', repos = 'http://cran.us.r-project.org')"
