FROM r-base:4.4.2

ENV R_PROGRAM=/usr/bin/R

RUN apt-get update && \
	apt-get install -y \
		bcftools \
		bowtie \
		cmake \
		libcurl4-openssl-dev \
		libssl-dev \
		mafft \
		python3 \
		python3-pip \
		samtools \
		xauth

RUN useradd hippo
ENV HOME_DIR=/home/hippo
WORKDIR $HOME_DIR

#
# Set up external dependencies
#

ENV SAMTOOLS_PATH=/usr/bin

ENV BOWTIE_VERSION=2.5.4

RUN mkdir $HOME_DIR/bin && \
	wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/$BOWTIE_VERSION/bowtie2-$BOWTIE_VERSION-linux-x86_64.zip/download -O $HOME_DIR/bin/bowtie2.zip && \
	unzip $HOME_DIR/bin/bowtie2.zip -d $HOME_DIR/bin/ && \
	rm $HOME_DIR/bin/bowtie2.zip && \
	mv $HOME_DIR/bin/bowtie2-$BOWTIE_VERSION-linux-x86_64/* $HOME_DIR/bin/ && \
	rmdir $HOME_DIR/bin/bowtie2-$BOWTIE_VERSION-linux-x86_64

ENV BOWTIE_PATH="$HOME_DIR/bin"
ENV BOWTIE_PROGRAM="${BOWTIE_PATH}/bowtie"
ENV BOWTIE2_PATH="$HOME_DIR/bin"
ENV BOWTIE2_PROGRAM="${BOWTIE2_PATH}/bowtie2"

ENV INDEXES_PATH="/indexes"
ENV BLASTINDEX_PATH="$INDEXES_PATH/Blast_Indexes"
ENV BOWTIE_INDEX_PATH="$INDEXES_PATH/BowtieIndexes"
ENV BOWTIE2_INDEX_PATH="$INDEXES_PATH/Bowtie2Indexes"

# RUN wget "https://github.com/dzerbino/velvet/archive/refs/heads/master.zip" -O /tmp/velvet.zip && \
# 	unzip /tmp/velvet.zip -d /tmp && \
# 	rm /tmp/velvet.zip && \
# 	cd /tmp/velvet-master && \
# 	make && \
# 	cd $HOME_DIR && \
# 	rm -rf /tmp/velvet-master
# couldn't get velvet to build this way... must be missing something. move on for now.

RUN pip3 install --break-system-packages cutadapt

ENV PYTHON2_VERSION=2.7.18

RUN wget "https://www.python.org/ftp/python/${PYTHON2_VERSION}/Python-${PYTHON2_VERSION}.tgz" -O /tmp/python2.tar.gz && \
	tar -x --directory=/tmp/ --file=/tmp/python2.tar.gz --gzip --verbose && \
	cd /tmp/Python-${PYTHON2_VERSION} && \
	./configure --enable-optimizations && \
	make && \
	make install && \
	cd $HOME_DIR && \
	rm -rf /tmp/Python-${PYTHON2_VERSION}

#
# Set up DuffyTools
#

RUN wget "https://github.com/ma-lab-cgidr/DuffyTools/archive/refs/heads/master.zip" -O $HOME_DIR/DuffyTools.zip && \
	unzip $HOME_DIR/DuffyTools.zip -d $HOME_DIR/ && \
	rm $HOME_DIR/DuffyTools.zip

RUN cd $HOME_DIR/DuffyTools-* && \
	sed -i 's;${R_PROGRAM};/usr/bin/R;' Makefile && \
	make && \
	cd $HOME_DIR

#
# Set up DuffyNGS
#

RUN R -e 'install.packages("BiocManager"); BiocManager::install("Biostrings", update=TRUE, ask=FALSE)'

RUN wget "https://github.com/ma-lab-cgidr/DuffyNGS/archive/refs/heads/master.zip" -O $HOME_DIR/DuffyNGS.zip && \
	unzip $HOME_DIR/DuffyNGS.zip -d . && \
	rm $HOME_DIR/DuffyNGS.zip

RUN cd $HOME_DIR/DuffyNGS-* && \
	sed -i 's;${R_PROGRAM};/usr/bin/R;' Makefile && \
	make && \
	cd $HOME_DIR

RUN R -e 'library(DuffyTools); checkDuffyPackages();' > $HOME_DIR/duffy_packages.txt
RUN R -e 'library(DuffyNGS); checkDuffyPrograms();' > $HOME_DIR/duffy_programs.txt

#
# Finalize
#

RUN apt-get install -y xvfb xfonts-base

ADD Options.txt Options.txt
ADD pipelineSetup.R pipelineSetup.R
ADD quickqc.R quickqc.R
ADD rnaseqbatch.R rnaseqbatch.R
ADD rnaseqfin.R rnaseqfin.R
ADD wrapper.py wrapper.py

ADD xvfb xvfb
RUN chmod +x xvfb

RUN mkdir /work
RUN chown hippo:hippo /work

RUN chown hippo:hippo $HOME_DIR/*
USER hippo

# VOLUME /data
# VOLUME /indexes
ENTRYPOINT ["python3", "/home/hippo/wrapper.py"]
CMD ["--help"]

# if we want to split this into 2 build stages, will need to copy over $HOME_DIR/bin, /usr/local/bin/, and /usr/local/lib/ (at least)