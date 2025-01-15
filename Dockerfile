FROM ubuntu:18.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV MICROC=Micro-C

RUN apt-get clean
RUN apt-get update
RUN echo 'Acquire::http::Pipeline-Depth 0;\nAcquire::http::No-Cache true;\nAcquire::BrokenProxy true;\n' > /etc/apt/apt.conf.d/99fixbadproxy

# Use old-releases mirrors
RUN apt-get update --fix-missing -y && \
    apt-get upgrade -y && \
    apt-get install -y \
    inetutils-ping \
    tabix \
    build-essential \
    python3 \
    python3-pip \
    zlib1g-dev \
    libcurl4 \
    bedtools \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    liblz4-tool \
    libjpeg-dev \
    zlib1g-dev \
    wget \
    git && \
    rm -rf /var/lib/apt/lists/*

# Set python and pip alternatives
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1 && \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1


# Install Python packages
RUN pip3 install \
    pysam \
    tabulate \
    numpy \
    scipy \
    py2bit \
    matplotlib \
    pyBigWig \
    deeptools \
    pandas
RUN pip3 install pairtools

# Install BWA
RUN git clone https://github.com/lh3/bwa.git && cd bwa && make -j $(nproc) && \
    cp bwa qualfa2fq.pl xa2multi.pl /usr/local/bin && cd .. && rm -rf bwa

# Install Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
    tar xvf samtools-1.11.tar.bz2 && cd samtools-1.11 && \
    ./configure && make -j $(nproc) && make install && \
    cd htslib-1.11 && make && make install && \
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib" && cd .. && rm -rf samtools-1.11*

# Install Preseq
RUN wget https://github.com/smithlabcode/preseq/releases/download/v3.1.2/preseq-3.1.2.tar.gz && \
    tar xvf preseq-3.1.2.tar.gz && cd preseq-3.1.2 && \
    ./configure --enable-hts CPPFLAGS='-I /usr/local/include' LDFLAGS='-L/usr/local/lib' && \
    make -j $(nproc) && make install && cd .. && rm -rf preseq-3.1.2*

# Set working directory
WORKDIR /.

# Clone Micro-C repository
RUN git clone https://github.com/dovetail-genomics/Micro-C 

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# To run directly the analysis until the qc
# the file I want to copy in the image, has to be in the directory I am building the container from
COPY Analysis.sh /Analysis.sh
# Ensure the script has execute permissions
# RUN chmod +x /Analysis.sh
ENTRYPOINT ["sh","/Analysis.sh"]

# Set entrypoint
#ENTRYPOINT ["/bin/bash"]ls
