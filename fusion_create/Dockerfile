FROM ubuntu:16.04

RUN apt-get update
RUN apt-get -y install wget tar g++ make python python-pip zlib1g-dev

WORKDIR /opt

RUN pip install pip --upgrade
RUN pip install numpy argparse gffutils biopython

COPY *.py* /opt
RUN chmod +x *.py*
ENV PATH /opt:$PATH
 
RUN wget https://github.com/deweylab/RSEM/archive/v1.2.31.tar.gz && \
    tar -zxvf v1.2.31.tar.gz && \
    cd RSEM-1.2.31/ && \
    make && \
    make install

RUN wget https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz && \
    tar -zxvf STAR_2.4.2a.tar.gz && \
    cp /opt/STAR-STAR_2.4.2a/bin/Linux_x86_64/* /usr/local/bin

