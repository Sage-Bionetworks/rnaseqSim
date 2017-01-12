FROM ubuntu:16.04

RUN apt-get update
RUN apt-get -y install git wget tar g++ make python python-pip zlib1g-dev

RUN pip install pip --upgrade
RUN pip install argparse


WORKDIR /opt

COPY *.py  /opt 
RUN chmod +x /opt/*
ENV PATH /opt:$PATH
