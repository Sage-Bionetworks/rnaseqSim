FROM ubuntu:16.04

RUN apt-get update

#RUN apt-get -y install software-properties-common gnupg

#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
#RUN add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -
RUN apt-get update
RUN apt-get -y install r-base git

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('argparse')"
RUN Rscript -e "install.packages('data.table')"

WORKDIR /opt

RUN git clone https://github.com/alliecreason/rnaseqSim.git  && \
    chmod +x model_isoforms/*.R

ENV PATH /opt/rnaseqSim/model_isoforms:$PATH
