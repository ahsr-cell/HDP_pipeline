FROM ubuntu:24.04

LABEL maintainer="Andrew Ramos" 
LABEL maintainer="ar39@sanger.ac.uk"

ARG DEBIAN_FRONTEND=noninteractive

ARG WHEN

# install basic dependencies, cleanup apt garbage.
RUN apt-get update \
    && apt-get upgrade -y \    
    && apt-get install -y \
        build-essential \
        curl \
        git \
        wget \
        libgit2-dev \
        r-base \
        r-cran-devtools \
        r-cran-biocmanager \
        r-cran-tidyverse \
        r-cran-argparse \
    && apt-get -y install git-all \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN R -e "devtools::install_github('NickWilliamsSanger/hdp')"

RUN R -e "library(dplyr); \
    library(hdp); \ 
    library(argparse); \
    print('All packages are ready to be used!')"
