#Start with base Ubuntu image and install python
FROM ubuntu:24.04

#Add maintainer information
LABEL maintainer="Andrew Ramos" 
LABEL maintainer="ar39@sanger.ac.uk"

#Avoid promtps from apt 
ARG DEBIAN_FRONTEND=noninteractive

#Add image creation date 
ARG WHEN

# install basic dependencies, cleanup apt garbage.
RUN apt-get update \
    && apt-get upgrade -y \    
    && apt-get install -y \
        build-essential \
	curl \
	git \
	wget \
	default-jdk \
	python3 \
	python3-pip \
	python3-venv \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set up a python environment at /opt/venv
RUN python3 -m venv /opt/venv

# Prepend environment bin to PATH (use environment python and pip)
ENV PATH="/opt/venv/bin:$PATH"

# Install python packages with pip
RUN pip install \
	SigProfilerPlotting==1.4.1

#Set working directory in the container
WORKDIR /usr/src/app

#Test package installation and access 
RUN python3 -c 'import sigProfilerPlotting; print("SigProfilerPlotting package is ready!")'
