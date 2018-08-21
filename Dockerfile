# Dockerfile for MEGAlib
# 
# Build the docker with:
# docker build -t megalib-experimental - < Dockerfile

FROM ubuntu:18.04

MAINTAINER Andreas Zoglauer <zoglauer@berkeley.edu>

# Install all the MEGAlib, ROOT, and Geant4 prerequisites
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -yq gosu vim nano less gzip git gawk dpkg-dev make g++ gcc gfortran gdb valgrind binutils libx11-dev libxpm-dev libxft-dev libxext-dev libssl-dev libpcre3-dev libglu1-mesa-dev libglew-dev libftgl-dev libmysqlclient-dev libfftw3-dev graphviz-dev libavahi-compat-libdnssd-dev libldap2-dev python-dev libxml2-dev libkrb5-dev libgsl-dev cmake libxmu-dev curl doxygen libblas-dev liblapack-dev expect dos2unix libncurses5-dev python-numpy python-scipy python-sklearn python-pip

# Add Mr. MEGAlib user
RUN groupadd -g 2468000 mrmegalib && useradd -u 2468000 -g 2468000 -ms /bin/bash mrmegalib

# Create ROOT directory
RUN mkdir /opt/MEGAlib && chown -R mrmegalib:mrmegalib /opt/MEGAlib

# Switch to Mr. MEGAlib
USER mrmegalib

# Setup MEGAlib
RUN cd /home/mrmegalib && git clone https://github.com/zoglauer/megalib.git MEGAlib && cd /home/mrmegalib/MEGAlib && /bin/bash setup.sh --ex=/opt/MEGAlib --release=dev --branch=experimental --clean=yes
RUN echo . /home/mrmegalib/MEGAlib/bin/source-megalib.sh >> /home/mrmegalib/.bashrc

# Switch back to ROOT
USER root

# Setup rootpy
RUN /bin/bash -c ". /home/mrmegalib/MEGAlib/bin/source-megalib.sh && pip install rootpy"

# Create entry-point script - changes the UID of mrmegalib to the local USER and group for full access to the exchange directory on all machines
RUN    cd /usr/local/bin \
    && echo '#!/bin/bash' >> entrypoint.sh \
    && echo 'USERID=${USERID:-9001}' >> entrypoint.sh \
    && echo 'GROUPID=${GROUPID:-9001}' >> entrypoint.sh \
    && echo 'usermod -u ${USERID} mrmegalib' >> entrypoint.sh \
    && echo 'groupmod -g ${GROUPID} mrmegalib' >> entrypoint.sh \
    && echo 'chown -R mrmegalib:mrmegalib /home/mrmegalib' >> entrypoint.sh \
    && echo 'gosu mrmegalib bash' >> entrypoint.sh \
    && chmod a+rx /usr/local/bin/entrypoint.sh

# The working directory is the home directory
WORKDIR /home/mrmegalib

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
