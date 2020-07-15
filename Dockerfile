FROM ubuntu:bionic

RUN echo 'cache'
RUN apt-get update && apt-get install -y apt-utils
RUN apt-get install -y bedtools python3 python3-pip
RUN pip3 install click pyfaidx primer3-py
RUN apt-get clean
COPY bin/mpcrutils /usr/local/
RUN ln -s /usr/local/mpcrutils.py3 /usr/local/bin/mpcrutils.py3
RUN chmod +x /usr/local/bin/mpcrutils.py3
RUN apt-get install -y tabix
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
RUN pip3 install biopython
RUN apt-get install -y wget && \
    wget https://github.com/quwubin/MFEprimer-3.0/releases/download/v3.2.0/mfeprimer-3.2.0-linux-amd64.gz && \
    gunzip mfeprimer-3.2.0-linux-amd64.gz && \
    mv mfeprimer-3.2.0-linux-amd64 /usr/local/bin/mfeprimer && \
    chmod +x /usr/local/bin/mfeprimer
