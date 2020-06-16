FROM ubuntu:bionic

RUN apt-get update && apt-get install -y apt-utils
RUN apt-get install -y bedtools python3 python3-pip
RUN pip3 install click pyfaidx primer3-py
RUN apt-get clean
COPY bin/mpcrutils /usr/local/
RUN ln -s /usr/local/mpcrutils.py3 /usr/local/bin/mpcrutils.py3
RUN chmod +x /usr/local/bin/mpcrutils.py3
RUN apt-get install -y tabix