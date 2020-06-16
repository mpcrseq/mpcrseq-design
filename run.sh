#!/bin/bash

nextflow run main.nf -resume -profile test,docker -with-timeline -with-trace 
