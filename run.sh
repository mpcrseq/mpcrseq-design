#!/bin/bash

nextflow run main.nf -resume -c nextflow.config -with-timeline -with-trace
