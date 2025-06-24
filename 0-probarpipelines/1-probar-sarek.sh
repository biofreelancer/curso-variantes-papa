#!/bin/bash
nextflow run nf-core/sarek \
        -profile test,docker \
        -r 3.5.1 \
        --outdir prueba-sarek \
        -c a-configfiles/low-res-machine.config \
        -resume
