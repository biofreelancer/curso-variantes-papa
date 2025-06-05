#!/bin/bash
nextflow run nf-core/sarek \
	-profile test,docker \
	-r 3.5.1 \
	--outdir prueba-sarek
