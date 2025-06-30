nextflow run nf-core/sarek \
    -r 3.5.1 \
    --input  3-samplesheets/1000genomas_samples.csv \
    --outdir results/1000genomas \
    --genome 'GATK.GRCh38' \
    --tools strelka \
    -profile docker \
    -resume
