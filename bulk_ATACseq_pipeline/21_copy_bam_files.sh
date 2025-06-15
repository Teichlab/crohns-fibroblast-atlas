
# Copy bam and index from filtered alignments output from ENCODE ATAC-seq pipeline
# in the 'call-filter' subdirectories and not a temporary 'glob' directory.
mkdir -p bams
find ./results -name "*nodup.no_chrM_MT.bam*" -path *call-filter* -not -path *glob* -exec cp {} ./bams \;