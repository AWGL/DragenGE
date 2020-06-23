#!/bin/bash
set -euo pipefail

# Usage: bash /data/pipelines/DragenGE/DragenGE-1.0.0/DragenGE.sh

version=1.0.0

. *.variables
. /data/pipelines/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/"$panel"/*.variables

cp /data/pipelines/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/"$panel"/*.variables ..


# make csv with fastqs in

echo "RGID,RGSM,RGLB,Lane,Read1File,Read2File" > fastqs.csv


for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do
   
   laneId=$(echo "$fastqPair" | cut -d_ -f3)
   read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
   read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

   echo "$seqId,$sampleId,$seqId,$laneId,$PWD/$read1Fastq,$PWD/$read2Fastq" >> fastqs.csv

done




/opt/edico/bin/dragen \
-r /staging/human/reference/GRCh37/ \
--output-directory . \
--output-file-prefix "$seqId"_"$sampleId" \
--output-format BAM \
--enable-map-align-output true \
--fastq-list fastqs.csv \
--fastq-list-sample-id $sampleId \
--enable-duplicate-marking true \
--enable-variant-caller true \
--vc-sample-name "$sampleId" \
--vc-target-bed /data/pipelines/$pipelineName/"$pipelineName"-"$pipelineVersion"/"$panel"/*.bed \
--vc-emit-ref-confidence BP_RESOLUTION \
--vc-target-bed-padding 100 \
--strict-mode true 

if [ -e "$seqId"_"$sampleId".hard-filtered.gvcf.gz ]; then
    echo $sampleId/"$seqId"_"$sampleId".hard-filtered.gvcf.gz >> ../gVCFList.txt
fi

# if all samples have been processed for the panel perform joint genotyping
# expected number
expGVCF=$(ls -d ../*/ | wc -l)

# observed number
obsGVCF=$(wc -l < ../gVCFList.txt)

if [ $expGVCF == $obsGVCF ]; then
    echo "$sampleId is the last sample"
    echo "performing joint genotyping"
    cd ..

    /opt/edico/bin/dragen \
        -r  /staging/human/reference/GRCh37/ \
        --output-directory . \
        --output-file-prefix "$seqId" \
        --enable-joint-genotyping true \
        --variant-list gVCFList.txt \
        --strict-mode true


    # delete gvcfs as we don't need anymore
    ls */*.gvcf.gz* | xargs rm

else
    echo "$sampleId is not the last sample"

fi
