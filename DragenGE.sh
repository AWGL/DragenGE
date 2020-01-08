#!/bin/bash
set -euo pipefail

# Usage: bash /data/pipelines/DragenGE/DragenGE-1.0.0/DragenGE.sh /staging/data/fastq/191010_D00501_0366_BH5JWHBCX3/Data/IlluminaTruSightOne/18M01315

version=1.0.0
sampleDir=$1

. *.variables
. /data/pipelines/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/"$panel"/*.variables


# make output dir for results
mkdir -p /staging/data/results/$seqId/$panel/$sampleId

cp *.variables /staging/data/results/$seqId/$panel/$sampleId
cp /data/pipelines/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/"$panel"/*.variables /staging/data/results/$seqId/$panel


# make csv with fastqs in

echo "RGID,RGSM,RGLB,Lane,Read1File,Read2File" > fastqs.csv


for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do
   
   laneId=$(echo "$fastqPair" | cut -d_ -f3)
   read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
   read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

   echo "$seqId,$sampleId,$seqId,$laneId,$PWD/$read1Fastq,$PWD/$read2Fastq" >> fastqs.csv

done




dragen \
-r /staging/human/reference/GRCh37/ \
--output-directory /staging/data/results/$seqId/$panel/$sampleId \
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

if [ -e /staging/data/results/$seqId/$panel/$sampleId/"$seqId"_"$sampleId".hard-filtered.gvcf.gz ]; then
    echo /staging/data/results/$seqId/$panel/$sampleId/"$seqId"_"$sampleId".hard-filtered.gvcf.gz >> /staging/data/results/$seqId/$panel/gVCFList.txt
fi

# if all samples have been processed for the panel perform joint genotyping
# expected number
expGVCF=$(ls -d /staging/data/fastq/$seqId/Data/$panel/*/ | wc -l)

# observed number
obsGVCF=$(wc -l < /staging/data/results/$seqId/$panel/gVCFList.txt)

if [ $expGVCF == $obsGVCF ]; then
    echo "$sampleId is the last sample"
    echo "performing joint genotyping"

    dragen \
        -r  /staging/human/reference/GRCh37/ \
        --output-directory /staging/data/results/$seqId/$panel \
        --output-file-prefix "$seqId" \
        --enable-joint-genotyping true \
        --variant-list /staging/data/results/$seqId/$panel/gVCFList.txt \
        --strict-mode true


    # delete gvcfs as we don't need anymore
    ls /staging/data/results/$seqId/$panel/*/*.gvcf.gz* | xargs rm

else
    echo "$sampleId is not the last sample"

fi
