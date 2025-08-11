#!/bin/bash

set -euo pipefail

# set max processes and open files as these differ between wren and head node
ulimit -S -u 16384
ulimit -S -n 65535


# Usage: cd /staging/data/results/$seqId/$panel/$sampleId && bash DragenWGS.sh 

version=2.7.0

##############################################
# SETUP                                      #
##############################################

pipeline_dir="/mnt/Data-MSA/diagnostics/pipelines/"
output_dir="/mnt/Data-MSA/results/"
bcftools_path="/mnt/Data-MSA/diagnostics/apps/miniconda3/envs/bcftools/bin/"


# load variables for sample and pipeline
. *.variables
. "$pipeline_dir"/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/config/"$panel"/*.variables


# copy relevant variables files to the results directory
cp "$pipeline_dir"/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/config/"$panel"/*.variables ..
cp -r "$pipeline_dir"/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/config .
cp -r "$pipeline_dir"/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/commands .

# make csv with fastqs in

echo "RGID,RGSM,RGLB,Lane,Read1File,Read2File" > fastqs.csv

for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do
   
   laneId=$(echo "$fastqPair" | cut -d_ -f3)
   read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
   read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

   echo "$seqId"_"$laneId","$sampleId","$seqId","$laneId","$read1Fastq","$read2Fastq" >> fastqs.csv

done


/opt/edico/bin/dragen \
-r $dragen_ref \
--output-directory . \
--output-file-prefix "$seqId"_"$sampleId" \
--output-format BAM \
--enable-map-align-output true \
--fastq-list fastqs.csv \
--fastq-list-sample-id $sampleId \
--enable-duplicate-marking true \
--enable-variant-caller true \
--vc-enable-joint-detection true \
--qc-cross-cont-vcf config/"$panel"/sample_cross_contamination_resource_"$genome_build".vcf \
--vc-sample-name "$sampleId" \
--vc-target-bed config/"$panel"/"$panel"_ROI_"$genome_build".bed \
--vc-emit-ref-confidence GVCF \
--vc-target-bed-padding 100 \
--strict-mode true \
--qc-coverage-region-1 config/"$panel"/"$panel"_ROI_"$genome_build".bed \
--qc-coverage-reports-1 cov_report \
--qc-coverage-filters-1 'mapq<20,bq<10' \
--enable-map-align true \
--alt-aware true

touch "$seqId"_"$sampleId".mapping_metrics.csv

if [ -e "$seqId"_"$sampleId".hard-filtered.gvcf.gz ]; then
    echo $output_dir/$seqId/$panel/raw_vcf/$sampleId/"$seqId"_"$sampleId".hard-filtered.gvcf.gz >> ../gVCFList.txt
fi

if [ -e "$seqId"_"$sampleId".bam ]; then
    echo "--bam-input "$output_dir/$seqId/$panel/alignments/$sampleId/$sampleId"/"$seqId"_"$sampleId".bam \\" >> ../BAMList.txt
fi

#Move over the sample results into the exact folder structure (would copying and then deleting be safer?)
#Alignments
if [ -d "$output_dir/$seqId/$panel/alignments/$sampleId/" ]; then
        echo "$output_dir/$seqId/$panel/alignments/$sampleId/ already exists"
else
	mkdir $output_dir/$seqId/$panel/alignments/$sampleId/
fi
mv ${seqId}_${sampleId}.bam* $output_dir/$seqId/$panel/alignments/$sampleId/
#VCFs
if [ -d "$output_dir/$seqId/$panel/raw_vcf/$sampleId/" ]; then
	echo "$output_dir/$seqId/$panel/raw_vcf/$sampleId/ already exists"
else
	mkdir $output_dir/$seqId/$panel/raw_vcf/$sampleId/
fi
mv ${seqId}_${sampleId}.hard-filtered.gvcf.gz* $output_dir/$seqId/$panel/raw_vcf/$sampleId/
#Variables
mv ${sampleId}.variables $output_dir/$seqId/$panel/variables
#Metrics
if [ -d "$output_dir/$seqId/$panel/metrics/$sampleId/" ]; then
	echo "$output_dir/$seqId/$panel/metrics/$sampleId/ already exists"
else
	mkdir $output_dir/$seqId/$panel/metrics/$sampleId/
fi
mv ${seqId}_${sampleId}.mapping_metrics.csv $output_dir/$seqId/$panel/metrics/$sampleId/
mv ${seqId}_${sampleId}.qc-coverage-region-1_coverage_metrics.csv $output_dir/$seqId/$panel/metrics/$sampleId/
mv ${seqId}_${sampleId}.vc_metrics.csv $output_dir/$seqId/$panel/metrics/$sampleId/
mv ${seqId}_${sampleId}.wgs_coverage_metrics.csv $output_dir/$seqId/$panel/metrics/$sampleId/
mv ${seqId}_${sampleId}.target_bed_coverage_metrics.csv $output_dir/$seqId/$panel/metrics/$sampleId/
#Tar up everything else
mkdir ${sampleId}_analysis
mv ${seqId}_* ${sampleId}_analysis
mv *_usage.txt ${sampleId}_analysis
mv fastqs.csv ${sampleId}_analysis
mv streaming_log_dragen.csv ${sampleId}_analysis
tar -czvf ${sampleId}_analysis.tar.gz ${sampleId}_analysis/
mv ${sampleId}_analysis.tar.gz $output_dir/$seqId/$panel/archive

# if all samples have been processed for the panel perform joint genotyping
# expected number
expGVCF=$(ls -d ../*/ | wc -l)

# observed number
obsGVCF=$(wc -l < ../gVCFList.txt)

if [ $expGVCF == $obsGVCF ]; then
    echo "$sampleId is the last sample"
    echo "performing joint genotyping"
    
    mv commands/joint_call_svs.sh ..
    mv commands/create_ped.py ..
    mv commands/by_family.py ..

    cd ..

    /opt/edico/bin/dragen \
        -r  $dragen_ref \
        --output-directory . \
        --output-file-prefix "$seqId" \
        --enable-joint-genotyping true \
        --vc-enable-joint-detection true \
        --variant-list gVCFList.txt \
        --strict-mode true


    
    if [ $callSV == true ]; then

        echo Joint Calling SVs

	python create_ped.py --variables '/mnt/Data-MSA/results/'"$seqId"'/'"$panel"'/variables/*.variables' > "$seqId".ped

        python by_family.py "$seqId".ped "$seqId" "$panel"

        mkdir sv_calling

        for family in *_for_sv.family; do        

            cp joint_call_svs.sh joint_call_svs.sh_"$family".sh
            cat $family >> joint_call_svs.sh_"$family".sh
            bash joint_call_svs.sh_"$family".sh $family $panel $dragen_ref $fasta
            rm joint_call_svs.sh_"$family".sh
        done

	#Combining family SV files - runs from conda environment installed on mount. If statement only runs bcftools if more than one family. bcftools merge crashes with a single vcf. 

        if [ `ls -1 sv_calling/*vcf.gz | wc -l` -eq 1 ]; then
            cp sv_calling/*.vcf.gz "$seqId".sv.vcf.gz
        else
            ${bcftools_path}/bcftools merge -m none -F x sv_calling/*.vcf.gz > "$seqId".sv.vcf
            ${bcftools_path}/bgzip "$seqId".sv.vcf
        fi
        
        ${bcftools_path}/tabix "$seqId".sv.vcf.gz

        md5sum "$seqId".sv.vcf.gz | cut -d" " -f 1 > "$seqId".sv.vcf.gz.md5sum

        rm -r sv_calling
        rm *.family
        rm create_ped.py
        rm by_family.py	
    fi   

    #Copy everything else
    #ped
    if [ -d "$output_dir/$seqId/$panel/ped/" ]; then
        echo "$output_dir/$seqId/$panel/ped/ already exists"
    else
        mkdir $output_dir/$seqId/$panel/ped/
    fi
    mv ${seqId}.ped $output_dir/$seqId/$panel/ped/
    #SV
    mv ${seqId}.sv.vcf.gz* $output_dir/$seqId/$panel/raw_sv_vcf/
    #VCFs
    mv ${seqId}.vcf.gz* $output_dir/$seqId/$panel/raw_vcf/
    mv ${seqId}.hard-filtered.vcf.gz* $output_dir/$seqId/$panel/raw_vcf/
    #Variables
    mv ${panel}.variables $output_dir/$seqId/$panel/
    #Metrics
    mv ${seqId}.time_metrics.csv $output_dir/$seqId/$panel/metrics
    mv ${seqId}.vc_hethom_ratio_metrics.csv $output_dir/$seqId/$panel/metrics
    mv ${seqId}.vc_metrics.csv $output_dir/$seqId/$panel/metrics


    # mark results as complete - do this first so post processing can start asap
    touch "$output_dir"/"$seqId"/"$panel"/dragen_complete.txt
    touch "$output_dir"/"$seqId"/"$panel"/post_processing_required.txt

    # clean up staging results
    rm -r /staging/data/results/"$seqId"/"$panel"
    # clean up fastq will depends where this has been written
    if [ -d /staging/data/fastq/"$seqId"/Data/"$panel" ]; then
	fastq_path=/staging/data/fastq/"$seqId"/
    else
	fastq_path=/mnt/Data-MSA/results/"$seqId"/fastq/
    fi
    rm -r "$fastq_path"/Data/"$panel"

    # clean up staging fastq if we have processed all panels
    if [ "$(ls -A "$fastq_path"/Data)" ]; then
        echo "Not all panels processed - keeping staging fastq"
    else
        echo "All panels processed - removing staging fastq directory"
        rm -r $fastq_path
    fi

    # clean up staging results if we have processed all panels
    if [ "$(ls -A /staging/data/results/"$seqId"/)" ]; then
        echo "Not all panels processed - keeping staging results"
    else
        echo "All panels processed - removing staging results directory"
        rm -r /staging/data/results/"$seqId"/

    fi

    # Remove lock file from dragen
    rm /mnt/Data-MSA/raw/dragen_markers/${seqId}_*_locked


else
    echo "$sampleId is not the last sample"

fi
