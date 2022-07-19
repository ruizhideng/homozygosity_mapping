#!/bin/bash

usage() { echo "Usage: $0 [-v <string>] [-o <string>] [-f <string>]" 1>&2; exit 1; }
help() { echo "
Version: 1.0.0

This sortware is used to search for run of homozygosity regions


usage: bash test.sh [-v <string>] [-o <string>] [-f <string>]

Options:
-v     vcf files, multiple files vcshould be seperated by comma;
       e.g A.vcf,B.vcf,C.vcf
-o     output directory
-f     You can provide a fam file for plink, including FID, IID, fatherID, motherID, sex; 
       If you don't have it, the software will define it. 
       See the fomart on the website https://www.cog-genomics.org/plink/1.9/formats#fam
" 1>&2; exit 1; }


while getopts ":d:v:o:hf::" opt; do
    case "$opt" in 
        v)vcf="$OPTARG"
        vcfs=(${vcf//,/ })
        numbervcf=${#vcfs[@]};;
        f) fam="$OPTARG";;
        o) out="$OPTARG";;
        h) help;; 
        *) usage;;
     esac
done


# prepare vcf in correct format
mkdir -p $out

if [ ! $fam ];then
    echo "--> No use of -f option, the phenotype set as default: ambiguous"
else
    echo "--> your phenotype file is ${fam}"
    cp $fam $out
fi

# Variants were fltered for a minimum depth of coverage of at least 10 reads and a genotype quality of at least 50.
# split multi-allelic records
for (( k=0; k<$numbervcf; k++ )); do
    vcf=${vcfs[$k]}
    bgzip -c $vcf > $out/$vcf.gz
    tabix -p vcf $out/$vcf.gz
    bcftools view -i 'MIN(FMT/DP)>10 & MIN(FMT/GQ)>50' $out/$vcf.gz -o  $out/$vcf.DP10GQ50.vcf.gz -O z
    bcftools norm -m-any $out/$vcf.DP10GQ50.vcf.gz -o $out/$vcf.DP10GQ50.split.vcf.gz -O z
    tabix -p vcf $out/$vcf.DP10GQ50.split.vcf.gz
done

# merge multiple vcfs into one 
if (( $numbervcf > 1 ));then
    bcftools merge $out/*DP10GQ50.split.vcf.gz -o $out/merged.vcf.gz -O z
else   
    mv  $out/*DP10GQ50.split.vcf.gz $out/merged.vcf.gz
fi

# PLINK detects ROH
cd $out
plink --vcf merged.vcf.gz --make-bed --out ex --allow-extra-chr # combine multiple vcfs

if [ -z $fam ];then
    plink --bfile ex --homozyg group --homozyg-gap 100000 --homozyg-kb 100 --homozyg-window-het 0 --homozyg-window-snp 5 --out out  --allow-extra-chr
    # plink --bfile ex --homozyg group --homozyg-gap 10000 --homozyg-kb 1000 --homozyg-window-het 1 --homozyg-window-missing 10 --homozyg-window-snp 20 --homozyg-density 10000 --homozyg-snp 10 --homozyg-window-threshold 0.05 --out out  --allow-extra-chr
else
    plink --bfile ex --fam $fam --homozyg group --homozyg-gap 100000 --homozyg-kb 100 --homozyg-window-het 0 --homozyg-window-snp 5 --out out  --allow-extra-chr 
    # plink --bfile ex --fam $fam --homozyg group --homozyg-gap 10000 --homozyg-kb 1000 --homozyg-window-het 1 --homozyg-window-missing 10 --homozyg-window-snp 20 --homozyg-density 10000 --homozyg-snp 10 --homozyg-window-threshold 0.05 --out out  --allow-extra-chr
fi

rm *gz* ex*

echo "***Finished***"