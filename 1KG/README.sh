#!/bin/bash

# from here: https://www.biostars.org/p/335605/

main() {
# requires >=v1.9
plinkbindir='/data/workspace/richard/software/plink-1.9-x86_64_20191028'

## 1, Download the files as VCF.gz (and tab-indices)
get1kgvcf;

## 2, Download 1000 Genomes PED file
get1kgped;

## 3, Download the GRCh37 / hg19 reference genome
getgrch37;

## 4, Convert the 1000 Genomes files to BCF
1kg2bcf;

## 4a, Convert the zika data to BCF
zika2bcf;

## 5, Convert the BCF files to PLINK format
bcf2plink;

## 5a, Convert the BCF files to PLINK format
zikabcf2plink;

## 6, Exclude variants not on the coding strand
# NB - This step is only for microarray studies where the probes may only target one strand or the other (sense or non-sense)

## 7, Prune variants from each chromosome
prune;

## 7a, Prune variants from each chromosome
zikaprune;

## 8, Get a list of all PLINK files
linkPlink;

## 9, Merge all projects into a single PLINK file
mergeFiles;

## 10, Perform PCA
pcaData;
}

get1kgvcf(){
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" ;
suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ;

for chr in {1..22};
do
	wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ;
done
}

get1kgped(){
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped ;
}

getgrch37(){
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ;
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai ;
gunzip human_g1k_v37.fasta.gz ;
}

1kg2bcf(){
for chr in {1..22}; do
    bcftools norm -m-any --check-ref w -f human_g1k_v37.fasta \
    ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | \

    bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \

    bcftools norm -Ob --rm-dup both \
    > ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;

    bcftools index ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
done
}

zika2bcf(){
( grep -e "^#" ../../tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.vcf | sed 's/=chr/=/'; grep -v "^#" ../../tmp/combined_HaplotypeCaller.d.n.vep.PASSonly.vcf | sed 's/^chr//';  ) > combined_HaplotypeCaller.d.n.vep.PASSonly.nochr.vcf

bcftools norm -m-any --check-ref w -f human_g1k_v37.fasta \
	combined_HaplotypeCaller.d.n.vep.PASSonly.nochr.vcf | \
	bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
	bcftools norm -Ob --rm-dup both \
	> combined_HaplotypeCaller.d.n.vep.PASSonly.nochr.bcf ;

bcftools index combined_HaplotypeCaller.d.n.vep.PASSonly.nochr.bcf ;
}

bcf2plink(){
for chr in {1..22}; do
    "${plinkbindir}"/plink --noweb --bcf ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf \
    --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed \
    --out ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;
done
}

zikabcf2plink(){
"${plinkbindir}"/plink --noweb --bcf combined_HaplotypeCaller.d.n.vep.PASSonly.nochr.bcf \
	--keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed \
	--out combined_HaplotypeCaller.d.n.vep.PASSonly.nochr ;
}

prune(){
# --maf 0.10, only retain SNPs with MAF greater than 10%
# --indep [window size] [step size/variant count)] [Variance inflation factor (VIF) threshold]
# e.g. indep 50 5 1.5, Generates a list of markers in approx. linkage equilibrium - takes 50 SNPs at a time and then shifts by 5 for the window. VIF (1/(1-r^2)) is the cut-off for linkage disequilibrium

mkdir Pruned ;

for chr in {1..22}; do
    "${plinkbindir}"/plink --noweb --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
    --maf 0.10 --indep 50 5 1.5 \
    --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;

    "${plinkbindir}"/plink --noweb --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
    --extract Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.prune.in --make-bed \
    --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;
done
}

zikaprune(){
# --maf 0.10, only retain SNPs with MAF greater than 10%
# --indep [window size] [step size/variant count)] [Variance inflation factor (VIF) threshold]
# e.g. indep 50 5 1.5, Generates a list of markers in approx. linkage equilibrium - takes 50 SNPs at a time and then shifts by 5 for the window. VIF (1/(1-r^2)) is the cut-off for linkage disequilibrium

"${plinkbindir}"/plink --noweb --bfile combined_HaplotypeCaller.d.n.vep.PASSonly.nochr \
--maf 0.10 --indep 50 5 1.5 \
--out Pruned/combined_HaplotypeCaller.d.n.vep.PASSonly.nochr ;

"${plinkbindir}"/plink --noweb --bfile combined_HaplotypeCaller.d.n.vep.PASSonly.nochr \
--extract Pruned/combined_HaplotypeCaller.d.n.vep.PASSonly.nochr.prune.in --make-bed \
--out Pruned/combined_HaplotypeCaller.d.n.vep.PASSonly.nochr ;
}

linkPlink(){
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list.zika ;
sed -i 's/.bim//g' ForMerge.list.zika ;
}

mergeFiles(){
"${plinkbindir}"/plink --merge-list ForMerge.list.zika --out Merge.zika ;
}

pcaData(){
"${plinkbindir}"/plink --bfile Merge.zika --pca
}

main;
