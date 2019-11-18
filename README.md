## admixture
To evaluate admixture in this Brazilian sample we carried out principal components analysis on an LD-pruned set of SNVs with minor allele frequencies >0.1.  Comparison with HapMap populations indicated predominant admixture between Caucasian and Negroid populations (Figure 2A), consistent with data from the ABraOM database of exome variants from 609 elderly Brazilians from Sao Paulo State.5

First convert your vcf file to plink file using for example vcftools (http://vcftools.sourceforge.net/)

```{bash}
./vcftools --vcf input_data.vcf --plink --out output_in_plink
```

Then download the hapmap data from plink:
```{bash}
http://pngu.mgh.harvard.edu/~purcell/plink/dist/hapmap_r23a.zip
```

Then extract snps lists from both datasets and filter based on the these snplist to get only overlapping snps.

```{bash}
plink --bfile fileA--write-snplist --out list1 --noweb
plink --bfile fileB --extract list1 --noweb --out fileB_filtered --make-bed
```

Then merge the files and make a mds plot

```{bash}
plink --bfile fileA_filtered --bmerge fileB_filtered.bed fileB_filtered.bim fileB_filtered.fam --noweb --out merged --make-bed
plink --bfile merged --mds-plot 2 --noweb --out mds
```

This can be easily plotted using R or even Excel... and then you can see which samples are derived from which ethic background.

