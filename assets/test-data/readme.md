## Test-data

Below is how the test data is created.

### Reference

Reference genome file

```bash
cp ~/icloud/Projects/2016-tiger-natera/panthera-10x-170301.fasta ./
bgzip -@8 panthera-10x-170301.fasta
samtools faidx panthera-10x-170301.fasta.gz
sort -grk2 panthera-10x-170301.fasta.gz.fai | head > chroms.txt
samtools faidx panthera-10x-170301.fasta.gz 790 833 1172 1310 > pantig.fasta
bgzip -f -@8 pantig.fasta
samtools faidx pantig.fasta.gz
rm panthera*
```


### VCF

VCF of all known variants in population. Used to mask the reference genome so that primers can be designed to not overlap known variants.

```bash
tabix ~/icloud/Projects/2016-tiger-wgs/vcfs/fb-10x-180612-s10-snp-q30-gq30-ma2-hwe-maf01.vcf.gz
bcftools view -r790 -r833 -r1172 -r1310 ~/icloud/Projects/2016-tiger-wgs/vcfs/fb-10x-180612-s10-snp-q30-gq30-ma2-hwe-maf01.vcf.gz > pantig.vcf
bcftools sort pantig.vcf > pantig.sorted.vcf
bgzip -f -@8 pantig.sorted.vcf
tabix pantig.sorted.vcf.gz
rm pantig.vcf
```


### Targets

Sample single nucleotide targets to design primers for.

```bash
bcftools view --type=snps --min-af=0.3 --max-af=0.7 pantig.vcf.gz > targets.vcf
bcftools query -f '%CHROM\t%POS\t%POS\n' targets.vcf | awk '{print($1, $2, $2+1, $1":"$2)}' > targets.bed
bgzip -@8 targets.bed
tabix -fp bed targets.bed.gz
```

