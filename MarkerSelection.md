# Instructions for creating custom marker files
 
## What you need:
- intervaltxt: Target intervals file (1-based) restricted to autosomes only. 
- kgpbed: 1000 genomes phase 3 PLINK files
  - Can be downloaded from [here](https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3)
  - Convert to PLINK bed file format by following instructions [here](https://cran.r-project.org/web/packages/plinkQC/vignettes/Genomes1000.pdf)
- [1000 genomes sites VCF and index files](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/):
  - ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
  - ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi
- gnomad.genomes.r2.1.1.sites.vcf.bgz and gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi
- gnomad.exomes.r2.1.1.sites.vcf.bgz and gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi
- plink
- bcftools
- gatk4

## Commands:

### Marker selection for ADMIXTURE:

#### Select autosomal bi-allelic SNPs with MAF>0.01 and within target intervals
`plink --bfile $kgpbed 
  --maf 0.01 
  --extract 'range' $intervaltxt 
  --out 1000genomes_maf0.01_snps_within_${intervalID} 
  --allow-extra-chr 
  --snps-only 
  --biallelic-only 'strict'`

#### Run LD pruning
`plink --bfile 1000genomes_maf0.01_snps_within_${intervalID} 
  --indep-pairwise 1000 100 0.2 
  --out 1000genomes_maf0.01_snps_within_${intervalID}`

#### Create PLINK bed file of markers in linkage equilibrium
`plink --bfile 1000genomes_maf0.01_snps_within_${intervalID} 
  --extract 1000genomes_maf0.01_snps_within_${intervalID}.prune.in 
  --out 1000genomes_maf0.01_snps_within_${intervalID}_ldpruned 
  --make-bed`

#### Create VCF file
`cut -f 2 1000genomes_maf0.01_snps_within_${intervalID}_ldpruned.bim > 1000genomes_maf0.01_snps_within_${intervalID}_ldpruned.snpIDs.txt`

`bcftools view -i ID=@1000genomes_maf0.01_snps_within_${intervalID}_ldpruned.snpIDs.txt ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz > ADMIXTURE_markers.vcf`

### Marker selection for Ashkenazi Jewish ancestry inference:

#### Select bi-allelic SNPs from gnomAD genomes

`bcftools view -f PASS -m2 -M2 -v snps '-iAF_asj>0.01 && AF_nfe < 0.001 && AF<0.001 && AF_eas < 0.001 && AF_afr < 0.001' -R ${intervaltxt} gnomad.genomes.r2.1.1.sites.vcf.bgz > ASJ_markers_temp1.vcf`

`grep -v "#" ASJ_markers_temp1.vcf | cut -f 3 > ASJ_markers.rsids.txt`

#### Extract MAFs from gnomAD exomes and filter

`bcftools view -i ID=@ASJ_markers.rsids.txt gnomad.exomes.r2.1.1.sites.vcf.bgz | bcftools query -f '%ID\t%INFO/AN\t%INFO/AF\t%INFO/AF_nfe\t%INFO/AF_asj\t%INFO/AF_sas\t%INFO/AF_eas\t%INFO/AF_afr\n' > ASJ_markers_genomadexome_maf.tsv`

- Filter variants where AF_asj > 0 and AF_asj/(max(AF, AF_asj, AF_nfe, AF_eas, AF_afr, AF_sas)) <= 2 and create a file of remaining IDs (ASJ_markers.rsids.filtered.txt)

#### LD prune (optional)
- Genotype these markers on samples that are known to be of Ashkenazi Jewish ancestry
- Use plink --indep-pairwise 1000 100 0.2 to LD prune and restrict to prune.in IDs (ASJ_markers.rsids.filtered.pruned.txt)

#### Create VCF file
`bcftools view -i ID=@ASJ_markers.rsids.filtered(.pruned).txt gnomad.genomes.r2.1.1.sites.vcf.bgz > ASJ_markers.vcf`

### Create final files and update config file

- Combine ADMIXTURE_markers.vcf and ASJ_markers.vcf. Use this as *markers_vcf* in config.yaml
- Use gatk IndexFeatureFile to index the above combined file
- Create a txt file from the above combined VCF file, which has the following columns: CHROM, POS, ID, REF, ALT and AF (Allele frequency). Use this file as *markers_txt* in config.yaml
- Use final ASJ marker IDs file (ASJ_markers.rsids.filtered(.pruned).txt) as *asj_markers* in config.yaml
- Use 1000genomes_maf0.01_snps_within_${intervalID}_ldpruned as *kgbed* in config.yaml
