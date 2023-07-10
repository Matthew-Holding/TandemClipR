# TandemClipR

![Gene Expression](/gene_plus_RNA.png "Annotation/Alignment Image")

Software required:

* [GNU Parallel](https://www.gnu.org/software/parallel/)
* [samtools](http://www.htslib.org/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* R
* RStudio
* For R scripts, install and load the following R Packages
  * biomaRt
  * biomartr
  * rtracklayer
  * tidyverse
  * plyr
  * ggpubr
  * nlme
  * ape
  * phytools
  * geiger
  * caper

## Contents

* [Get Ensembl Genome Files](#get-ensembl-genome-files)
* [Set target gene family](#set-target-gene-family)
* [Get orthologous gene arrays](#get-orthologous-bookend-gene-coordinates)
* [Write clipped gff and bedfiles for gene arrays](#write-clipped-gff-and-bedfiles-for-gene-arrays)
* [Calculate genome stats with bedtools](#calculate-genome-stats-with-bedtools)

### Get Ensembl Genome Files
Run in bash command line

```
#download genomes and gff3 annotations from Ensembl with rsync command line

mkdir genomes
rsync -av rsync://https://ftp.ensembl.org/pub/release-107/fasta/*/dna/*fa.gz genomes/

mkdir gffs
rsync -av rsync://https://ftp.ensembl.org/pub/release-107/gff3/*/*gff3* gffs/EnsemblVertebrates/

```

### Set target gene family
Run in R

```
#make a vector of target gene family headers
genefam=c("KRTAP","HOXA","HOXD","KRTAP_2","IL1R","GABR","MYH","THEM_LCE","SERPINA","CYP2", "HBB_OR")

#edit below to specify tandem array bookends (with official gene names) that provide microsynteny across species
#make list order the same as for the genefam object above

genes <- list(c("CLDN8","TIAM1"),c("SKAP2","EVX1"),
              c("EVX2","MTX2"),c("TNS4","EIF1"),
              c("MAP4K4","SLC9A4"),c("ATP10B","NUDCD2"),
              c("GAS7","SCO1"),c("RORC","CHTOP"),
              c("PPP4R4","GSC"),c("EGLN2","AXL"),
              c("RRM1","FHIP1B"))

#Set g# below to access individual members of the genefam list for code below. Could be changed to a for loop

g=1     # will analyze KRTAP bookened by CLDN8 and TIAM1
      
```


### Get orthologous bookend gene coordinates 
Run in R

```
#link to Ensembl mart
ensembl = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", mirror = "useast")

#retrieved genomes list with below code, but then edited with metadata and reread
genomes <- listDatasets(ensembl)
genomes <- read.table("genomes.txt", header=TRUE, sep = "\t") #find as file in this repository

#get human data
human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "https://dec2021.archive.ensembl.org/")


#Now we will link chromosome locations for the human bookends with the ortholog locations in other species genomes using biomaRt getLDS() function
print(paste0("Time for ",genefam[g]))

locs.out <- data.frame(matrix(ncol = 10, nrow = 0))

for(i in 1:length(genomes$dataset)){
#tryCatch will keep loop going when species have no ortholog known and so toss an error back
tryCatch({
print(paste("Working on",genomes$dataset[i]))
tbl <- getLDS(attributes = c("hgnc_symbol"), 
       filters = "hgnc_symbol", values = genes[[g]], mart = human, 
       attributesL = c("chromosome_name","start_position","end_position","ensembl_gene_id_version"), martL = useMart("ENSEMBL_MART_ENSEMBL", dataset = genomes$dataset[i], verbose = FALSE, host = "https://dec2021.archive.ensembl.org/"), bmHeader=FALSE)
tbl1 <- tbl[tbl$hgnc_symbol==genes[[g]][1],]
tbl2 <- tbl[tbl$hgnc_symbol==genes[[g]][2],]
df <- data.frame(matrix(ncol = 10, nrow = 1))
#get just the data we want in a nice table
colnames(df) <- c("species","scientific_name","genome","version","gene1_chromosome","gene1_start","gene1_end","gene2_chromosome","gene2_start","gene2_end")
df$species[1] <- genomes$Common.name[i]
df$scientific_name[1] <- genomes$Scientific.name[i]
df$genome[1] <- genomes$dataset[i]
df$version[1] <- genomes$version[i]
df$gene1_chromosome[1] <- tbl1$chromosome_name[1]
df$gene2_chromosome[1] <- tbl2$chromosome_name[1]
df$gene1_start[1] <- min(tbl1$start_position)
df$gene2_start[1] <- min(tbl2$start_position)
df$gene1_end[1] <- max(tbl1$end_position)
df$gene2_end[1] <- max(tbl2$end_position)
locs.out <- rbind(locs.out,df)
}, error=function(e){writeLines(paste0("at index/step ", i, " occurred following error ", as.character(e) ))})
}
```

### Write clipped gff and bedfiles for gene arrays
Run in R

```
#Make a subdirectory for the clipped gff files
gffpath="/gffs/EnsemblVertebrates/"
dir.create(paste0(gffpath,genefam[g]))
complete_locs <- locs.out[complete.cases(locs.out[ , c(5,10)]),]
intact <- complete_locs[complete_locs$gene1_chromosome == complete_locs$gene2_chromosome,]

#clip full genome gffs by chromosome and between only the bookend genes and make new gff file and array bedfile
for(i in 1:length(intact$species)){
tryCatch({
print(paste("Working on",intact$scientific_name[i]))
chr <- intact$gene1_chromosome[i]
my_filter <- list(type=c("gene", "mRNA","CDS","exon","five_prime_UTR","three_prime_UTR"), seqid=as.character(chr))
f <- list.files(gffpath, pattern=intact$version[i]) 
gffgr <- readGFFAsGRanges(file=paste0(gffpath,f), filter = my_filter)
start <- 2 + sort(c(intact$gene1_start[i],intact$gene1_end[i],intact$gene2_start[i],intact$gene2_end[i]))[2]
end <- sort(c(intact$gene1_start[i],intact$gene1_end[i],intact$gene2_start[i],intact$gene2_end[i]))[3] - 2
query=GRanges(seqnames=as.character(chr),
          ranges=IRanges(start = start, end = end))
sub <- subsetByOverlaps(gffgr,query)
export.gff(sub, con = paste0(gffpath,genefam[g],"/",intact$version[i],"_",genefam[g],".gff3"))
write.table(cbind(chr,start,end), file = paste0(gffpath,genefam[g],"/",intact$version[i],"_cluster.bed"), 
            row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
}, error=function(e){writeLines(paste0("at index/step ", i, " occurred following error ", as.character(e) ))})
}


#get human coords and human clipped gff and bedfile
tbl <- getBM(attributes = c('hgnc_symbol',"chromosome_name","start_position","end_position"), 
      filters = "hgnc_symbol", 
      values = genes[[g]],
      mart = human)
ID <- grepl('^-?[0-9.]+$', tbl$chromosome_name)
tbl <- tbl[ID,]
tbl1 <- tbl[tbl$hgnc_symbol==genes[[g]][1],]
tbl2 <- tbl[tbl$hgnc_symbol==genes[[g]][2],]
df <- data.frame(matrix(ncol = 10, nrow = 1))
colnames(df) <- c("species","scientific_name","genome","version","gene1_chromosome","gene1_start","gene1_end","gene2_chromosome","gene2_start","gene2_end")
df$species[1] <- "Human"
df$scientific_name[1] <- "Homo sapiens"
df$genome[1] <- genomes$dataset[72]
df$version[1] <- genomes$version[72]
df$gene1_chromosome[1] <- tbl1$chromosome_name[1]
df$gene2_chromosome[1] <- tbl2$chromosome_name[1]
df$gene1_start[1] <- min(tbl1$start_position)
df$gene2_start[1] <- min(tbl2$start_position)
df$gene1_end[1] <- max(tbl1$end_position)
df$gene2_end[1] <- max(tbl2$end_position)
chr <- df$gene1_chromosome[1]
my_filter <- list(type=c("gene", "mRNA","CDS","exon","five_prime_UTR","three_prime_UTR"), seqid=as.character(chr))
f <- list.files(gffpath, pattern=df$version[1]) 
gffgr <- readGFFAsGRanges(file=paste0(gffpath,f), filter = my_filter)
start <- 2 + sort(c(df$gene1_start[1],df$gene1_end[1],df$gene2_start[1],df$gene2_end[1]))[2]
end <- sort(c(df$gene1_start[1],df$gene1_end[1],df$gene2_start[1],df$gene2_end[1]))[3] - 2
query=GRanges(seqnames=as.character(chr),ranges=IRanges(start = start, end = end))
sub <- subsetByOverlaps(gffgr,query)
export.gff(sub, con = paste0(gffpath,genefam[g],"/",df$version[1],"_",genefam[g],".gff3"))
write.table(cbind(chr,start,end), file = paste0(gffpath,genefam[g],"/",df$version[1],"_cluster.bed"), 
               row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
```


### Calculate genome stats with bedtools

```
###This script takes a genome fasta and gff annotation files and extracts the following:
###coordinates for the chromosome, intergenic, genic, exonic, intronic, and promoter regions
###GC/AT content data for the above coordinates

#######################
###set gene family info
#######################

#provide a dirlist.txt file with a \n delimited list of gene family ids correpsonding to the directories produced by the prior ensembl data pulldown:

for gene_fam in  $(cat /gffs/EnsemblVertebrates/dirlist.txt)
do

workdir=/gffs/EnsemblVertebrates/$gene_fam

genome_dir=/genomes/

genome_gc_dir=/genomes/GenomeGC

#will loop over all files in directory
#directory needs to contain fasta files paired with gff files (same prefix), and a single "species.txt"
cd $workdir
mkdir analyses
#make 'species.txt', which should contain shared file prefixes, which will then be appended to other files
ls *gff3 > species.txt
sed -i "s/_${gene_fam}.gff3//g" species.txt


for name in $(cat species.txt)
do

echo "working on $name"

bedtools sort -i ${name}_${gene_fam}.gff3 > ${name}_sorted.gff

#make bedtools genome file and chromosome bed
samtools faidx $genome_dir/*${name}*.fa 
awk -v OFS='\t' {'print $1,$2'} $genome_dir/*${name}*.fa.fai > ${name}.gn
awk 'BEGIN {OFS="\t"}; {print $1,"0",$2}' < ${name}.gn > ${name}_genome.bed

chromosome=$(awk 'NR==1{print $1}' ${name}_sorted.gff)
grep -w "^${chromosome}" ${genome_gc_dir}/${name}_genome.bed > ${name}_chromosome.bed

#get exon bed
awk '$3 == "exon"' ${name}_sorted.gff | gff2bed > tmp1
bedtools merge -i tmp1 > ${name}_exon.bed
rm tmp1

#get CDS bed
awk '$3 == "CDS"' ${name}_sorted.gff | gff2bed > tmp1
bedtools merge -i tmp1 > ${name}_CDS.bed
rm tmp1

#get gene bed
awk '$3 == "gene"' ${name}_sorted.gff | gff2bed > tmp1
bedtools merge -i tmp1 > ${name}_gene.bed
rm tmp1

#subtract gene and exon bed to get introns
bedtools subtract -a ${name}_gene.bed -b ${name}_exon.bed > ${name}_intron.bed

#subtract gene and cluster bed to get intergenic
bedtools subtract -a ${name}_cluster.bed -b ${name}_gene.bed > ${name}_intergenic.bed

#do math to get promoter region based on transcription start site
awk 'BEGIN {OFS="\t"}; {
if($3 == "five_prime_UTR" && $7 == "+")
{
print $1,$4 - 751,$4 + 249
}
else if($3 == "five_prime_UTR" && $7 == "-")
{
print $1,$5 - 251,$5 + 749
}
}' analyses/${name}_sorted.gff > ${name}_promoter.bed





#Get GC and AT values for each interval

bedtools nuc -fi $genome_dir/*${name}*.fa -bed ${name}_genome.bed | awk '{print $0,FS,a}' FS='\t' a="$name" > ${name}_genome_GC.tsv
cp ${genome_gc_dir}/${name}_genome_GC.tsv .

bedtools nuc -fi $genome_dir/*${name}*.fa -bed ${name}_cluster.bed | awk '{print $0,FS,a}' FS='\t' a="$name" > ${name}_cluster_GC.tsv
sed -i 's/ \+//g' ${name}_cluster_GC.tsv

bedtools nuc -fi $genome_dir/*${name}*.fa -bed ${name}_intergenic.bed | awk '{print $0,FS,a}' FS='\t' a="$name" > ${name}_intergenic_GC.tsv
sed -i 's/ \+//g' ${name}_intergenic_GC.tsv

bedtools nuc -fi $genome_dir/*${name}*.fa -bed ${name}_gene.bed | awk '{print $0,FS,a}' FS='\t' a="$name" > ${name}_gene_GC.tsv
sed -i 's/ \+//g' ${name}_gene_GC.tsv

bedtools nuc -fi $genome_dir/*${name}*.fa -bed ${name}_exon.bed | awk '{print $0,FS,a}' FS='\t' a="$name" > ${name}_exon_GC.tsv
sed -i 's/ \+//g' ${name}_exon_GC.tsv

bedtools nuc -fi $genome_dir/*${name}*.fa -bed ${name}_CDS.bed | awk '{print $0,FS,a}' FS='\t' a="$name" > ${name}_CDS_GC.tsv
sed -i 's/ \+//g' ${name}_CDS_GC.tsv

bedtools nuc -fi $genome_dir/*${name}*.fa -bed ${name}_intron.bed | awk '{print $0,FS,a}' FS='\t' a="$name" > ${name}_intron_GC.tsv
sed -i 's/ \+//g' ${name}_intron_GC.tsv

bedtools nuc -fi $genome_dir/*${name}*.fa -bed ${name}_chromosome.bed | awk '{print $0,FS,a}' FS='\t' a="$name" > ${name}_chromosome_GC.tsv
sed -i 's/ \+//g' ${name}_chromosome_GC.tsv

bedtools nuc -fi $genome_dir/*${name}*.fa -bed ${name}_promoter.bed | awk '{print $0,FS,a}' FS='\t' a="$name" > ${name}_promoter_GC.tsv
sed -i 's/ \+//g' ${name}_promoter_GC.tsv

sed -i 's/ \+//g' ${name}_genome_GC.tsv

done


mv *tsv analyses/
mv *gff analyses/
mv *gn analyses/
mv *bed analyses/

#make a compressed zip file of the results for export from server to your laptop for R graphing
dir=$(basename "$PWD")
mkdir analyses/tsvs
cp analyses/*tsv analyses/tsvs
zip -r $dir analyses/tsvs
cat dirlist.txt | parallel "zip -r {} {}"

```




