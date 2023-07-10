# TandemClipR


This code is used in the manuscript 'Emergence and influence of sequence bias in evolutionarily malleable, mammalian tandem
arrays'
Authors: Margarita V Brovkina; Margaret A. Chapman; Matthew L. Holding; E. Josephine Clowney

The following code, run with a combination of command-line tools and R, will:
 1. Download ensembl genome fastas and corresponding gffs for full Ensembl vertebrate database
 2. Tabulate locations for human genes and their orthologs based on a list, where the named genes bookend a focal region of the genome
 3. Calculate GC content for defined regions and subregions (genic, intergenic, exonic, intronic, promoter)
 4. Aggregate GC content and gene count data for defined regions and subregions across species for downstream statistical analyses



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
* [Collating GC and Gene Count Data](#collating-gc-and-gene-count-data)
* [Get Promoter GC data](#get-promoter-gc-data)


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

####################################################################################################################
#Now we will link chromosome locations for the human bookends with the ortholog locations in other species genomes using biomaRt getLDS() function
####################################################################################################################

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

####################################################################################################################
#clip full genome gffs by chromosome and between only the bookend genes and make new gff file and array bedfile
####################################################################################################################

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

####################################################################################################################
#get human coords and human clipped gff and bedfile with getBM() 
####################################################################################################################

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
Run in command-line with bedtools and awk installed

```
####################################################################################################################
###This script takes a genome fasta and gff annotation files and extracts the following:
###coordinates for the chromosome, intergenic, genic, exonic, intronic, and promoter regions
###GC/AT content data for the above coordinates
####################################################################################################################


####################################################################################################################
###set gene family info
####################################################################################################################

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




####################################################################################################################
#Get GC and AT values for each interval
####################################################################################################################

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

####################################################################################################################
### make a compressed zip file of the results for export from server to your laptop for R graphing
####################################################################################################################
dir=$(basename "$PWD")
mkdir analyses/tsvs
cp analyses/*tsv analyses/tsvs
zip -r $dir analyses/tsvs
cat dirlist.txt | parallel "zip -r {} {}"

```


### Collating GC and Gene Count Data
Run in R

```
#set up target family list and set g# to access one family
genefam=c("KRTAP","HOXA","KRTAP_2","THEM_LCE","SERPINA", "CYP2", "HBB_OR")

#genefam=c("KRTAP","HOXA","HOXD","KRTAP_2","IL1R","GABR","MYH","THEM_LCE","SERPINA", "CYP2", "HBB_OR")
for(g in 1:length(genefam)){
  ## get list of species' names
  genomes <- read.table("~/Desktop/families/genomes.txt", header=TRUE, sep = "\t")
  list <- sub('_chromosome_GC.tsv$', '', list.files(paste0("~/Desktop/families/",genefam[g],"/analyses/"), pattern="chromosome_GC.tsv"))
  print(paste("There are",length(list),"species for",genefam[g]))


####################################################################################################################
  #### Collate results of ```bedtools nuc``` runs for chromosome, intergenic regions, genes, exons, and introns
####################################################################################################################
  ## read in whole genome GC data
  rm(genome_data)
  genome_data <- data.frame(matrix(ncol = 10, nrow = 0))
  
  Chr <- 1
  Start <- 2
  End <- 3
  Acol <- 6
  Ccol <- 7
  Gcol <- 8
  Tcol <- 9
  Ncol <- 10
  length.col <- 12
  version.col <- 13 
  
  for(i in 1:length(list)) {
    tryCatch({
      table <- read.table(paste0("~/Desktop/families/",genefam[g],"/analyses/",list[i],"_genome_GC.tsv"))
      table <- table[,c(Chr,Start,End,Acol,Ccol,Gcol,Tcol,Ncol,length.col,version.col)]
      median_genome <- data.frame(t(data.frame(apply(table[,c(4:9)],2,sum))))
      median_genome$Version <- table$V13[1]
      median_genome$feature <-"genome"
      rownames(median_genome) <- i
      init <- data.frame(t(data.frame(c("NA","NA","NA"))))
      colnames(init) <- c("Chr", "Start", "End")
      rownames(init) <- i
      median_genome <- cbind(init,median_genome)
      genome_data <- rbind(genome_data,median_genome)
    }, error=function(e){writeLines(paste0("at index/step ", i, " occurred following error ", as.character(e) ))})
  }
  
  colnames(genome_data) <- c("Chr", "Start", "End","A.count","C.count","G.count","T.count","N.count","Length","Version","feature")
  
  
  ## read in whole chromosome GC data
  rm(chromosome_data)
  chromosome_data <- data.frame(matrix(ncol = 10, nrow = 0))
  
  Chr <- 1
  Start <- 2
  End <- 3
  Acol <- 6
  Ccol <- 7
  Gcol <- 8
  Tcol <- 9
  Ncol <- 10
  length.col <- 12
  version.col <- 13 
  
  for(i in 1:length(list)) {
    tryCatch({
      table <- read.table(paste0("~/Desktop/families/",genefam[g],"/analyses/",list[i],"_chromosome_GC.tsv"))
      table <- table[,c(Chr,Start,End,Acol,Ccol,Gcol,Tcol,Ncol,length.col,version.col)]
      chromosome_data <- rbind(chromosome_data,table)
    }, error=function(e){writeLines(paste0("at index/step ", i, " occurred following error ", as.character(e) ))})
  }
  colnames(chromosome_data) <- c("Chr", "Start", "End","A.count","C.count","G.count","T.count","N.count","Length","Version")
  chromosome_data$feature <- "chromosome"
  
  
  ## REMOVE BAD GENOMES
  thediff <- setdiff(genome_data$Version,chromosome_data$Version)
  ID <- genome_data$Version %in% thediff
  genome_data <- genome_data[!ID,]
  rownames(genome_data) <- c(1:length(genome_data$Chr))
  
  
  ## read in whole intergenic GC data
  rm(intergenic_data)
  intergenic_data <- data.frame(matrix(ncol = 10, nrow = 0))
  
  Chr <- 1
  Start <- 2
  End <- 3
  Acol <- 6
  Ccol <- 7
  Gcol <- 8
  Tcol <- 9
  Ncol <- 10
  length.col <- 12
  version.col <- 13 
  
  for(i in 1:length(list)) {
    tryCatch({
      table <- read.table(paste0("~/Desktop/families/",genefam[g],"/analyses/",list[i],"_intergenic_GC.tsv"))
      table <- table[,c(Chr,Start,End,Acol,Ccol,Gcol,Tcol,Ncol,length.col,version.col)]
      #chr_minus_cluster <- table[c(1,nrow(table)),] #first and last rows are chromosome wings outside cluster
      #table <- table[c(2:(nrow(table)-1)),]
      intergenic_data <- rbind(intergenic_data,table)
    }, error=function(e){writeLines(paste0("at index/step ", i, " occurred following error ", as.character(e) ))})
  }
  
  colnames(intergenic_data) <- c("Chr", "Start", "End","A.count","C.count","G.count","T.count","N.count","Length","Version")
  intergenic_data$feature <- "intergenic"
  intergenic_data <- intergenic_data[complete.cases(intergenic_data[ ,10]),]
  
  
  ## read in whole cluster GC data
  rm(cluster_data)
  cluster_data <- data.frame(matrix(ncol = 10, nrow = 0))
  
  Chr <- 1
  Start <- 2
  End <- 3
  Acol <- 6
  Ccol <- 7
  Gcol <- 8
  Tcol <- 9
  Ncol <- 10
  length.col <- 12
  version.col <- 13 
  
  for(i in 1:length(list)) {
    tryCatch({
      table <- read.table(paste0("~/Desktop/families/",genefam[g],"/analyses/",list[i],"_cluster_GC.tsv"))
      table <- table[,c(Chr,Start,End,Acol,Ccol,Gcol,Tcol,Ncol,length.col,version.col)]
      cluster_data <- rbind(cluster_data,table)
    }, error=function(e){writeLines(paste0("at index/step ", i, " occurred following error ", as.character(e) ))})
  }
  colnames(cluster_data) <- c("Chr", "Start", "End","A.count","C.count","G.count","T.count","N.count","Length","Version")
  cluster_data$feature <- "cluster"
  
  
  
  ## read in gene GC data
  rm(gene_data)
  gene_data <- data.frame(matrix(ncol = 10, nrow = 0))
  
  Chr <- 1
  Start <- 2
  End <- 3
  Acol <- 6
  Ccol <- 7
  Gcol <- 8
  Tcol <- 9
  Ncol <- 10
  length.col <- 12
  version.col <- 13
  
  for(i in 1:length(list)) {
    tryCatch({
      table <- read.table(paste0("~/Desktop/families/",genefam[g],"/analyses/",list[i],"_gene_GC.tsv"), sep = "\t")
      table <- table[,c(Chr,Start,End,Acol,Ccol,Gcol,Tcol,Ncol,length.col,version.col)]
      gene_data <- rbind(gene_data,table)
    }, error=function(e){writeLines(paste0("at index/step ", i, " occurred following error ", as.character(e) ))})
  }
  colnames(gene_data) <- c("Chr", "Start", "End","A.count","C.count","G.count","T.count","N.count","Length","Version")
  gene_data$feature <- "gene"
  
  
  
  ## read in CDS GC data
  rm(CDS_data)
  CDS_data <- data.frame(matrix(ncol = 10, nrow = 0))
  
  Chr <- 1
  Start <- 2
  End <- 3
  Acol <- 6
  Ccol <- 7
  Gcol <- 8
  Tcol <- 9
  Ncol <- 10
  length.col <- 12
  version.col <- 13 
  
  for(i in 1:length(list)) {
    tryCatch({
      table <- read.table(paste0("~/Desktop/families/",genefam[g],"/analyses/",list[i],"_CDS_GC.tsv"), sep = "\t")
      table <- table[,c(Chr,Start,End,Acol,Ccol,Gcol,Tcol,Ncol,length.col,version.col)]
      CDS_data <- rbind(CDS_data,table)
    }, error=function(e){writeLines(paste0("at index/step ", i, " occurred following error ", as.character(e) ))})
  }
  colnames(CDS_data) <- c("Chr", "Start", "End","A.count","C.count","G.count","T.count","N.count","Length","Version")
  CDS_data$feature <- "CDS"
  
  
  
  ## read in intron GC data
  rm(intron_data)
  intron_data <- data.frame(matrix(ncol = 10, nrow = 0))
  
  Chr <- 1
  Start <- 2
  End <- 3
  Acol <- 6
  Ccol <- 7
  Gcol <- 8
  Tcol <- 9
  Ncol <- 10
  length.col <- 12
  version.col <- 13 
  
  for(i in 1:length(list)) {
    tryCatch({
      table <- read.table(paste0("~/Desktop/families/",genefam[g],"/analyses/",list[i],"_intron_GC.tsv"), sep = "\t")
      table <- table[,c(Chr,Start,End,Acol,Ccol,Gcol,Tcol,Ncol,length.col,version.col)]
      intron_data <- rbind(intron_data,table)
    }, error=function(e){writeLines(paste0("at index/step ", i, " occurred following error ", as.character(e) ))})
  }
  colnames(intron_data) <- c("Chr", "Start", "End","A.count","C.count","G.count","T.count","N.count","Length","Version")
  intron_data$feature <- "intron"
  

####################################################################################################################
  #### Summarize into species' counts and GC-content for plotting purposes
####################################################################################################################

  ## Aggregate nucleotide counts by species and annotation type, and calculate GC%
  data <- rbind(genome_data,chromosome_data,intergenic_data,cluster_data,gene_data,CDS_data,intron_data)
  df <- aggregate(cbind(A.count, C.count,G.count,T.count,N.count,Length) ~ Version + feature, data = data, FUN = sum, na.rm = FALSE)
  df$GC <-(df$G.count + df$C.count)/(df$Length-df$N.count)
  
  #test <- (plyr::count(data, c("feature", "Version")))
  #test$Vfeat <- paste0(test$feature,"_",test$Version)
  #df$Vfeat <- paste0(df$feature,"_",df$Version)
  #setdiff(test$Vfeat,df$Vfeat)
  
  df$feature.count <- (plyr::count(data, c("feature", "Version")))[,3] 
  
  ## have a look at the head of df and gene features to check that it came together correctly
  print(head(subset(df, feature=="gene")))
  
  ## Get additional focal metadata comparisons for graphs
  for(i in 1:length(df$Version)){
    df$species[i] <- genomes$species[genomes$version==df$Version[i]]
  }
  
  for(i in 1:length(df$Version)){
    df$common_name[i] <- genomes$Common.name[genomes$version==df$Version[i]]
  }
  
  for(i in 1:length(df$Version)){
    df$class[i] <- genomes$Class[genomes$version==df$Version[i]]
  }
  
  intergenic_check <- unique(sort(intergenic_data$Version))
  df <- df[df$Version %in% intergenic_check,]
  
  genic_check <- unique(sort(gene_data$Version))
  df <- df[df$Version %in% genic_check,]
  
  CDS_check <- unique(sort(CDS_data$Version))
  df <- df[df$Version %in% CDS_check,]
  
  intronic_check <- unique(sort(intron_data$Version))
  df <- df[df$Version %in% intronic_check,]
  
  chromosomes <- subset(df,feature=="chromosome")
  species <- chromosomes$species
  class <- chromosomes$class
  clusters <- subset(df,feature=="cluster")
  genome <- subset(df,feature=="genome")
  thediff <- setdiff(genome$Version,clusters$Version)
  ID <- genome$Version %in% thediff
  genome <- genome[!ID,]
  genes <-  subset(df,feature=="gene")
  GenGC <- genome$GC
  ClusterGC <- clusters$GC
  GenGC_vs_ClustGC <- clusters$GC - genome$GC 
  length_kb <- (clusters$Length - clusters$N.count)/1000
  introns <- subset(df,feature=="intron")
  CDSs <- subset(df,feature=="CDS")
  intergenic <- subset(df,feature=="intergenic")
  genicGC <- genes$GC
  CDSGC <- CDSs$GC
  intronicGC <- introns$GC
  intergenicGC <- intergenic$GC
  
  plot_data <- as.data.frame(cbind(species,class,genes$feature.count,GenGC,ClusterGC,GenGC_vs_ClustGC, length_kb, genicGC, CDSGC, intronicGC, intergenicGC))
  
  colnames(plot_data) <- c("Species","Class", "Gene.count","GenGC", "ClusterGC","GenGC_vs_ClustGC","LengthKB","genicGC", "CDSGC", "intronicGC", "intergenicGC")
  
  plot_data$Species <- factor(plot_data$Species)
  plot_data$Class <- factor(plot_data$Class)
  plot_data$Gene.count <- as.numeric(plot_data$Gene.count)
  plot_data$logGenes <- log(plot_data$Gene.count)
  plot_data$GenGC_vs_ClustGC <- as.numeric(plot_data$GenGC_vs_ClustGC)
  plot_data$GenGC <- as.numeric(plot_data$GenGC)
  plot_data$ClusterGC <- as.numeric(plot_data$ClusterGC)
  plot_data$LengthKB <- as.numeric(plot_data$LengthKB)
  plot_data$genicGC <- as.numeric(plot_data$genicGC)
  plot_data$CDSGC <- as.numeric(plot_data$CDSGC)
  plot_data$intronicGC <- as.numeric(plot_data$intronicGC)
  plot_data$intergenicGC <- as.numeric(plot_data$intergenicGC)

####################################################################################################################
### write aggregated plot data for each genefam
####################################################################################################################
  write.table(plot_data, file = paste0("~/Desktop/families/plots/",genefam[g],"_data.txt"), sep="\t", quote = FALSE,row.names = FALSE)
  
```


### Get Promoter GC data
Run in R

```

####################################################################################################################
### repeats loop and aggregation as above for other regions, but for subset of genomes with annotated transcription start sites TSS
####################################################################################################################

for(g in 1:length(genefam)){
  
  ## get list of species' names
  genomes <- read.table("~/Desktop/families/genomes.txt", header=TRUE, sep = "\t")
  list <- sub('_promoter_GC.tsv$', '', list.files(paste0("~/Desktop/families/",genefam[g],"/analyses/"), pattern="promoter_GC.tsv"))
  print(paste("There are",length(list),"species for",genefam[g]))
  
  ## read in promoter GC data
  rm(promoter_data)
  promoter_data <- data.frame(matrix(ncol = 10, nrow = 0))
  
  Chr <- 1
  Start <- 2
  End <- 3
  Acol <- 6
  Ccol <- 7
  Gcol <- 8
  Tcol <- 9
  Ncol <- 10
  length.col <- 12
  version.col <- 13 
  
  for(i in 1:length(list)) {
    tryCatch({
      table <- read.table(paste0("~/Desktop/families/",genefam[g],"/analyses/",list[i],"_promoter_GC.tsv"), sep = "\t")
      table <- table[,c(Chr,Start,End,Acol,Ccol,Gcol,Tcol,Ncol,length.col,version.col)]
      promoter_data <- rbind(promoter_data,table)
    }, error=function(e){writeLines(paste0("at index/step ", i, " occurred following error ", as.character(e) ))})
  }
  colnames(promoter_data) <- c("Chr", "Start", "End","A.count","C.count","G.count","T.count","N.count","Length","Version")
  promoter_data$feature <- "promoter"
  
  #prep tss data
  df <- aggregate(cbind(A.count, C.count,G.count,T.count,N.count,Length) ~ Version + feature, data = promoter_data, FUN = sum, na.rm = FALSE)
  df$promoterGC <-(df$G.count + df$C.count)/(df$Length-df$N.count)
  
  ## Get additional focal metadata comparisons for graphs
  for(i in 1:length(df$Version)){
    df$species[i] <- genomes$species[genomes$version==df$Version[i]]
  }
  
  plot_data <- read.table(paste0("~/Desktop/families/plots/",genefam[g],"_data.txt"),header = TRUE,
                          sep = "\t")
  
  promoter_check <- df$species
  plot_data <- plot_data[plot_data$Species %in% promoter_check,]
  
  for(i in 1:length(plot_data$Species)){
    plot_data$promoterGC[i] <- df$promoterGC[df$species==plot_data$Species[i]]
  }
  
  write.table(plot_data, file = paste0("~/Desktop/families/plots/",genefam[g],"_promoter_reduced_data.txt"), sep="\t", quote = FALSE,row.names = FALSE)
  ```



