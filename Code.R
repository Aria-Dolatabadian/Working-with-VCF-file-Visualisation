#import files

pkg <- "pinfsc50"
vcf_file <- ("pinf_sc50.vcf.gz")
dna_file <- ("pinf_sc50.fasta")
gff_file <- ("pinf_sc50.gff")

#Then read in the VCF file with vcfR.


library(vcfR)
vcf <- read.vcfR( vcf_file, verbose = FALSE )


dna <- ape::read.dna(dna_file, format = "fasta")

gff <- read.table(gff_file, sep="\t", quote="")

#Creating chromR objects

library(vcfR)
chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)

#Processing chromR objects

plot(chrom)

chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, min_MQ = 59.9,  max_MQ = 60.1)
plot(chrom)

chrom <- proc.chromR(chrom, verbose=TRUE)

plot(chrom)

#Visualizing data

chromoqc(chrom, dp.alpha=20)

chromoqc(chrom, xlim=c(5e+05, 6e+05))

#https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html
