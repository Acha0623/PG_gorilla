library(snpMatrix)

#Gbb
gbb <- read.plink("Gbb_block")
ldgbb <- ld.snp(gbb, dep=791)
plot.snp.dprime(ldgbb, filename="Gbb_LDblock.eps", res=100)
#Gbg
gbg <- read.plink("Gbg_block")
ldgbg <- ld.snp(gbg, dep=1286)
plot.snp.dprime(ldgbg, filename="Gbg_LDblock.eps", res=100)
#Ggg
ggg <- read.plink("Ggg_block")
ldggg <- ld.snp(ggg, dep=2166)
plot.snp.dprime(ldggg, filename="Ggg_LDblock.eps", res=100)

print("Done LDblock plotting")
