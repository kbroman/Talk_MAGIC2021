# figures showing data files

library(broman)
library(data.table)

url <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/ArabMAGIC/"

pfile <- paste0(url, "arabmagic_pheno.csv")
gfile <- paste0(url, "arabmagic_geno.csv")
fgfile <- paste0(url, "arabmagic_foundergeno.csv")
pmapfile <- paste0(url, "arabmagic_pmap_tair9.csv")

# phenotype file
phe <- fread(pfile, data.table=FALSE, skip=3)
phe[is.na(phe)] <- "NA"
excel_fig(phe[1:15,1:7], "../Figs/phefile.pdf", cellwidth=110)

# genotype file
gen <- fread(gfile, data.table=FALSE, skip=3)
excel_fig(gen[1:15,1:7], "../Figs/genfile.pdf", cellwidth=110)

# founder genotype file
fg <- fread(fgfile, data.table=FALSE, skip=3)
excel_fig(fg[1:21,], "../Figs/fgfile.pdf", cellwidth=c(80, 110, rep(50, 19)))

# physical map file
pmap <- fread(pmapfile, data.table=FALSE, skip=3)
excel_fig(pmap[1:19,], "../Figs/pmapfile.pdf", cellwidth=110)
