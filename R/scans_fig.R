library(qtl2)

url <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/ArabMAGIC/"

arab <- read_cross2(paste0(url, "arabmagic_tair9.zip"))

gmap <- insert_pseudomarkers(arab$gmap, step=0.2, stepwidth="max")
pmap <- interp_map(gmap, arab$gmap, arab$pmap)

file <- "_cache/scans.RData"
if(file.exists(file)) {
    load(file)
} else {
    pr <- calc_genoprob(arab, gmap, error_prob=0.002, cores=0)

    set.seed(13914402)
    out_hk <- scan1(pr, pheno=arab$pheno[,"fruit_length"], cores=0)
    operm_hk <- scan1perm(pr, pheno=arab$pheno[,"fruit_length"], n_perm=1000, cores=0)

    k <- calc_kinship(pr, cores=0)
    out_lmm <- scan1(pr, arab$pheno[,"fruit_length"], k, cores=0)
#    operm_lmm <- scan1perm(pr, arab$pheno[,"fruit_length"], k, n_perm=1000, cores=0)

    kloco <- calc_kinship(pr, "loco", cores=0)
    out_loco <- scan1(pr, arab$pheno[,"fruit_length"], kloco, cores=0)
#    operm_loco <- scan1perm(pr, arab$pheno[,"fruit_length"], kloco, n_perm=1000, cores=0)

    save(out_hk, operm_hk,
         out_lmm, #operm_lmm,
         out_loco, #operm_loco,
         file=file)
}

pdf("../Figs/scan_hk.pdf", height=5.5, width=11, pointsize=16)
par(mar=c(4.1,4.1,0.6,0.6),
    fg="white", col.lab="white", col.axis="white")
plot(out_hk, pmap, ylim=c(0, max(out_loco)*1.05))
abline(h=summary(operm_hk, 0.05), lty=2,col="darkslateblue")
u <- par("usr")
text(u[2]-diff(u[1:2])*0.02, u[4]-diff(u[1:2])*0.02, "Fruit length", col="black", adj=c(1, 0.5))
dev.off()
