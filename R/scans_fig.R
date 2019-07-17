# genome scans by H-K, LMM, LMM w/ LOCO

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
    out_hk <- scan1(pr, pheno=arab$pheno, cores=0)
    operm_hk <- scan1perm(pr, pheno=arab$pheno[,"fruit_length"], n_perm=1000, cores=0)

    k <- calc_kinship(pr, cores=0)
    out_lmm <- scan1(pr, arab$pheno, k, cores=0)

    kloco <- calc_kinship(pr, "loco", cores=0)
    out_loco <- scan1(pr, arab$pheno, kloco, cores=0)

    save(out_hk, operm_hk, out_lmm, out_loco, file=file)
}

plot_scan <-
    function(version=1, lodcolumn="fruit_length", title="fruit length",
             color=c("darkslateblue", "#49A56E", "#F15A24"), legend=TRUE,
             ymx=max(cbind(out_hk[,lodcolumn], out_lmm[,lodcolumn], out_loco[,lodcolumn])))
{
    par(mar=c(4.1,4.1,0.6,0.6),
        fg="white", col.lab="white", col.axis="white")

    plot(out_hk, pmap, lodcolumn=lodcolumn, ylim=c(0, ymx*1.05), col=color[1])

    if(version==1) abline(h=summary(operm_hk, 0.05), lty=2, col=color[1])

    u <- par("usr")
    text(u[2]-diff(u[1:2])*0.015, u[4]-diff(u[3:4])*0.02, title, col="black", adj=c(1, 1))

    if(version>=2) plot(out_lmm, pmap, lodcolumn=lodcolumn, col=color[2], add=TRUE)
    if(version>=3) plot(out_loco, pmap, lodcolumn=lodcolumn, col=color[3], add=TRUE)

    if(version>=2 && legend) {
        legend("topleft", col=color[1:version], lwd=2, c("haley-knott", "lmm", "lmm w/loco")[1:version],
               text.col="black", box.col="black", bg="gray90")
    }
}

pdf("../Figs/scan_hk.pdf", height=5.5, width=11, pointsize=16)
plot_scan(1, ymx=41.55)
dev.off()

pdf("../Figs/scan_lmm.pdf", height=5.5, width=11, pointsize=16)
plot_scan(2, ymx=41.55)
dev.off()

pdf("../Figs/scan_loco.pdf", height=5.5, width=11, pointsize=16)
plot_scan(3, ymx=41.55)
dev.off()

pdf("../Figs/scan_seedwt.pdf", height=5.5, width=11, pointsize=16)
plot_scan(3, "seed_weight", "seed weight", ymx=41.55)
dev.off()
