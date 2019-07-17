# SNP association scan

library(qtl2)

url <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/ArabMAGIC/"

arab <- read_cross2(paste0(url, "arabmagic_tair9.zip"))
arab <- drop_nullmarkers(arab)

file <- "_cache/snp_asso.rds"
if(file.exists(file)) {
    out_snps <- readRDS(file)
} else {
    pr <- calc_genoprob(arab, error_prob=0.002, cores=0)

    snp_pr <- genoprob_to_snpprob(pr, arab)

    out_snps <- scan1(snp_pr, arab$pheno, cores=0)

    saveRDS(out_snps, file)
}

load("_cache/scans.RData")


plot_snp_asso <-
    function(version=1, lodcolumn="fruit_length", title="fruit length",
             ymx=max(out_hk[,"fruit_length"]))
{
    par(mar=c(4.1,4.1,0.6,0.6),
    fg="white", col.lab="white", col.axis="white")

    plot(out_snps, arab$pmap, lod=lodcolumn, type="p", pch=16, altcol="green3", gap=0,
         ylim=c(0, ymx*1.05), cex=0.6)

    if(version==2) {
        plot(out_hk, arab$pmap, lod=lodcolumn, gap=0, altcol="green3", add=TRUE)
    }

    u <- par("usr")
    text(u[2]-diff(u[1:2])*0.015, u[4]-diff(u[3:4])*0.02, title, col="black", adj=c(1, 1))
}



pdf("../Figs/snp_asso.pdf", height=5.5, width=11, pointsize=16)
plot_snp_asso()
dev.off()

pdf("../Figs/snp_asso_B.pdf", height=5.5, width=11, pointsize=16)
plot_snp_asso(2)
dev.off()

pdf("../Figs/snp_asso_C.pdf", height=5.5, width=11, pointsize=16)
plot_snp_asso(2, "seed_weight", "seed weight")
dev.off()
