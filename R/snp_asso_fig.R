# SNP association scan

library(qtl2)

url <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/ArabMAGIC/"
url <- "~/Code/Rqtl2/qtl2data/ArabMAGIC/"

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
             ymx=max(out_hk[,"fruit_length"]), logp=FALSE)
{
    par(mar=c(4.1,4.1,0.6,0.6), fg="white", col.lab="white", col.axis="white")

    green <- "#49A56E"

    if(logp) {
        ylab <- expression(paste(-log[10], "  p-value"))
        out_snps[,lodcolumn] <- -pchisq(out_snps[,lodcolumn]*2*log(10),
                                        1, log=TRUE, lower=FALSE)/log(10)
        out_hk[,lodcolumn] <- -pchisq(out_hk[,lodcolumn]*2*log(10),
                                      18, log=TRUE, lower=FALSE)/log(10)
        ymx <- max(c(out_snps[,lodcolumn], out_hk[,lodcolumn]))
    } else {
        ylab <- "LOD score"
    }

    plot(out_snps, arab$pmap, lod=lodcolumn, type="p", pch=16, altcol=green, gap=0,
         ylim=c(0, ymx*1.05), cex=0.6, ylab=ylab)

    if(version==2) {
        plot(out_hk, arab$pmap, lod=lodcolumn, gap=0, altcol=green, add=TRUE)
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

pdf("../Figs/snp_asso_B_logp.pdf", height=5.5, width=11, pointsize=16)
plot_snp_asso(2, logp=TRUE)
dev.off()

pdf("../Figs/snp_asso_C_logp.pdf", height=5.5, width=11, pointsize=16)
plot_snp_asso(2, "seed_weight", "seed weight", logp=TRUE)
dev.off()
