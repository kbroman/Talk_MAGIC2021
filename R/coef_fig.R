# plot coefficients, simple and BLUP

library(qtl2)
library(broman)
library(glue)

url <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/ArabMAGIC/"

arab <- read_cross2(paste0(url, "arabmagic_tair9.zip"))

file <- "_cache/coef.RData"
if(file.exists(file)) {
    load(file)
} else {
    pr <- calc_genoprob(arab, error_prob=0.002, cores=0)

    load("_cache/scans.RData")

    fl_peak <- max(out_hk, arab$pmap, lodcolumn="fruit_length")

    peaks <- find_peaks(out_hk, arab$pmap, threshold=20, peakdrop=5)
    sw_peak <- subset(peaks, lodcolumn=="seed_weight" & chr==1 & pos > 20)

    fl_pr <- pull_genoprobpos(pr, arab$pmap, fl_peak$chr, fl_peak$pos)
    fit1_fl <- fit1(fl_pr, arab$pheno[,"fruit_length"])
    blup_fl <- fit1(fl_pr, arab$pheno[,"fruit_length"], blup=TRUE)

    sw_pr <- pull_genoprobpos(pr, arab$pmap, sw_peak$chr, sw_peak$pos)
    fit1_sw <- fit1(sw_pr, arab$pheno[,"seed_weight"])
    blup_sw <- fit1(sw_pr, arab$pheno[,"seed_weight"], blup=TRUE)
    # LL is quite different, because there are no lines with LL genotype

    m_fl <- maxmarg(pr, arab$pmap, chr=fl_peak$chr, pos=fl_peak$pos, return_char=TRUE)
    m_sw <- maxmarg(pr, arab$pmap, chr=sw_peak$chr, pos=sw_peak$pos, return_char=TRUE)

    save(fl_peak, fit1_fl, blup_fl,
         sw_peak, fit1_sw, blup_sw,
         m_fl, m_sw, file=file)
}


fit1_fl_lo <- (fit1_fl$coef - 2*fit1_fl$SE)[1:19]
fit1_fl_hi <- (fit1_fl$coef + 2*fit1_fl$SE)[1:19]

blup_fl_lo <- (blup_fl$coef - 2*blup_fl$SE)[1:19]
blup_fl_hi <- (blup_fl$coef + 2*blup_fl$SE)[1:19]

fit1_sw_lo <- (fit1_sw$coef - 2*fit1_sw$SE)[1:19]
fit1_sw_hi <- (fit1_sw$coef + 2*fit1_sw$SE)[1:19]

blup_sw_lo <- (blup_sw$coef - 2*blup_sw$SE)[1:19]
blup_sw_hi <- (blup_sw$coef + 2*blup_sw$SE)[1:19]


xd <- 0.2


pdf("../Figs/coef_fl.pdf", height=5.5, width=11, pointsize=16)
par(mar=c(2.6, 3.1, 2.6, 0.6))
par(col.lab="white", col.axis="white", col.main="white")
ciplot(fit1_fl$coef[1:19], fit1_fl$SE[1:19], labels=rep("", 19), ci_endseg=0.1,
       ylab="", xlab="", mgp=c(1.8,0.3,0))
main <- glue("Fruit length (chr {chr} @ {pos} Mbp)", chr=fl_peak$chr, pos=myround(fl_peak$pos, 1))
title(ylab="QTL effect", , mgp=c(1.8,0.3,0), main=main)
axis(side=1, at=1:19, rownames(arab$founder_geno[[1]]), las=3, tick=FALSE, mgp=c(0, 0.3, 0))
dev.off()

pdf("../Figs/blup_fl.pdf", height=5.5, width=11, pointsize=16)
par(mar=c(2.6, 3.1, 2.6, 0.6))
par(col.lab="white", col.axis="white", col.main="white")
ciplot(fit1_fl$coef[1:19], fit1_fl$SE[1:19], labels=rep("", 19), ci_endseg=0.1,
       ylab="", xlab="", mgp=c(1.8,0.3,0))
title(ylab="QTL effect", , mgp=c(1.8,0.3,0), main=main)
axis(side=1, at=1:19, rownames(arab$founder_geno[[1]]), las=3, tick=FALSE, mgp=c(0, 0.3, 0))

segments(1:19+xd, blup_fl_lo, 1:19+xd, blup_fl_hi, lwd=2)
points(1:19+xd, blup_fl$coef[1:19], pch=21, bg="violetred")
segments(1:19+xd-0.1, blup_fl_lo, 1:19+xd+0.1, blup_fl_lo, lwd=2)
segments(1:19+xd-0.1, blup_fl_hi, 1:19+xd+0.1, blup_fl_hi, lwd=2)
legend("bottomright", lwd=2, pch=21, pt.bg=c("slateblue", "violetred"), c("least squares", "BLUP"),
       bg="gray90")
dev.off()

pdf("../Figs/blup_sw.pdf", height=5.5, width=11, pointsize=16)
par(mar=c(2.6, 3.1, 2.6, 0.6))
par(col.lab="white", col.axis="white", col.main="white")
ciplot(fit1_sw$coef[1:19], fit1_sw$SE[1:19], labels=rep("", 19), ci_endseg=0.1,
       ylab="", xlab="", mgp=c(1.8,0.3,0), ylim=c(-4.3, 4.3))
main <- glue("Seed weight (chr {chr} @ {pos} Mbp)", chr=sw_peak$chr, pos=myround(sw_peak$pos, 1))
title(ylab="QTL effect", , mgp=c(1.8,0.3,0), main=main)
axis(side=1, at=1:19, rownames(arab$founder_geno[[1]]), las=3, tick=FALSE, mgp=c(0, 0.3, 0))

segments(1:19+xd, blup_sw_lo, 1:19+xd, blup_sw_hi, lwd=2)
points(1:19+xd, blup_sw$coef[1:19], pch=21, bg="violetred")
segments(1:19+xd-0.1, blup_sw_lo, 1:19+xd+0.1, blup_sw_lo, lwd=2)
segments(1:19+xd-0.1, blup_sw_hi, 1:19+xd+0.1, blup_sw_hi, lwd=2)
legend("topright", lwd=2, pch=21, pt.bg=c("slateblue", "violetred"), c("least squares", "BLUP"),
       bg="gray90")
dev.off()
