# plot coefficients, simple and BLUP

url <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/ArabMAGIC/"

arab <- read_cross2(paste0(url, "arabmagic_tair9.zip"))

file <- "_cache/coef.RData"
if(file.exists(file)) {
    load(file)
} else {
    pr <- calc_genoprob(arab, error_prob=0.002, cores=0)

    load("_cache/scans.RData")
    peaks <- find_peaks(out_hk, arab$pmap, threshold=20, peakdrop=5)

    fl_peak <- subset(peaks, lodcolumn=="fruit_length" & chr==2)
    sw_peak <- subset(peaks, lodcolumn=="seed_weight" & chr==1 & pos > 20)

    fl_marker <- find_marker(arab$pmap, fl_peak$chr, fl_peak$pos)
    fl_pr <- pull_genoprobpos(pr, fl_marker)
    fit1_fl <- fit1(fl_pr, arab$pheno[,"fruit_length"])
    blup_fl <- fit1(fl_pr, arab$pheno[,"fruit_length"], blup=TRUE)

    sw_marker <- find_marker(arab$pmap, sw_peak$chr, sw_peak$pos)
    sw_pr <- pull_genoprobpos(pr, sw_marker)
    fit1_sw <- fit1(sw_pr, arab$pheno[,"seed_weight"])
    blup_sw <- fit1(sw_pr, arab$pheno[,"seed_weight"], blup=TRUE)
    # LL is quite different

    m_fl <- maxmarg(pr, arab$pmap, chr=fl_peak$chr, pos=fl_peak$pos, return_char=TRUE)
    m_sw <- maxmarg(pr, arab$pmap, chr=sw_peak$chr, pos=sw_peak$pos, return_char=TRUE)

    save(fl_peak, fit1_fl, blup_fl,
         sw_peak, fit1_sw, blup_sw,
         m_fl, m_sw, file=file)
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
