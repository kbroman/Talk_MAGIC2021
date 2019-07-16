library(qtl2)

url <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/ArabMAGIC/"
arab <- read_cross2(paste0(url, "arabmagic_tair9.zip"))
arab <- recode_snps(arab)

arab <- pull_markers(arab, unlist(lapply(reduce_markers(arab$pmap, min_distance=0.02), names)))

pr <- calc_genoprob(arab["MAGIC.244",3], arab$pmap, error_prob=0.002, cores=24)

m <- maxmarg(pr, minprob=0.5)

f_to_show <- 1:19


maxy_founder <- 5.3
y_do <- 6.1
y_doinf <- 6.5
pt_cex <- 1.8
pt_cex_small <- 0.8


geno_reconstruct <- function(version=1)
{
    geno_colors <- c("#CBE4F0", "gray", "#0071BC")
    hap_colors <- c("#49A56E", "#F15A24")
    fg_color <- "white"

    par(mar=c(3.1, 5.6, 0.6, 2.6), bty="n")
    par(fg=fg_color, col.axis=fg_color, col.lab=fg_color)

    plot(0,0,type="n", xlab="Chr 3 position (Mbp)", ylab="",
         yaxt="n", xlim=c(2, 6), xaxs="i", ylim=c(6.8, 0),
         xaxt="n", mgp=c(1.8, 0, 0))
    axis(side=1, tick=FALSE, mgp=c(0, 0.4, 0))
    abline(v=seq(45, 46, by=0.2), col="gray80")

    y <- seq(0, maxy_founder, length=length(f_to_show))
    abline(h=y)

    fg <- arab$founder_geno[["3"]][f_to_show,]
    magic_ind <- "MAGIC.244"
    g <- arab$geno[["3"]][magic_ind,]
    pmap <- arab$pmap[["3"]]

    for(i in seq_along(f_to_show)) {
        typed <- (fg[i,]!=0)
        yval <- rep(y[i], ncol(fg))[typed]
        xval <- pmap[typed]


        points(xval, yval, pch=21, col=fg_color,
               bg=geno_colors[fg[i, typed]],
               cex=pt_cex)
    }

    axis(side=2, at=y, rownames(fg), tick=FALSE, mgp=c(0, 0.3, 0), las=1)

    abline(h=y_do)

    typed <- (g != 0)
    yval <- rep(y_do, ncol(fg))[typed]
    xval <- pmap[typed]
    points(xval, yval, pch=21, col=fg_color,
           bg=geno_colors[g[typed]],
           cex=pt_cex)
    axis(side=2, at=y_do, magic_ind, las=1, mgp=c(0, 0.3, 0), tick=FALSE)

    if(version < 2) return()

    u <- par("usr")
    yd <- 0.2

    m <- maxmarg(pr, minprob=0.8)[[1]][1,]

    index_left <- max(which(!is.na(m) & pmap > 3 & pmap < 4))
    index_right <- min(which(!is.na(m) & pmap > 4 & pmap < 5))
    par_left <- rownames(arab$founder_geno[[1]])[m[index_left]]
    par_right <- rownames(arab$founder_geno[[1]])[m[index_right]]
    rec_left <- pmap[index_left]
    rec_right <- pmap[index_right]

    rect(u[1], y_doinf+yd*0.1, rec_left, y_doinf+yd*2, border=NA, col=hap_colors[1])
    rect(rec_right, y_doinf+yd*0.1, u[2], y_doinf+yd*2, border=NA, col=hap_colors[2])

    axis(side=2, at=y_doinf+yd*1.05, par_left, tick=FALSE, mgp=c(0, 0.3, 0), las=1)
    axis(side=4, at=y_doinf+yd*1.05, par_right, tick=FALSE, mgp=c(0, 0.3, 0), las=1)

}

pdf("../Figs/geno_reconstruct.pdf", height=6.5, width=11, pointsize=14)
geno_reconstruct()
dev.off()

pdf("../Figs/geno_reconstruct_B.pdf", height=6.5, width=11, pointsize=14)
geno_reconstruct(2)
dev.off()
