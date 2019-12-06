R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

STEM = msu2019

FIGS = Figs/rqtl_lines_code.pdf \
	   Figs/phefile.pdf \
	   Figs/geno_reconstruct.pdf \
	   Figs/scan_hk.pdf \
	   Figs/snp_asso.pdf \
	   Figs/coef_fl.pdf \
	   Figs/intercross.pdf \
	   Figs/lodcurve_insulin_with_effects.pdf \
	   Figs/congenic.pdf \
	   Figs/ail.pdf \
	   Figs/rilines.pdf \
	   Figs/ri8.pdf \
	   Figs/hs.pdf

docs/$(STEM).pdf: $(STEM).pdf
	cp $< $@

$(STEM).pdf: $(STEM).tex header.tex $(FIGS)
	xelatex $<

Figs/rqtl_lines_code.pdf: R/rqtl_lines_code.R R/colors.R Data/lines_code_by_version.csv
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

Figs/phefile.pdf: R/data_files_fig.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

Figs/geno_reconstruct.pdf: R/geno_reconstruct_fig.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

Figs/scan_hk.pdf: R/scans_fig.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

Figs/snp_asso.pdf: R/snp_asso_fig.R Figs/scan_hk.pdf
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

Figs/coef_fl.pdf: R/coef_fig.R Figs/scan_hk.pdf
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

Data/lines_code_by_version.csv: Perl/grab_lines_code.pl Data/versions.txt
	cd Perl;grab_lines_code.pl

Figs/magic19_scan.pdf: R/magic19_figs.R R/colors.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

Figs/intercross.pdf: R/intercross.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

Figs/lodcurve_insulin_with_effects.pdf: R/lodcurve_insulin.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

Figs/congenic.pdf: R/congenic_fig.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

Figs/ail.pdf: R/ail_fig.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

Figs/rilines.pdf: R/rilines_fig.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

Figs/ri8.pdf: R/ri8_fig.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

Figs/hs.pdf: R/hs_fig.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"
