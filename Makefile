R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

STEM = magic2019

FIGS = Figs/rqtl_lines_code.pdf \
	   Figs/phefile.pdf \
	   Figs/geno_reconstruct.pdf \
	   Figs/scan_hk.pdf \
	   Figs/snp_asso.pdf

docs/$(STEM).pdf: $(STEM).pdf
	cp $< $@

$(STEM).pdf: $(STEM).tex header.tex $(FIGS)
	xelatex $<

web: $(STEM).pdf
	scp $(STEM).pdf adhara.biostat.wisc.edu:Website/presentations/$(STEM).pdf

Figs/rqtl_lines_code.pdf: R/rqtl_lines_code.R R/colors.R Data/lines_code_by_version.csv
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

Figs/phefile.pdf: R/data_files_fig.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

Figs/geno_reconstruct.pdf: R/geno_reconstruct_fig.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

Figs/scan_hk.pdf: R/scans_fig.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

Figs/snp_asso.pdf: R/snp_asso_fig.R
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

Data/lines_code_by_version.csv: Perl/grab_lines_code.pl Data/versions.txt
	cd Perl;grab_lines_code.pl

Figs/magic19_scan.pdf: R/magic19_figs.R R/colors.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"
