SHELL   := /bin/bash
PACKAGE := $(shell perl -aF: -ne 'print, exit if s/^Package:\s+//' DESCRIPTION)
VERSION := $(shell perl -aF: -ne 'print, exit if s/^Version:\s+//' DESCRIPTION)
BUILD   := $(PACKAGE)_$(VERSION).tar.gz

.PHONY: doc build install test vignette check-bioc $(BUILD)

check: $(BUILD)
	R CMD check --as-cran $<

check-rhub: $(BUILD)
	Rscript -e 'rhub::check("$(BUILD)", platform = c("ubuntu-gcc-devel", "windows-x86_64-devel", "solaris-x86-patched", "solaris-x86-patched-ods", "macos-m1-bigsur-release"))'

check-bioc: $(BUILD)
	Rscript -e "BiocCheck::BiocCheck('$(BUILD)')"

check-solaris: $(BUILD)
	Rscript -e 'rhub::check("$(BUILD)", platform = c("solaris-x86-patched", "solaris-x86-patched-ods"))'

check-m1: $(BUILD)
	Rscript -e 'rhub::check("$(BUILD)", platform = c("macos-m1-bigsur-release"))'

reconf:
	autoreconf -fi

build:
	find . -type d -exec chmod 755 {} \;
	find . -type f -exec chmod 644 {} \;
	chmod 755 cleanup
	chmod 755 configure
	./configure
	./cleanup
	Rscript -e "library(Rcpp); compileAttributes('.');"
	Rscript -e "devtools::load_all(); roxygen2::roxygenise('.');"
	find . -iname "*.a" -exec rm {} \;
	find . -iname "*.o" -exec rm {} \;
	find . -iname "*.so" -exec rm {} \;
	R CMD build .

install:
	find . -type d -exec chmod 755 {} \;
	find . -type f -exec chmod 644 {} \;
	chmod 755 cleanup
	chmod 755 configure
	./configure
	./cleanup
	Rscript -e "library(Rcpp); compileAttributes('.');"
	Rscript -e "devtools::load_all(); roxygen2::roxygenise('.');"
	find . -iname "*.a" -exec rm {} \;
	find . -iname "*.o" -exec rm {} \;
	find . -iname "*.so" -exec rm {} \;
	R CMD build .
	R CMD INSTALL $(BUILD)

install-no-vignette:
	find . -type d -exec chmod 755 {} \;
	find . -type f -exec chmod 644 {} \;
	chmod 755 cleanup
	chmod 755 configure
	./configure
	./cleanup
	Rscript -e "library(Rcpp); compileAttributes('.');"
	Rscript -e "devtools::load_all(); roxygen2::roxygenise('.');"
	find . -iname "*.a" -exec rm {} \;
	find . -iname "*.o" -exec rm {} \;
	find . -iname "*.so" -exec rm {} \;
	R CMD build . --no-build-vignettes
	R CMD INSTALL $(BUILD)

vignette:
	Rscript -e "rmarkdown::render('vignettes/fraq_getting_started.Rmd', output_format = BiocStyle::html_document())"
	Rscript -e "rmarkdown::render('vignettes/fraq_getting_started.Rmd', output_format=rmarkdown::github_document(html_preview=FALSE))"
	mv vignettes/fraq_getting_started.md README.md
	perl -pi -e 's{plots/fig-qc-quicklook-1\.png}{vignettes/plots/fig-qc-quicklook-1.png}g' README.md
	perl -pi -e 's{plots/fraq_flow_graph_mermaid\.png}{vignettes/plots/fraq_flow_graph_mermaid.png}g' README.md
	perl -0pi -e 's/\n## Session information.*//s' README.md

test: $(BUILD)
	FRAQ_EXTENDED_TESTS=1 Rscript tests/fraq_prebuilt_tests.R
	FRAQ_EXTENDED_TESTS=1 Rscript tests/fraq_format_tests.R
	FRAQ_EXTENDED_TESTS=1 Rscript tests/fraq_fifo_tests.R
	FRAQ_EXTENDED_TESTS=1 Rscript tests/fraq_concat_tests.R
	FRAQ_EXTENDED_TESTS=1 Rscript tests/fraq_utils_tests.R
