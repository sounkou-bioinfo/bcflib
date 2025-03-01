# dev makefile inspired from @Zilong-Li https://github.com/Zilong-Li/vcfppR/blob/main/Makefile 
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: rcppCompile rd check clean


rcppCompile:
  Rscript -e 'Rcpp::compileAttributes()'
rd: rcppCompile
	Rscript -e 'roxygen2::roxygenise(".")'

#readme: rd
#	Rscript -e 'rmarkdown::render("README.Rmd")'
#	Rscript -e 'pkgdown::build_site(examples=FALSE)'
#	Rscript -e 'pkgdown::build_articles()'

#readme2: rd
#	Rscript -e 'rmarkdown::render("README.Rmd", "html_document")'

build: rd
	cd ..;\
	R CMD build $(PKGSRC)

install: build
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check:
	cd ..;\
	Rscript -e 'devtools::check()'

check2: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

