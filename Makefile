# dev makefile inspired from @Zilong-Li https://github.com/Zilong-Li/vcfppR/blob/main/Makefile 
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)
CPP_SRCS := $(wildcard src/*.cpp src/*.c src/*.h)
all: rcppCompile rd check clean

rcppCompile: $(CPP_SRCS)
	echo "Making RcppExport files..." && \
	echo "C++ source files: $(CPP_SRCS)" && \
	Rscript -e 'Rcpp::compileAttributes()'

rd: rcppCompile
	Rscript -e 'roxygen2::roxygenise(".")'

build: rd
	R CMD build .

install: build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check:
	Rscript -e 'devtools::check()'

check2: build
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

clean:
	$(RM) -r $(PKGNAME).Rcheck/

test: install
	Rscript -e 'devtools::test()'

