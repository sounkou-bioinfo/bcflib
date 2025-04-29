#include <Rcpp.h>
extern "C" {
#include "bcftools.h"
}
using namespace Rcpp;

//' Get BCFtools version
//'
//' Returns the version string of the linked BCFtools library.
//' @return Character string of BCFtools version.
//' @export
// [[Rcpp::export]]
std::string getBcftoolsVersion() {
    return std::string(bcftools_version());
}