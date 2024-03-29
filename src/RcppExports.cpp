// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// BigRowSums
arma::vec BigRowSums(SEXP pBigMat);
RcppExport SEXP _GCMSpro_BigRowSums(SEXP pBigMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    rcpp_result_gen = Rcpp::wrap(BigRowSums(pBigMat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GCMSpro_BigRowSums", (DL_FUNC) &_GCMSpro_BigRowSums, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_GCMSpro(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
