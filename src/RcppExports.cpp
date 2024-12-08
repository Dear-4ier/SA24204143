// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// karatsubaMultiplyRcpp
std::string karatsubaMultiplyRcpp(const std::string& num1, const std::string& num2);
RcppExport SEXP _SA24204143_karatsubaMultiplyRcpp(SEXP num1SEXP, SEXP num2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type num1(num1SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type num2(num2SEXP);
    rcpp_result_gen = Rcpp::wrap(karatsubaMultiplyRcpp(num1, num2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA24204143_karatsubaMultiplyRcpp", (DL_FUNC) &_SA24204143_karatsubaMultiplyRcpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA24204143(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
