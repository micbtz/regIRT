// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// nominallikRcppA
arma::mat nominallikRcppA(arma::vec par, arma::mat data, arma::vec nodes, arma::vec weights, arma::mat numpatt, arma::uvec itemsselect, double lambda);
RcppExport SEXP _regIRT_nominallikRcppA(SEXP parSEXP, SEXP dataSEXP, SEXP nodesSEXP, SEXP weightsSEXP, SEXP numpattSEXP, SEXP itemsselectSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type numpatt(numpattSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type itemsselect(itemsselectSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(nominallikRcppA(par, data, nodes, weights, numpatt, itemsselect, lambda));
    return rcpp_result_gen;
END_RCPP
}
// gradnominallikRcppA
arma::mat gradnominallikRcppA(arma::vec par, arma::mat data, arma::vec nodes, arma::vec weights, arma::mat numpatt, arma::uvec itemsselect, double lambda);
RcppExport SEXP _regIRT_gradnominallikRcppA(SEXP parSEXP, SEXP dataSEXP, SEXP nodesSEXP, SEXP weightsSEXP, SEXP numpattSEXP, SEXP itemsselectSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type numpatt(numpattSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type itemsselect(itemsselectSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(gradnominallikRcppA(par, data, nodes, weights, numpatt, itemsselect, lambda));
    return rcpp_result_gen;
END_RCPP
}
// nominallikaugRcppA
arma::mat nominallikaugRcppA(arma::vec par, arma::mat data, arma::vec nodes, arma::vec weights, List gamma, List v, double c, arma::mat numpatt, arma::uvec itemsselect, double lambda);
RcppExport SEXP _regIRT_nominallikaugRcppA(SEXP parSEXP, SEXP dataSEXP, SEXP nodesSEXP, SEXP weightsSEXP, SEXP gammaSEXP, SEXP vSEXP, SEXP cSEXP, SEXP numpattSEXP, SEXP itemsselectSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< List >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type numpatt(numpattSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type itemsselect(itemsselectSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(nominallikaugRcppA(par, data, nodes, weights, gamma, v, c, numpatt, itemsselect, lambda));
    return rcpp_result_gen;
END_RCPP
}
// gradnominallikaugRcppA
arma::mat gradnominallikaugRcppA(arma::vec par, arma::mat data, arma::vec nodes, arma::vec weights, List gamma, List v, double c, arma::mat numpatt, arma::uvec itemsselect, double lambda);
RcppExport SEXP _regIRT_gradnominallikaugRcppA(SEXP parSEXP, SEXP dataSEXP, SEXP nodesSEXP, SEXP weightsSEXP, SEXP gammaSEXP, SEXP vSEXP, SEXP cSEXP, SEXP numpattSEXP, SEXP itemsselectSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< List >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< List >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type numpatt(numpattSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type itemsselect(itemsselectSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(gradnominallikaugRcppA(par, data, nodes, weights, gamma, v, c, numpatt, itemsselect, lambda));
    return rcpp_result_gen;
END_RCPP
}
// nominallik_fpenRcppA
arma::mat nominallik_fpenRcppA(arma::vec par, arma::mat data, arma::vec nodes, arma::vec weights, arma::mat numpatt, arma::uvec itemsselect, List w, double lambda, double eps);
RcppExport SEXP _regIRT_nominallik_fpenRcppA(SEXP parSEXP, SEXP dataSEXP, SEXP nodesSEXP, SEXP weightsSEXP, SEXP numpattSEXP, SEXP itemsselectSEXP, SEXP wSEXP, SEXP lambdaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type numpatt(numpattSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type itemsselect(itemsselectSEXP);
    Rcpp::traits::input_parameter< List >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(nominallik_fpenRcppA(par, data, nodes, weights, numpatt, itemsselect, w, lambda, eps));
    return rcpp_result_gen;
END_RCPP
}
// gradnominallik_fpenRcppA
arma::mat gradnominallik_fpenRcppA(arma::vec par, arma::mat data, arma::vec nodes, arma::vec weights, arma::mat numpatt, arma::uvec itemsselect, List w, double lambda, double eps);
RcppExport SEXP _regIRT_gradnominallik_fpenRcppA(SEXP parSEXP, SEXP dataSEXP, SEXP nodesSEXP, SEXP weightsSEXP, SEXP numpattSEXP, SEXP itemsselectSEXP, SEXP wSEXP, SEXP lambdaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type numpatt(numpattSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type itemsselect(itemsselectSEXP);
    Rcpp::traits::input_parameter< List >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(gradnominallik_fpenRcppA(par, data, nodes, weights, numpatt, itemsselect, w, lambda, eps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_regIRT_nominallikRcppA", (DL_FUNC) &_regIRT_nominallikRcppA, 7},
    {"_regIRT_gradnominallikRcppA", (DL_FUNC) &_regIRT_gradnominallikRcppA, 7},
    {"_regIRT_nominallikaugRcppA", (DL_FUNC) &_regIRT_nominallikaugRcppA, 10},
    {"_regIRT_gradnominallikaugRcppA", (DL_FUNC) &_regIRT_gradnominallikaugRcppA, 10},
    {"_regIRT_nominallik_fpenRcppA", (DL_FUNC) &_regIRT_nominallik_fpenRcppA, 9},
    {"_regIRT_gradnominallik_fpenRcppA", (DL_FUNC) &_regIRT_gradnominallik_fpenRcppA, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_regIRT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
