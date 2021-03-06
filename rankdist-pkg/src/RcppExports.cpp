// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// KendallNeighbour
NumericMatrix KendallNeighbour(NumericVector rank);
RcppExport SEXP rankdist_KendallNeighbour(SEXP rankSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rank(rankSEXP);
    rcpp_result_gen = Rcpp::wrap(KendallNeighbour(rank));
    return rcpp_result_gen;
END_RCPP
}
// CayleyNeighbour
NumericMatrix CayleyNeighbour(NumericVector rank);
RcppExport SEXP rankdist_CayleyNeighbour(SEXP rankSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rank(rankSEXP);
    rcpp_result_gen = Rcpp::wrap(CayleyNeighbour(rank));
    return rcpp_result_gen;
END_RCPP
}
// LogC
double LogC(NumericVector fai);
RcppExport SEXP rankdist_LogC(SEXP faiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type fai(faiSEXP);
    rcpp_result_gen = Rcpp::wrap(LogC(fai));
    return rcpp_result_gen;
END_RCPP
}
// CWeightGivenPi
NumericVector CWeightGivenPi(NumericVector r1, NumericVector r2);
RcppExport SEXP rankdist_CWeightGivenPi(SEXP r1SEXP, SEXP r2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type r1(r1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r2(r2SEXP);
    rcpp_result_gen = Rcpp::wrap(CWeightGivenPi(r1, r2));
    return rcpp_result_gen;
END_RCPP
}
// FindV
NumericMatrix FindV(NumericMatrix obs, NumericVector pi0);
RcppExport SEXP rankdist_FindV(SEXP obsSEXP, SEXP pi0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi0(pi0SEXP);
    rcpp_result_gen = Rcpp::wrap(FindV(obs, pi0));
    return rcpp_result_gen;
END_RCPP
}
// LogC_Component
double LogC_Component(NumericVector fai);
RcppExport SEXP rankdist_LogC_Component(SEXP faiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type fai(faiSEXP);
    rcpp_result_gen = Rcpp::wrap(LogC_Component(fai));
    return rcpp_result_gen;
END_RCPP
}
// cycle_decomp
int cycle_decomp(NumericVector comp);
RcppExport SEXP rankdist_cycle_decomp(SEXP compSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type comp(compSEXP);
    rcpp_result_gen = Rcpp::wrap(cycle_decomp(comp));
    return rcpp_result_gen;
END_RCPP
}
// FindCayley
NumericVector FindCayley(NumericMatrix obs, NumericVector pi0);
RcppExport SEXP rankdist_FindCayley(SEXP obsSEXP, SEXP pi0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi0(pi0SEXP);
    rcpp_result_gen = Rcpp::wrap(FindCayley(obs, pi0));
    return rcpp_result_gen;
END_RCPP
}
// Wtau
NumericMatrix Wtau(NumericMatrix obs, NumericVector pi0);
RcppExport SEXP rankdist_Wtau(SEXP obsSEXP, SEXP pi0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi0(pi0SEXP);
    rcpp_result_gen = Rcpp::wrap(Wtau(obs, pi0));
    return rcpp_result_gen;
END_RCPP
}
// AllPerms
NumericMatrix AllPerms(int nobj);
RcppExport SEXP rankdist_AllPerms(SEXP nobjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nobj(nobjSEXP);
    rcpp_result_gen = Rcpp::wrap(AllPerms(nobj));
    return rcpp_result_gen;
END_RCPP
}
