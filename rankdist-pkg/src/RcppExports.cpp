// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// KendallNeighbour
NumericMatrix KendallNeighbour(NumericVector rank);
RcppExport SEXP rankdist_KendallNeighbour(SEXP rankSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type rank(rankSEXP );
        NumericMatrix __result = KendallNeighbour(rank);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// CayleyNeighbour
NumericMatrix CayleyNeighbour(NumericVector rank);
RcppExport SEXP rankdist_CayleyNeighbour(SEXP rankSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type rank(rankSEXP );
        NumericMatrix __result = CayleyNeighbour(rank);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// LogC
double LogC(NumericVector fai);
RcppExport SEXP rankdist_LogC(SEXP faiSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type fai(faiSEXP );
        double __result = LogC(fai);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// CWeightGivenPi
NumericVector CWeightGivenPi(NumericVector r1, NumericVector r2);
RcppExport SEXP rankdist_CWeightGivenPi(SEXP r1SEXP, SEXP r2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type r1(r1SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type r2(r2SEXP );
        NumericVector __result = CWeightGivenPi(r1, r2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// FindV
NumericMatrix FindV(NumericMatrix obs, NumericVector pi0);
RcppExport SEXP rankdist_FindV(SEXP obsSEXP, SEXP pi0SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type obs(obsSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type pi0(pi0SEXP );
        NumericMatrix __result = FindV(obs, pi0);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// LogC_Component
double LogC_Component(NumericVector fai);
RcppExport SEXP rankdist_LogC_Component(SEXP faiSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type fai(faiSEXP );
        double __result = LogC_Component(fai);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cycle_decomp
int cycle_decomp(NumericVector comp);
RcppExport SEXP rankdist_cycle_decomp(SEXP compSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type comp(compSEXP );
        int __result = cycle_decomp(comp);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// FindCayley
NumericVector FindCayley(NumericMatrix obs, NumericVector pi0);
RcppExport SEXP rankdist_FindCayley(SEXP obsSEXP, SEXP pi0SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type obs(obsSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type pi0(pi0SEXP );
        NumericVector __result = FindCayley(obs, pi0);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Wtau
NumericMatrix Wtau(NumericMatrix obs, NumericVector pi0);
RcppExport SEXP rankdist_Wtau(SEXP obsSEXP, SEXP pi0SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type obs(obsSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type pi0(pi0SEXP );
        NumericMatrix __result = Wtau(obs, pi0);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// AllPerms
NumericMatrix AllPerms(int nobj);
RcppExport SEXP rankdist_AllPerms(SEXP nobjSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type nobj(nobjSEXP );
        NumericMatrix __result = AllPerms(nobj);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}