// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// assureBoundsCPP
NumericVector assureBoundsCPP(NumericVector ind, NumericVector g, NumericVector lower, NumericVector upper);
RcppExport SEXP _mogsa_assureBoundsCPP(SEXP indSEXP, SEXP gSEXP, SEXP lowerSEXP, SEXP upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ind(indSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type upper(upperSEXP);
    rcpp_result_gen = Rcpp::wrap(assureBoundsCPP(ind, g, lower, upper));
    return rcpp_result_gen;
END_RCPP
}
// computeVectorLengthCPP
double computeVectorLengthCPP(NumericVector vec);
RcppExport SEXP _mogsa_computeVectorLengthCPP(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(computeVectorLengthCPP(vec));
    return rcpp_result_gen;
END_RCPP
}
// normalizeVectorCPP
NumericVector normalizeVectorCPP(NumericVector vec, double prec);
RcppExport SEXP _mogsa_normalizeVectorCPP(SEXP vecSEXP, SEXP precSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< double >::type prec(precSEXP);
    rcpp_result_gen = Rcpp::wrap(normalizeVectorCPP(vec, prec));
    return rcpp_result_gen;
END_RCPP
}
// computeAngleCPP
double computeAngleCPP(NumericVector vec1, NumericVector vec2, double prec);
RcppExport SEXP _mogsa_computeAngleCPP(SEXP vec1SEXP, SEXP vec2SEXP, SEXP precSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec1(vec1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec2(vec2SEXP);
    Rcpp::traits::input_parameter< double >::type prec(precSEXP);
    rcpp_result_gen = Rcpp::wrap(computeAngleCPP(vec1, vec2, prec));
    return rcpp_result_gen;
END_RCPP
}
// findNextCellCPP
IntegerVector findNextCellCPP(NumericVector gradient);
RcppExport SEXP _mogsa_findNextCellCPP(SEXP gradientSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type gradient(gradientSEXP);
    rcpp_result_gen = Rcpp::wrap(findNextCellCPP(gradient));
    return rcpp_result_gen;
END_RCPP
}
// convertIndices2CellIDCPP
int convertIndices2CellIDCPP(IntegerVector indices, IntegerVector dims);
RcppExport SEXP _mogsa_convertIndices2CellIDCPP(SEXP indicesSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(convertIndices2CellIDCPP(indices, dims));
    return rcpp_result_gen;
END_RCPP
}
// convertCellID2IndicesCPP
IntegerVector convertCellID2IndicesCPP(int cellID, IntegerVector dims);
RcppExport SEXP _mogsa_convertCellID2IndicesCPP(SEXP cellIDSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type cellID(cellIDSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(convertCellID2IndicesCPP(cellID, dims));
    return rcpp_result_gen;
END_RCPP
}
// cumulateGradientsCPP
NumericVector cumulateGradientsCPP(NumericMatrix centers, NumericMatrix gradients, double precVectorLength, double precNorm);
RcppExport SEXP _mogsa_cumulateGradientsCPP(SEXP centersSEXP, SEXP gradientsSEXP, SEXP precVectorLengthSEXP, SEXP precNormSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type centers(centersSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gradients(gradientsSEXP);
    Rcpp::traits::input_parameter< double >::type precVectorLength(precVectorLengthSEXP);
    Rcpp::traits::input_parameter< double >::type precNorm(precNormSEXP);
    rcpp_result_gen = Rcpp::wrap(cumulateGradientsCPP(centers, gradients, precVectorLength, precNorm));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mogsa_assureBoundsCPP", (DL_FUNC) &_mogsa_assureBoundsCPP, 4},
    {"_mogsa_computeVectorLengthCPP", (DL_FUNC) &_mogsa_computeVectorLengthCPP, 1},
    {"_mogsa_normalizeVectorCPP", (DL_FUNC) &_mogsa_normalizeVectorCPP, 2},
    {"_mogsa_computeAngleCPP", (DL_FUNC) &_mogsa_computeAngleCPP, 3},
    {"_mogsa_findNextCellCPP", (DL_FUNC) &_mogsa_findNextCellCPP, 1},
    {"_mogsa_convertIndices2CellIDCPP", (DL_FUNC) &_mogsa_convertIndices2CellIDCPP, 2},
    {"_mogsa_convertCellID2IndicesCPP", (DL_FUNC) &_mogsa_convertCellID2IndicesCPP, 2},
    {"_mogsa_cumulateGradientsCPP", (DL_FUNC) &_mogsa_cumulateGradientsCPP, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_mogsa(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
