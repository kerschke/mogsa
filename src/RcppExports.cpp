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
// crossProductCPP
NumericVector crossProductCPP(NumericVector ab, NumericVector ac);
RcppExport SEXP _mogsa_crossProductCPP(SEXP abSEXP, SEXP acSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ab(abSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ac(acSEXP);
    rcpp_result_gen = Rcpp::wrap(crossProductCPP(ab, ac));
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
// gridBasedGradientCPP
NumericMatrix gridBasedGradientCPP(NumericVector fnVec, IntegerVector dims, NumericVector stepSizes, double precNorm, double precAngle);
RcppExport SEXP _mogsa_gridBasedGradientCPP(SEXP fnVecSEXP, SEXP dimsSEXP, SEXP stepSizesSEXP, SEXP precNormSEXP, SEXP precAngleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type fnVec(fnVecSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stepSizes(stepSizesSEXP);
    Rcpp::traits::input_parameter< double >::type precNorm(precNormSEXP);
    Rcpp::traits::input_parameter< double >::type precAngle(precAngleSEXP);
    rcpp_result_gen = Rcpp::wrap(gridBasedGradientCPP(fnVec, dims, stepSizes, precNorm, precAngle));
    return rcpp_result_gen;
END_RCPP
}
// cumulateGradientsCPP
NumericVector cumulateGradientsCPP(NumericMatrix centers, NumericMatrix gradients, double precVectorLength, double precNorm, bool fixDiagonals);
RcppExport SEXP _mogsa_cumulateGradientsCPP(SEXP centersSEXP, SEXP gradientsSEXP, SEXP precVectorLengthSEXP, SEXP precNormSEXP, SEXP fixDiagonalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type centers(centersSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gradients(gradientsSEXP);
    Rcpp::traits::input_parameter< double >::type precVectorLength(precVectorLengthSEXP);
    Rcpp::traits::input_parameter< double >::type precNorm(precNormSEXP);
    Rcpp::traits::input_parameter< bool >::type fixDiagonals(fixDiagonalsSEXP);
    rcpp_result_gen = Rcpp::wrap(cumulateGradientsCPP(centers, gradients, precVectorLength, precNorm, fixDiagonals));
    return rcpp_result_gen;
END_RCPP
}
// getBiObjGradientCPP
NumericVector getBiObjGradientCPP(NumericVector g1, NumericVector g2, double precNorm, double precAngle);
RcppExport SEXP _mogsa_getBiObjGradientCPP(SEXP g1SEXP, SEXP g2SEXP, SEXP precNormSEXP, SEXP precAngleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type g1(g1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g2(g2SEXP);
    Rcpp::traits::input_parameter< double >::type precNorm(precNormSEXP);
    Rcpp::traits::input_parameter< double >::type precAngle(precAngleSEXP);
    rcpp_result_gen = Rcpp::wrap(getBiObjGradientCPP(g1, g2, precNorm, precAngle));
    return rcpp_result_gen;
END_RCPP
}
// getTriObjGradientCPP
NumericVector getTriObjGradientCPP(NumericVector g1, NumericVector g2, NumericVector g3, double precNorm, double precAngle);
RcppExport SEXP _mogsa_getTriObjGradientCPP(SEXP g1SEXP, SEXP g2SEXP, SEXP g3SEXP, SEXP precNormSEXP, SEXP precAngleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type g1(g1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g2(g2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g3(g3SEXP);
    Rcpp::traits::input_parameter< double >::type precNorm(precNormSEXP);
    Rcpp::traits::input_parameter< double >::type precAngle(precAngleSEXP);
    rcpp_result_gen = Rcpp::wrap(getTriObjGradientCPP(g1, g2, g3, precNorm, precAngle));
    return rcpp_result_gen;
END_RCPP
}
// getBiObjGradientGridCPP
NumericMatrix getBiObjGradientGridCPP(NumericMatrix gradMat1, NumericMatrix gradMat2, double precNorm, double precAngle);
RcppExport SEXP _mogsa_getBiObjGradientGridCPP(SEXP gradMat1SEXP, SEXP gradMat2SEXP, SEXP precNormSEXP, SEXP precAngleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type gradMat1(gradMat1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gradMat2(gradMat2SEXP);
    Rcpp::traits::input_parameter< double >::type precNorm(precNormSEXP);
    Rcpp::traits::input_parameter< double >::type precAngle(precAngleSEXP);
    rcpp_result_gen = Rcpp::wrap(getBiObjGradientGridCPP(gradMat1, gradMat2, precNorm, precAngle));
    return rcpp_result_gen;
END_RCPP
}
// getTriObjGradientGridCPP
NumericMatrix getTriObjGradientGridCPP(NumericMatrix gradMat1, NumericMatrix gradMat2, NumericMatrix gradMat3, double precNorm, double precAngle);
RcppExport SEXP _mogsa_getTriObjGradientGridCPP(SEXP gradMat1SEXP, SEXP gradMat2SEXP, SEXP gradMat3SEXP, SEXP precNormSEXP, SEXP precAngleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type gradMat1(gradMat1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gradMat2(gradMat2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gradMat3(gradMat3SEXP);
    Rcpp::traits::input_parameter< double >::type precNorm(precNormSEXP);
    Rcpp::traits::input_parameter< double >::type precAngle(precAngleSEXP);
    rcpp_result_gen = Rcpp::wrap(getTriObjGradientGridCPP(gradMat1, gradMat2, gradMat3, precNorm, precAngle));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mogsa_assureBoundsCPP", (DL_FUNC) &_mogsa_assureBoundsCPP, 4},
    {"_mogsa_crossProductCPP", (DL_FUNC) &_mogsa_crossProductCPP, 2},
    {"_mogsa_computeVectorLengthCPP", (DL_FUNC) &_mogsa_computeVectorLengthCPP, 1},
    {"_mogsa_normalizeVectorCPP", (DL_FUNC) &_mogsa_normalizeVectorCPP, 2},
    {"_mogsa_computeAngleCPP", (DL_FUNC) &_mogsa_computeAngleCPP, 3},
    {"_mogsa_findNextCellCPP", (DL_FUNC) &_mogsa_findNextCellCPP, 1},
    {"_mogsa_convertIndices2CellIDCPP", (DL_FUNC) &_mogsa_convertIndices2CellIDCPP, 2},
    {"_mogsa_convertCellID2IndicesCPP", (DL_FUNC) &_mogsa_convertCellID2IndicesCPP, 2},
    {"_mogsa_gridBasedGradientCPP", (DL_FUNC) &_mogsa_gridBasedGradientCPP, 5},
    {"_mogsa_cumulateGradientsCPP", (DL_FUNC) &_mogsa_cumulateGradientsCPP, 5},
    {"_mogsa_getBiObjGradientCPP", (DL_FUNC) &_mogsa_getBiObjGradientCPP, 4},
    {"_mogsa_getTriObjGradientCPP", (DL_FUNC) &_mogsa_getTriObjGradientCPP, 5},
    {"_mogsa_getBiObjGradientGridCPP", (DL_FUNC) &_mogsa_getBiObjGradientGridCPP, 4},
    {"_mogsa_getTriObjGradientGridCPP", (DL_FUNC) &_mogsa_getTriObjGradientGridCPP, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_mogsa(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
