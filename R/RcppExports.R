# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

assureBoundsCPP <- function(ind, g, lower, upper) {
    .Call(`_mogsa_assureBoundsCPP`, ind, g, lower, upper)
}

computeVectorLengthCPP <- function(vec) {
    .Call(`_mogsa_computeVectorLengthCPP`, vec)
}

normalizeVectorCPP <- function(vec, prec) {
    .Call(`_mogsa_normalizeVectorCPP`, vec, prec)
}

computeAngleCPP <- function(vec1, vec2, prec) {
    .Call(`_mogsa_computeAngleCPP`, vec1, vec2, prec)
}

findNextCellCPP <- function(gradient) {
    .Call(`_mogsa_findNextCellCPP`, gradient)
}

convertIndices2CellIDCPP <- function(indices, dims) {
    .Call(`_mogsa_convertIndices2CellIDCPP`, indices, dims)
}

convertCellID2IndicesCPP <- function(cellID, dims) {
    .Call(`_mogsa_convertCellID2IndicesCPP`, cellID, dims)
}

cumulateGradientsCPP <- function(centers, gradients, precVectorLength, precNorm) {
    .Call(`_mogsa_cumulateGradientsCPP`, centers, gradients, precVectorLength, precNorm)
}

