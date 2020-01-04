#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <math.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector assureBoundsCPP(NumericVector ind, NumericVector g, NumericVector lower, NumericVector upper) {
  // ensures that a gradient step for an individual ind into direction g
  // ends in a position, where it remains within the lower and upper bounds

  int d = ind.size();                  // problem dimension 
  int i;                               // counter for the for-loops
  NumericVector child = ind + g;       // possible child (= location resulting from the gradient step)
  double stepMax = std::numeric_limits<double>::infinity(); // maximum step length
  double step;                         // step length to ensure to stay in bounds

  // check for lower bounds
  int stepPos = d;
  LogicalVector flagLow = (child < lower);
  if (is_true(any(flagLow))) {
    // adjust to stay within ALL lower bounds
    for (i = 0; i < d; i++) {
      bool fl = flagLow[i];
      if (fl) {
        step = (ind[i] - lower[i]) / abs(g[i]);
        if (step < stepMax) {
          stepMax = step;
          stepPos = i;
        }
      }
    }
    if (stepPos < d) {
      child = ind + g / abs(g[stepPos]) * (ind[stepPos] - lower[stepPos]);
    }
  }

  // check for upper bounds
  stepPos = d;
  stepMax = std::numeric_limits<double>::infinity(); 
  LogicalVector flagUpp = (child > upper);
  if (is_true(any(flagUpp))) {
    // adjust to stay within ALL upper bounds
    for (i = 0; i < d; i++) {
      bool fu = flagUpp[i];
      if (fu) {
        step = (upper[i] - ind[i]) / abs(g[i]);
        if (step < stepMax) {
          stepMax = step;
          stepPos = i;
        }
      }
    }
    if (stepPos < d) {
      child = ind + g / abs(g[stepPos]) * (upper[stepPos] - ind[stepPos]);
    }
  }

  return child ;
}

// [[Rcpp::export]]
NumericVector crossProductCPP(NumericVector ab, NumericVector ac) {
  double abci = ab(1) * ac(2) - ac(1) * ab(2);
  double abcj = ac(0) * ab(2) - ab(0) * ac(2);
  double abck = ab(0) * ac(1) - ac(0) * ab(1);
  return NumericVector::create(abci,abcj,abck);
}

// [[Rcpp::export]]
double computeVectorLengthCPP(NumericVector vec) {
  // computes length of the vector
  return sqrt(sum(vec * vec)) ;
}

// [[Rcpp::export]]
NumericVector normalizeVectorCPP(NumericVector vec, double prec) {
  // normalizes vector vec by its length
  NumericVector vec2 = ifelse(abs(vec) < prec, 0, vec);
  if (is_false(all(vec2 == 0))) {
    vec2 = vec2 / computeVectorLengthCPP(vec2);
  }
  return vec2 ;
}

// [[Rcpp::export]]
double computeAngleCPP(NumericVector vec1, NumericVector vec2, double prec) {
  // computes angle spanned by vectors vec1 and vec2
  // prec = max absolute difference (per element) between the vectors

  double angle = 0;
  if (max(abs(vec1 + vec2)) < prec) {
    // if vec1 and vec2 show in opposite directions
    angle = 180;
  } else {
    if (max(abs(vec1 - vec2)) >= prec) {
      // if vec1 and vec2 are (almost) identical
      double enumerator = sum(vec1 * vec2);
      double denominator = computeVectorLengthCPP(vec1) * computeVectorLengthCPP(vec2);
      if (enumerator < denominator) {
        // the enumerator can only be larger than the denominator, in case of numerical imprecision
        double rad = acos(enumerator / denominator);
        double pi = atan(1) * 4;
        angle = 180 * rad / pi;
      }
    } 
  }
  return angle ;
}


// // [[Rcpp::export]]
// IntegerVector findNextCellCPP(double angle) {
//   // which is the next cell one has to move to (based on the angle)?
//   
//   IntegerVector direction(2, 3);
//   // ensure to only consider angles between 0 and 360 degrees
//   while (angle < 360) {
//     angle += 360.0;
//   };
//   angle = std::fmod(angle, 360.0);
//   
//   if ((angle > 22.5) & (angle <= 67.5)) {
//     direction[0] = 1; // step in x-direction
//     direction[1] = 1; // step in y-direction
//   } else if ((angle > 67.5) & (angle <= 112.5)) {
//     direction[0] = 1;
//     direction[1] = 0;
//   } else if ((angle > 112.5) & (angle <= 157.5)) {
//     direction[0] = 1;
//     direction[1] = -1;
//   } else if ((angle > 157.5) & (angle <= 202.5)) {
//     direction[0] = 0;
//     direction[1] = -1;
//   } else if ((angle > 202.5) & (angle <= 247.5)) {
//     direction[0] = -1;
//     direction[1] = -1;
//   } else if ((angle > 247.5) & (angle <= 292.5)) {
//     direction[0] = -1;
//     direction[1] = 0;
//   } else if ((angle > 292.5) & (angle <= 337.5)) {
//     direction[0] = -1;
//     direction[1] = 1;
//   } else {
//     direction[0] = 0;
//     direction[1] = 1;
//   }
//   
//   return direction ;
// }


// [[Rcpp::export]]
IntegerVector findNextCellCPP(NumericVector gradient) {
// IntegerVector findNextCellHighDimensionCPP(NumericVector gradient) {
  // which is the next cell one has to move to (based on the high-dimensional gradient)?

  int d = gradient.size();
  IntegerVector direction(d, 0.0);

  // we need the largest absolute value to stretch the gradient to the boundary
  double maxval = max(abs(gradient));

  // stretch gradient vector such that its longest component touches the boundary
  // of [-3, 3]^d; for each component with a (stretched) value >= 1 go in positive
  // direction, for all components with a value <= -1 go in the negative direction
  for (int j = 0; j < d; j++) {
    double vec = 1.5 * gradient[j] / maxval;
    if (vec >= 1.0) {
      direction[j] = 1; // step in positive direction of x[j]
    } else if (vec <= -1.0) {
      direction[j] = -1; // step in negative direction of x[j]
    }
  }
  return direction;
}


// // convert rowIndex (1, ..., nRow) and columnIndex (1, ..., nColumns) to cell ID (1, ..., nRows * nColumns)
// // [[Rcpp::export]]
// int convertIndices2CellIDCPP(int rowIndex, int columnIndex, int nRows = 100, int nColumns = 100) {
//   
//   int cellID = -1;
//   if ((columnIndex > nColumns) | (columnIndex < 1) | (rowIndex > nRows) | (rowIndex < 1)) {
//     // if either the row- or columnIndex are located outside the boundaries, return -1 as ID
//     cellID = -1;
//   } else {
//     cellID = (rowIndex - 1) * nColumns + columnIndex;
//   }
//   return cellID;
// }


// convert index per dimension [(1, ..., rows), (1, ..., columns), ...] to cell ID (1, ..., prod(dims))
// [[Rcpp::export]]
int convertIndices2CellIDCPP(IntegerVector indices, IntegerVector dims) {
// int convertHighDimensionalIndices2CellIDCPP(IntegerVector indices, IntegerVector dims) {

  int cellID = -1;
  if (is_true(any(indices > dims)) | is_true(any(indices < 1))) {
    // if either the row- or columnIndex are located outside the boundaries, return -1 as ID
    cellID = -1;
  } else {
    int d = dims.size();
    IntegerVector cumcells(d, 1);
    cumcells[0] = 1;
    for (int j = 1; j < d; j++) {
      cumcells[j] = cumcells[j - 1] * dims[j - 1];
    }
    cellID = 1;
    for (int j = 0; j < d; j++) {
      cellID += (indices[j] - 1) * cumcells[j];
    }
  }
  return cellID;
}


// // convert cellID (1, ..., nRows * nColumns) to rowIndex (1, ..., nRow) and columnIndex (1, ..., nColumns)
// // [[Rcpp::export]]
// IntegerVector convertCellID2IndicesCPP(int cellID, int nRows = 100, int nColumns = 100) {
// 
//   IntegerVector indexVector(2, -1);
//   if ((cellID > (nRows * nColumns)) | (cellID < 1)) {
//     // if the cell ID has an irregular value, return -1 for the row and column indices
//     indexVector = rep(-1, 2);
//   } else {
//     int x = cellID - 1;
//     int y = floor(x / nColumns);
//     indexVector[0] = y + 1;              // rowIndex
//     indexVector[1] = (x % nColumns) + 1; // columnIndex
//   }
//   return indexVector;
// }


// convert cellID (1, ..., prod(dims)) to index per dimension [(1, ..., rows), (1, ..., columns), ...]
// [[Rcpp::export]]
IntegerVector convertCellID2IndicesCPP(int cellID, IntegerVector dims) {
// IntegerVector convertHighDimensionalCellID2IndicesCPP(int cellID, IntegerVector dims) {

  int d = dims.size();
  IntegerVector indexVector(d, -1);

  IntegerVector cumcells(d, 1);
  cumcells[0] = dims[0];
  for (int j = 1; j < d; j++) {
    cumcells[j] = cumcells[j - 1] * dims[j];
  }

  if ((cellID <= cumcells[d - 1]) & (cellID >= 1)) {
    indexVector[0] = ((cellID - 1) % cumcells[0]) + 1;
    for (int j = 1; j < d; j++) {
      int tmp = floor((cellID - 1) / cumcells[j - 1]);
      indexVector[j] = (tmp % dims[j]) + 1;
    }
  }
  return indexVector;
}


// // compute the cumulated gradient matrix based on a grid of points and their corresponding gradients;
// // precVectorLength defines the threshold of the multi-objective gradient length, for which the
// // corresponding point is considered to be locally efficient (suggest: 0.001 or smaller);
// // precNorm is the threshold used when normalizing a vector (within the angle-computation)
// // [[Rcpp::export]]
// NumericVector cumulateGradientsCPP(NumericMatrix centers, NumericMatrix gradients, double precVectorLength, double precNorm) {
// 
//   // FIXME: so far, the code only supports 2D-problems
//   int d = centers.ncol();                   // dimensionality of search space
//   int n = centers.nrow();                   // number of grid points
// 
//   // number of grid points per (search space) dimension
//   NumericVector ctr1 = centers(_,0);
//   ctr1 = unique(ctr1);
//   int nColumns = ctr1.size();               // #columns = number of grid points in x1-dimension
//   NumericVector ctr2 = centers(_,1);
//   ctr2 = unique(ctr2);
//   int nRows = ctr2.size();                  // #rows = number of grid points in x2-dimension
// 
//   NumericVector gradientLengths(n);         // length vector for the multi-objective gradients
//   IntegerVector cellPointer(n, -999);       // vector, which indicates per cell the successor cell
//   NumericVector gradFieldVector(n, -999.0); // the "final" result vector
//   LogicalVector visited = rep(false, n);    // which cells have already been "visited" / "processed"
// 
//   // helper variables
//   double angle;
//   double vectorLength = -1.0;
//   NumericVector baseVector = NumericVector::create( 1.0, 0.0 );
//   NumericVector currentGradients;
//   int visitCounter = 0;                     // counter, enabling early stop of algorithm once all cells are processed
//   IntegerVector currentCell(d);             // row and column index of the current cell
//   IntegerVector nextCell(d);                // row and column index of the next (= successor) cell
//   int nextCellID;
// 
//   for (int i = 0; i < n; i++) {
//     // iterate over all cells, extract their gradients and store the gradient lengths
//     currentGradients = gradients(i,_);
//     vectorLength = computeVectorLengthCPP(currentGradients);
//     gradientLengths[i] = vectorLength;
//     if (vectorLength < precVectorLength) {
//       // if the multi-objective gradient is "short enough", the current cell is considered
//       // to be locally efficient (--> move on to the next cell)
//       visited[i] = true;
//       visitCounter++;
//       gradFieldVector[i] = vectorLength;
//     } else {
//       // compute the angle of the gradient and identify the succeeding cell
//       angle = computeAngleCPP(baseVector, currentGradients, precNorm);
//       if (currentGradients[1] < 0) {
//         angle = 360 - angle;
//       }
//       nextCell = findNextCellCPP(angle);
//       currentCell = convertCellID2IndicesCPP(i + 1, nRows, nColumns);
//       nextCell += currentCell;
//       if ((nextCell[0] > nRows) | (nextCell[0] < 1) | (nextCell[1] > nColumns) | (nextCell[1] < 1)) {
//         // if the next cell is located outside the boundaries, one can not move further
//         // into the direction of its gradient
//         visited[i] = true;
//         visitCounter++;
//         gradFieldVector[i] = vectorLength;
//       } else {
//         // otherwise, store the successor cell
//         nextCellID = convertIndices2CellIDCPP(nextCell[0], nextCell[1], nRows, nColumns);
//         cellPointer[i] = nextCellID - 1;
//       }
//     }
//   }
// 
//   /* below, the actual cumulative part begins */
// 
//   // helper variables
//   int currentCellID;
//   int pathLength;
//   IntegerVector path(n, -999);
//   for (int i = 0; i < n; i++) {
//     // iterate over all cells
//     if (visitCounter == n) {
//       // leave the for-loop if all cells have already been visited
//       break;
//     }
//     if (!visited[i]) {
//       // if the i-th cell has not yet been visited, follow the path
//       // from that cell towards a previously "visited" cell (either
//       // from an earlier path, a boundary point pointing out of bounds,
//       // or a local efficient point)
// 
//       // initialize the path with cell i
//       currentCellID = i;
//       path[0] = currentCellID;
//       pathLength = 1;
//       visited[currentCellID] = true;
//       nextCellID = cellPointer[currentCellID];
//       while (!visited[nextCellID]) {
//         // as long as the succeeding cell has not been visited, append
//         // such a cell to the current path
//         currentCellID = nextCellID;
//         path[pathLength] = currentCellID;
//         visited[currentCellID] = true;
//         visitCounter++;
//         nextCellID = cellPointer[currentCellID];
//         pathLength++;
//       }
// 
//       double a = 0.0;
//       if (gradFieldVector[nextCellID] > -999) {
//         // if the path stopped, because its next cell was already
//         // part of an earlier loop, add the value of the next
//         // cell to the cumulated sum (a) -- and ensure that
//         // this cell is not counted twice
//         a = gradFieldVector[nextCellID];
//         visitCounter--;
//       }
//       if (pathLength > 0) {
//         // iterate backwards over the path and add the cumulated
//         // sum of gradient lengths to the cells of the path
//         for (int j = pathLength; j > 0; j--) {
//           int index = path[j - 1];
//           a += gradientLengths[index];
//           gradFieldVector[index] = a;
//         }
//       }
//     }
//   }
// 
//   return gradFieldVector;
// }

// [[Rcpp::export]]
IntegerVector locallyNondominatedCPP(NumericMatrix fnMat, IntegerVector dims) {
  int n = fnMat.nrow();
  int d = dims.size();
  
  IntegerVector locallyNondominated;
  
  Rcpp::Function expandGrid("expand.grid");
  List dimensionSteps;
  
  for (int dim = 0; dim < d; dim++) {
    IntegerVector v = IntegerVector::create(-2,-1,0,1,2);
    dimensionSteps.push_back(v);
  }
  
  DataFrame deltasDF = expandGrid(dimensionSteps);
  IntegerMatrix deltas = internal::convert_using_rfunction(deltasDF, "as.matrix");
  int nDeltas = deltas.nrow();
  
  // helper
  IntegerVector indices;
  IntegerVector neighbourIndices;
  bool dominated;
  int neighbourID;
  
  for (int id = 1; id <= n; id++) {
    indices = convertCellID2IndicesCPP(id, dims);
    dominated = false;
    
    for (int r = 0; r < nDeltas; r++) {
      neighbourIndices = indices + deltas(r, _);
      if (is_true(any(neighbourIndices < 1)) ||
          is_true(any(neighbourIndices > dims))) {
        continue;
      }
      
      neighbourID = convertIndices2CellIDCPP(neighbourIndices, dims);
      
      if (is_true(all(fnMat(neighbourID - 1, _) <= fnMat(id - 1, _))) &&
          is_true(any(fnMat(neighbourID - 1, _) < fnMat(id - 1, _)))) {
        dominated = true;
        break;
      }
    }
    
    if (!dominated) {
      locallyNondominated.push_back(id);
    }
  }
  
  return locallyNondominated;
}

// [[Rcpp::export]]
NumericMatrix gridBasedGradientCPP(NumericVector fnVec, IntegerVector dims, NumericVector stepSizes, double precNorm, double precAngle) {
  int n = fnVec.size();
  int d = dims.size();
  bool scaling = (stepSizes != nullptr);
  NumericMatrix gradMat(n, d);
  
  // helper variables
  double diff;
  IntegerVector indices;
  IntegerVector iPlus;
  IntegerVector iMinus;
  double f;
  double fPlus;
  double fMinus;
  
  for (int id = 1; id <= n; id++) {
    indices = convertCellID2IndicesCPP(id, dims);
    
    for (int dim = 0; dim < d; dim++) {
      // vector access is zero-indexed!
      if (indices(dim) == 1) {
        iPlus = clone(indices);
        iPlus(dim)++;
        
        fPlus = fnVec[convertIndices2CellIDCPP(iPlus, dims) - 1];
        f = fnVec[convertIndices2CellIDCPP(indices, dims) - 1];
        
        diff = fPlus - f;
      } else if (indices(dim) == dims(dim)) {
        iMinus = clone(indices);
        iMinus(dim)--;
        
        fMinus = fnVec[convertIndices2CellIDCPP(iMinus, dims) - 1];
        f = fnVec[convertIndices2CellIDCPP(indices, dims) - 1];
        
        diff = f - fMinus;
      } else {
        iPlus = clone(indices);
        iPlus(dim)++;
        iMinus = clone(indices);
        iMinus(dim)--;
        
        fPlus = fnVec[convertIndices2CellIDCPP(iPlus, dims) - 1];
        fMinus = fnVec[convertIndices2CellIDCPP(iMinus, dims) - 1];
        
        diff = (fPlus - fMinus) / 2;
      }
      
      if (scaling) {
        diff = diff / stepSizes(dim);
      }
      
      gradMat(id-1, dim) = diff;
    }
  }
  
  return gradMat;
}


// [[Rcpp::export]]
NumericVector cumulateGradientsCPP(NumericMatrix centers, NumericMatrix gradients, double precVectorLength, double precNorm, bool fixDiagonals, bool cumulateGradientLength) {
// NumericVector cumulateHighDimensionalGradientsCPP(NumericMatrix centers, NumericMatrix gradients, double precVectorLength, double precNorm) {

  // FIXME: so far, the code only supports 2D-problems
  int d = centers.ncol();                   // dimensionality of search space
  int n = centers.nrow();                   // number of grid points
  
  // number of grid points per (search space) dimension
  IntegerVector dims(d, 0);
  for (int j = 0; j < d; j++) {
    NumericVector ctr = centers(_,j);
    ctr = unique(ctr);
    dims[j] = ctr.size();
  }

  NumericVector gradientLengths(n);         // length vector for the multi-objective gradients
  IntegerVector cellPointer(n, -999);       // vector, which indicates per cell the successor cell
  NumericVector gradFieldVector(n, -999.0); // the "final" result vector
  LogicalVector visited = rep(false, n);    // which cells have already been "visited" / "processed"
  
  // helper variables
  double vectorLength = -1.0;
  NumericVector currentGradients;
  int visitCounter = 0;                     // counter, enabling early stop of algorithm once all cells are processed
  IntegerVector currentCell(d);             // row and column index of the current cell
  IntegerVector nextCell(d);                // row and column index of the next (= successor) cell
  int nextCellID;

  for (int i = 0; i < n; i++) {
    // iterate over all cells, extract their gradients and store the gradient lengths
    currentGradients = gradients(i,_);
    vectorLength = computeVectorLengthCPP(currentGradients);
    gradientLengths[i] = vectorLength;
    if (vectorLength < precVectorLength) {
      // if the multi-objective gradient is "short enough", the current cell is considered
      // to be locally efficient (--> move on to the next cell)
      visited[i] = true;
      visitCounter++;
      gradFieldVector[i] = vectorLength;
    } else {
      // nextCell = findNextCellHighDimensionCPP(currentGradients);
      // currentCell = convertHighDimensionalCellID2IndicesCPP(i + 1, dims);
      nextCell = findNextCellCPP(currentGradients);
      currentCell = convertCellID2IndicesCPP(i + 1, dims);
      nextCell += currentCell;

      if (is_true(any(nextCell > dims)) | is_true(any(nextCell < 1))) {
        // if the next cell is located outside the boundaries, one can not move further
        // into the direction of its gradient
        visited[i] = true;
        visitCounter++;
        gradFieldVector[i] = vectorLength;
      } else {
        // otherwise, store the successor cell
        // nextCellID = convertHighDimensionalIndices2CellIDCPP(nextCell, dims);
        nextCellID = convertIndices2CellIDCPP(nextCell, dims);
        cellPointer[i] = nextCellID - 1;
      }
    }
  }

  /* below, the actual cumulative part begins */

  // helper variables
  int currentCellID;
  int pathLength;
  IntegerVector path(n, -999);
  for (int i = 0; i < n; i++) {
    // iterate over all cells
    if (visitCounter == n) {
      // leave the for-loop if all cells have already been visited
      break;
    }
    if (!visited[i]) {
      // if the i-th cell has not yet been visited, follow the path
      // from that cell towards a previously "visited" cell (either
      // from an earlier path, a boundary point pointing out of bounds,
      // or a local efficient point)

      // initialize the path with cell i
      currentCellID = i;
      path[0] = currentCellID;
      pathLength = 1;
      visited[currentCellID] = true;
      nextCellID = cellPointer[currentCellID];
      while (!visited[nextCellID]) {
        // as long as the succeeding cell has not been visited, append
        // such a cell to the current path
        currentCellID = nextCellID;
        path[pathLength] = currentCellID;
        visited[currentCellID] = true;
        visitCounter++;
        nextCellID = cellPointer[currentCellID];
        pathLength++;
      }

      double a = 0.0;
      if (gradFieldVector[nextCellID] > -999) {
        // if the path stopped, because its next cell was already
        // part of an earlier loop, add the value of the next
        // cell to the cumulated sum (a) -- and ensure that
        // this cell is not counted twice
        a = gradFieldVector[nextCellID];
        visitCounter--;
      }
      if (pathLength > 0) {
        // iterate backwards over the path and add the cumulated
        // sum of gradient lengths to the cells of the path
        for (int j = pathLength; j > 0; j--) {
          int index = path[j - 1];
          double length;
          
          if (cumulateGradientLength) {
            length = gradientLengths[index];
          } else {
            length = 1.0;
          }
          
          if(fixDiagonals) {
            NumericVector delta = as<NumericVector>(findNextCellCPP(gradients(index,_)));
            double stepSize = computeVectorLengthCPP(delta);
            a += length * stepSize;
          } else {
            a += length;
          }
          
          gradFieldVector[index] = a;
        }
      }
    }
  }

  return gradFieldVector;
}

// Multiobjective gradient for two objectives
// [[Rcpp::export]]
NumericVector getBiObjGradientCPP(NumericVector g1, NumericVector g2, double precNorm, double precAngle) {
  int len = g1.length();
  NumericVector zeros (len);
  zeros.attr("dim") = R_NilValue;
  
  g1 = normalizeVectorCPP(g1, precNorm);
  if (is_true(all(g1 == 0))) {
    // if the gradient of fn1 is zero, this has to be a local efficient point
    return(zeros);
  }
  
  g2 = normalizeVectorCPP(g2, precNorm);
  if (is_true(all(g2 == 0))) {
    // if the gradient of fn2 is zero, this has to be a local efficient point
    return(zeros);
  }
  
  double angle1 = computeAngleCPP(g1, g2, precNorm);
  
  if (abs(180 - angle1) < precAngle) {
    // if the angle between both gradients is (approximately) 180 degree,
    // this has to be a local efficient point
    return(zeros);
  }
  
  return(0.5 * (g1 + g2));
}

// Multiobjective gradient for three objectives
// [[Rcpp::export]]
NumericVector getTriObjGradientCPP(NumericVector g1, NumericVector g2, NumericVector g3, double precNorm, double precAngle) {
  int len = g1.length();
  NumericVector zeros (len);

  g1 = normalizeVectorCPP(g1, precNorm);
  g2 = normalizeVectorCPP(g2, precNorm);
  g3 = normalizeVectorCPP(g3, precNorm);
  
  if (is_true(all(g1 == 0)) ||
      is_true(all(g2 == 0)) ||
      is_true(all(g3 == 0))) {
    // if the gradient of any objective is zero, this has to be a local efficient point
    return(zeros);
  }
  
  double angle1 = computeAngleCPP(g1, g2, precNorm);
  double angle2 = computeAngleCPP(g1, g3, precNorm);
  double angle3 = computeAngleCPP(g2, g3, precNorm);
  
  if (abs(180 - angle1) < precAngle ||
      abs(180 - angle2) < precAngle ||
      abs(180 - angle3) < precAngle) {
    // if the angle between any two gradients is (approximately) 180 degree,
    // this has to be a local efficient point
    return(zeros);
  }
  
  if (abs(angle1 + angle2 + angle3 - 360) < precAngle) {
    // if all gradients show in "opposite" directions, this has to be a local effient point
    return(zeros);
  }
  
  if (len >= 3) {
    NumericVector n = normalizeVectorCPP(crossProductCPP(g1-g2,g1-g3), 0);
    NumericVector moNorm = n * sum(n * g1); // n * (n dot g1) = n * (distance to plane created by g1,g2,g3)
    
    // useful?
    // if (is_true(all(abs(moNorm) < precNorm))) {
    //   // if the norm gradient is zero, this has to be a local efficient point
    //   return(zeros);
    // }
    
    NumericMatrix gMat = cbind(g1,g2,g3);
    
    NumericVector weights;
    
    try {
      Function solve = Environment::base_env()["solve"];
      weights = solve(gMat, moNorm);
    } catch(const std::exception& e) {
      // ignore
    }
    
    if (weights != nullptr) {
      if (is_true(any(weights < 0))) {
        weights[weights < 0] = 0.0;
        weights[weights > 0] = 1.0 / sum(weights > 0);
        return weights[0] * g1 + weights[1] * g2 + weights[2] * g3;
      } else {
        return(moNorm);
      }
    } else {
      return(moNorm);
    }
  } else {
    double maxAngle = max(NumericVector::create(angle1, angle2, angle3));
    if (angle1 == maxAngle) {
      return(0.5 * (g1 + g2));
    } else if (angle2 == maxAngle) {
      return(0.5 * (g1 + g3));
    } else {
      return(0.5 * (g2 + g3));
    }
  }
  
}

// [[Rcpp::export]]
NumericMatrix getBiObjGradientGridCPP(NumericMatrix gradMat1, NumericMatrix gradMat2, double precNorm, double precAngle) {
  int n = gradMat1.rows();
  int d = gradMat1.cols();
  NumericMatrix moGradMat(n, d);
  
  for (int i = 0; i < n; i++) {
    moGradMat(i,_) = getBiObjGradientCPP(gradMat1(i,_), gradMat2(i,_), precNorm, precAngle);
  }
  
  return moGradMat;
}

// [[Rcpp::export]]
NumericMatrix getTriObjGradientGridCPP(NumericMatrix gradMat1, NumericMatrix gradMat2, NumericMatrix gradMat3, double precNorm, double precAngle) {
  int n = gradMat1.rows();
  int d = gradMat1.cols();
  NumericMatrix moGradMat(n, d);
  
  for (int i = 0; i < n; i++) {
    moGradMat(i,_) = getTriObjGradientCPP(gradMat1(i,_), gradMat2(i,_), gradMat3(i,_), precNorm, precAngle);
  }
  
  return moGradMat;
}

// [[Rcpp::export]]
NumericVector calculateMaxDisplayHeightCPP(NumericMatrix centers, NumericVector heights, bool includeDiagonals) {
  // points are "dominated" by their neighbours from some point on
  // calculate this point here
  int n = centers.nrow();
  int d = centers.ncol();
  NumericVector maxHeights(n);
  
  IntegerVector dims(d, 0);
  for (int j = 0; j < d; j++) {
    NumericVector ctr = centers(_,j);
    ctr = unique(ctr);
    dims[j] = ctr.size();
  }
  
  IntegerMatrix deltas;
  if (includeDiagonals) {
    Rcpp::Function expandGrid("expand.grid");

    IntegerVector v = IntegerVector::create(-1,0,1);
    
    DataFrame deltasDF = expandGrid(v, v, v);
    IntegerMatrix deltas = internal::convert_using_rfunction(deltasDF, "as.matrix");
  } else {
    IntegerVector v = {0,0,0,1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1};
    v.attr("dim") = Dimension(3, 7);
    deltas = as<IntegerMatrix>(v);
    deltas = transpose(deltas);
  }
  
  int nDeltas = deltas.nrow();
  
  IntegerVector indices;
  IntegerVector neighbourIndices;
  int neighbourID;
  double maxNeighbour;
  
  for (int id = 1; id <= n; id++) {
    indices = convertCellID2IndicesCPP(id, dims);
    maxNeighbour = 0.0;
    
    for (int r = 0; r < nDeltas; r++) {
      neighbourIndices = indices + deltas(r, _);
      if (is_true(any(neighbourIndices < 1)) ||
          is_true(any(neighbourIndices > dims))) {
        // keep Infinity as maxHeight when at boundaries
        maxNeighbour = std::numeric_limits<double>::infinity();
        break;
      }
      
      neighbourID = convertIndices2CellIDCPP(neighbourIndices, dims);
      
      maxNeighbour = std::max(maxNeighbour, heights(neighbourID-1));
    }
    
    maxHeights(id-1) = maxNeighbour;
  }
  
  return maxHeights;
}
