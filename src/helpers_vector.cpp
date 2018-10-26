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


// [[Rcpp::export]]
IntegerVector findNextCellCPP(double angle) {
  // which is the next cell one has to move to (based on the angle)?
  
  IntegerVector direction(2, 3);
  // ensure to only consider angles between 0 and 360 degrees
  while (angle < 360) {
    angle += 360.0;
  };
  angle = std::fmod(angle, 360.0);
  
  if ((angle > 22.5) & (angle <= 67.5)) {
    direction[0] = 1; // step in x-direction
    direction[1] = 1; // step in y-direction
  } else if ((angle > 67.5) & (angle <= 112.5)) {
    direction[0] = 1;
    direction[1] = 0;
  } else if ((angle > 112.5) & (angle <= 157.5)) {
    direction[0] = 1;
    direction[1] = -1;
  } else if ((angle > 157.5) & (angle <= 202.5)) {
    direction[0] = 0;
    direction[1] = -1;
  } else if ((angle > 202.5) & (angle <= 247.5)) {
    direction[0] = -1;
    direction[1] = -1;
  } else if ((angle > 247.5) & (angle <= 292.5)) {
    direction[0] = -1;
    direction[1] = 0;
  } else if ((angle > 292.5) & (angle <= 337.5)) {
    direction[0] = -1;
    direction[1] = 1;
  } else {
    direction[0] = 0;
    direction[1] = 1;
  }
  
  return direction ;
}



// convert rowIndex (1, ..., nRow) and columnIndex (1, ..., nColumns) to cell ID (1, ..., nRows * nColumns)
// [[Rcpp::export]]
int convertIndices2CellIDCPP(int rowIndex, int columnIndex, int nRows = 100, int nColumns = 100) {
  
  int cellID = -1;
  if ((columnIndex > nColumns) | (columnIndex < 1) | (rowIndex > nRows) | (rowIndex < 1)) {
    // if either the row- or columnIndex are located outside the boundaries, return -1 as ID
    cellID = -1;
  } else {
    cellID = (rowIndex - 1) * nColumns + columnIndex;
  }
  return cellID;
}


// convert cellID (1, ..., nRows * nColumns) to rowIndex (1, ..., nRow) and columnIndex (1, ..., nColumns)
// [[Rcpp::export]]
IntegerVector convertCellID2IndicesCPP(int cellID, int nRows = 100, int nColumns = 100) {

  IntegerVector indexVector(2, -1);
  if ((cellID > (nRows * nColumns)) | (cellID < 1)) {
    // if the cell ID has an irregular value, return -1 for the row and column indices
    indexVector = rep(-1, 2);
  } else {
    int x = cellID - 1;
    int y = floor(x / nColumns);
    indexVector[0] = y + 1;              // rowIndex
    indexVector[1] = (x % nColumns) + 1; // columnIndex
  }
  return indexVector;
}


// compute the cumulated gradient matrix based on a grid of points and their corresponding gradients;
// precVectorLength defines the threshold of the multi-objective gradient length, for which the
// corresponding point is considered to be locally efficient (suggest: 0.001 or smaller);
// precNorm is the threshold used when normalizing a vector (within the angle-computation)
// [[Rcpp::export]]
NumericVector cumulateGradientsCPP(NumericMatrix centers, NumericMatrix gradients, double precVectorLength, double precNorm) {

  // FIXME: so far, the code only supports 2D-problems
  int d = centers.ncol();                   // dimensionality of search space
  int n = centers.nrow();                   // number of grid points

  // number of grid points per (search space) dimension
  NumericVector ctr1 = centers(_,0);
  ctr1 = unique(ctr1);
  int nColumns = ctr1.size();               // #columns = number of grid points in x1-dimension
  NumericVector ctr2 = centers(_,1);
  ctr2 = unique(ctr2);
  int nRows = ctr2.size();                  // #rows = number of grid points in x2-dimension

  NumericVector gradientLengths(n);         // length vector for the multi-objective gradients
  IntegerVector cellPointer(n, -999);       // vector, which indicates per cell the successor cell
  NumericVector gradFieldVector(n, -999.0); // the "final" result vector
  LogicalVector visited = rep(false, n);    // which cells have already been "visited" / "processed"

  // helper variables
  double angle;
  double vectorLength = -1.0;
  NumericVector baseVector = NumericVector::create( 1.0, 0.0 );
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
      // compute the angle of the gradient and identify the succeeding cell
      angle = computeAngleCPP(baseVector, currentGradients, precNorm);
      if (currentGradients[1] < 0) {
        angle = 360 - angle;
      }
      nextCell = findNextCellCPP(angle);
      currentCell = convertCellID2IndicesCPP(i + 1, nRows, nColumns);
      nextCell += currentCell;
      if ((nextCell[0] > nRows) | (nextCell[0] < 1) | (nextCell[1] > nColumns) | (nextCell[1] < 1)) {
        // if the next cell is located outside the boundaries, one can not move further
        // into the direction of its gradient
        visited[i] = true;
        visitCounter++;
        gradFieldVector[i] = vectorLength;
      } else {
        // otherwise, store the successor cell
        nextCellID = convertIndices2CellIDCPP(nextCell[0], nextCell[1], nRows, nColumns);
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
          a += gradientLengths[index];
          gradFieldVector[index] = a;
        }
      }
    }
  }

  return gradFieldVector;
}

