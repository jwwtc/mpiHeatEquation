#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#ifndef SINGLE
typedef double REAL;
typedef int   INT;
#else
typedef float REAL;
typedef int    INT;
#endif

void meshGrid (REAL *x, REAL *y);
void decomposeMesh(const INT N, const INT n_procs, const INT rank, INT *start, INT *end, const INT n_ghost_layers);

void writeOutput(const REAL *x, const REAL *y, const REAL *u);
void writeOutputRank(const REAL *x, const REAL *y, const REAL *u, const REAL *u_exact, INT start, INT end, INT rank);
void writeOutputComparison(const REAL *x, const REAL *y, const REAL *u, const REAL *u_exact);
void writeOutputLine(const REAL *x, const REAL *y, const REAL *u, const REAL *u_exact);

REAL exactSolutionSum(const REAL x, const REAL y);
void exactSolution(REAL *u_exact, const REAL *x, const REAL *y);

void heatEquation(REAL *unew, REAL *u, const INT start, const INT end, const REAL *x, const REAL *y);
void exchangeSendRecv(REAL *u, const INT start, const INT end, INT source, INT dest, const INT rank, const INT n_procs);

#endif /* FUNCTIONS_H */