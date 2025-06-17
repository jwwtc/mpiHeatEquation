#include "functions.h"
#include "constants.h"

#include <mpi.h>
#include <math.h>
#include <stdio.h>

/* Use PI definition from constants.h */

// Function to construct the mesh.
void meshGrid (REAL *x, REAL *y)
{
    INT i, j;
    for (i=0; i<NX; i++) {
        for (j=0; j<NY; j++){
            x[i*NY+j] =  DX * ((REAL) i);
            y[i*NY+j] =  DY * ((REAL) j);
        }
    }
}

// Function to write the final results.
void writeOutput(const REAL *x, const REAL *y, const REAL *u)
{
    INT   i, j;
    FILE *output;
    output = fopen("./../data/par_temp.csv", "w+");
    fprintf(output,"X,Y,Temperature\n");

    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            fprintf(output, "%8f,%8f,%8f\n",
            x[i*NY+j],y[i*NY+j],u[i*NY+j]);
        }
    }
    fclose(output);
}

// Function to write the final results of a rank.
void writeOutputRank(const REAL *x, const REAL *y, 
                    const REAL *u, const REAL *u_exact, 
                    INT start, INT end, INT rank)
{
    INT   i, j;
    FILE *output;
    output = fopen("./../data/par_temp_rank_comp.csv", "w+");
    fprintf(output,"X,Y,Numerical,Exact,Rank    %d\n", rank);

    for (i = start; i < end; i++) {
        for (j = 0; j < NY; j++) {
            fprintf(output, "%8f,%8f,%8f,%f\n",
            x[i*NY+j],y[i*NY+j],u[i*NY+j],u_exact[i*NY+j]);
        }
    }
    fclose(output);
}

// Function to write the comparison results.
void writeOutputComparison(const REAL *x, const REAL *y, const REAL *u, const REAL *u_exact)
{
    INT   i, j;
    FILE *output;
    output = fopen("./../data/par_temp_comp.csv", "w+");
    fprintf(output,"X,Y,Temp_Num,Temp_Exact\n");

    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            //printf("(X, Y, u_num): %d %d %d \n\n", x[i*NY+j], y[i*NY+j], u[i*NY+j]);
            fprintf(output, "%8f,%8f,%8f,%8f\n", 
            x[i*NY+j],y[i*NY+j],u[i*NY+j],u_exact[i*NY+j]);
        }
    }
    fclose(output);
}

// Function to write the results for the line plots.
void writeOutputLine(const REAL *x, const REAL *y, const REAL *u, const REAL *u_exact)
{
    INT   i, j;
    FILE *output;

    output = fopen("./../data/par_temp_line.csv", "w+");
    fprintf(output,"X,Y,Temp_Num,Temp_Exact\n");

    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {

        if (y[i*NY+j] == LY / 2)
            {fprintf(output, "%8f,%8f,%8f,%8f\n", 
            x[i*NY+j],y[i*NY+j],u[i*NY+j],u_exact[i*NY+j]);}

        else if (x[i*NY+j] == LX / 2)
            {fprintf(output, "%8f,%8f,%8f,%8f\n", 
            x[i*NY+j],y[i*NY+j],u[i*NY+j],u_exact[i*NY+j]);}
        }
    }
    fclose(output);
}

// Function to decompose the mesh.
void decomposeMesh(const INT N, const INT n_procs, const INT rank, INT *start, INT *end,
                      const INT n_ghost_layers)
{
    INT remainder = N % n_procs;
    if (remainder == 0) {
        *start = 0;
        *end   = (N / n_procs) + n_ghost_layers - 1;
    } else {
        *start = 0;
        INT pointsPerProcess = (N - remainder) / n_procs + 1;
        if (rank == (n_procs - 1))
            *end = (N - pointsPerProcess * (n_procs - 1)) + n_ghost_layers - 1;
        else
            *end = pointsPerProcess + n_ghost_layers - 1;
    }
}

// Complimentary function to the exactSolution function.
REAL exactSolutionSum(const REAL x, const REAL y)
{
    REAL top_temp = 100.0f;
    INT n_prec = 95;

    REAL sum = 0;
    REAL result = 0;

    for (INT n = 1; n < n_prec; n += 2){

        result = (4.0f * top_temp) / (n * M_PI) * sin( (n * M_PI * x) / LX )
        * sinh ( (n * M_PI * y) / LY ) / sinh(n * M_PI);

        sum = sum + result;
    }

    return sum;

}

// Function to calculate the analytical solution.
void exactSolution(REAL *u_exact, const REAL *x, const REAL *y)
{
    for (INT i = 0; i < NX - 1; i++) {
        for (INT j = 0; j < NY - 1; j++) {

            REAL x = i * DX;
            REAL y = j * DY;
            REAL y2 = (j+1) * H;

            if (y2 == LY)
                { u_exact[i*NY+(j+1)] = 100.0f;
                u_exact[i*NY+j] = exactSolutionSum(x, y); }

            else{ u_exact[i*NY+j] = exactSolutionSum(x, y); }

        }
    }
}

// Function to predict the new temperature values.
void heatEquation(REAL *unew, REAL *u, const INT start, const INT end, const REAL *x, const REAL *y)
{           
    for (INT i = start; i < end; i++) {
        for (INT j = 1; j < NY; j++) {

            // Set boundary conditions.
            if (y[i*NY+j] == LY)
                {u[i*NY+j] = 100; unew[i*NY+j] = 100;}

            else if (x[i*NY+j] == 0 || x[i*NY+j] == LX || y[i*NY+j] == 0)
                {u[i*NY+j] = 0; unew[i*NY+j] = 0;}

            // Solve for the other coordinates.
            else {
                unew[i*NY+j] =
                u[i*NY+j] + (DT / (H * H)) *
                ( u[(i+1)*NY+j] + u[(i-1)*NY+j] +
                u[i*NY+(j+1)] + u[i*NY+(j-1)] - 4.0f*u[i*NY+j] );
            }
        }
    }
}

// Function for sendrecv exchange.
void exchangeSendRecv(REAL *u, const INT start, const INT end, INT source, INT dest, const INT rank, const INT n_procs)
{
    INT tag0 = 0;
    INT tag1 = 1;
    INT e = (end - 1) * NX;
    INT s = 0;

    MPI_Sendrecv(&u[ e ], NX, MPI_DOUBLE, dest, tag0, &u[ s ], NX, MPI_DOUBLE, source, 
    tag0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    s = (start + 1) * NX;
    e = end * NX;

    // swaping the source and destination.
    INT tmp = dest;
    dest = source;
    source = tmp;

    MPI_Sendrecv(&u[ s ], NX, MPI_DOUBLE, dest, tag1, &u[ e ], NX, MPI_DOUBLE,
    source, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
