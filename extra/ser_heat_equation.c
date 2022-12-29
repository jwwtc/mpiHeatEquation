/*
 * Serial Implementation of the 2D Heat Conduction Equation.
 *
 * Author: Michael Paleos
 * Date: 12/07/2022
 *
 * to compile: gcc -O2 -std=c99 ser_heat_equation.c -lm -o ser_heat.exe
 * to execute: ./ser_heat.exe
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>

#define LX 20.0f
#define LY 20.0f

#define NX 11
#define NY 11

#define DX LX / ((REAL)(NX - 1.0))
#define DY LY / ((REAL)(NY - 1.0))
#define H DX

#define K   1.0f
#define DT  0.20f * DX * DX / K

#define STEPS (INT) ((double) (TIME) / (double) (DT))
#define TIME 200

#ifndef SINGLE
typedef double REAL;
typedef int   INT;
#else
typedef float REAL;
typedef int    INT;
#endif

#define PI M_PI
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

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
    output = fopen("results/temp.csv", "w+");
    fprintf(output,"X,Y,Temperature\n");

    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            printf("(X, Y, u_num): %d %d %d \n\n", x[i*NX+j], y[i*NX+j], u[i*NX+j]);
            fprintf(output, "%8f,%8f,%8f\n",
            x[i*NX+j],y[i*NX+j],u[i*NX+j]);
        }
    }
    fclose(output);
}

// Function to write the comparison results.
void writeOutputComparison(const REAL *x, const REAL *y, const REAL *u, const REAL *uExact)
{
    INT   i, j;
    FILE *output;
    output = fopen("results/temp_comp.csv", "w+");
    fprintf(output,"X,Y,Temp_Num,Temp_Exact\n");

    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            fprintf(output, "%8f,%8f,%8f,%8f\n", 
            x[i*NY+j],y[i*NY+j],u[i*NY+j],uExact[i*NY+j]);
        }
    }
    fclose(output);
}

// Function to write the comparison results.
void writeOutputLine(const REAL *x, const REAL *y, const REAL *u, const REAL *uExact)
{
    INT   i, j;
    FILE *output;

    output = fopen("results/temp_line.csv", "w+");
    fprintf(output,"X,Y,Temp_Num,Temp_Exact\n");

    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {

        if (y[i*NX+j] == LY / 2)
            {fprintf(output, "%8f,%8f,%8f,%8f\n", 
            x[i*NY+j],y[i*NY+j],u[i*NY+j],uExact[i*NY+j]);}

        else if (x[i*NX+j] == LX / 2)
            {fprintf(output, "%8f,%8f,%8f,%8f\n", 
            x[i*NY+j],y[i*NY+j],u[i*NY+j],uExact[i*NY+j]);}
        }
    }
    fclose(output);
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
void exactSolution(REAL *uExact, const REAL *x, const REAL *y)
{
    for (INT i = 0; i < NX - 1; i++) {
        for (INT j = 0; j < NY - 1; j++) {

            REAL x = i * DX;
            REAL y = j * DY;
            REAL y2 = (j+1) * H;

            if (y2 == LY)
                { uExact[i*NY+(j+1)] = 100.0f;
                uExact[i*NY+j] = exactSolutionSum(x, y); }

            else{ uExact[i*NY+j] = exactSolutionSum(x, y); }

        }
    }
}

// Function to predict the new temperature values.
void heatEquation(REAL *unew, REAL *u, INT i, INT j, const REAL *x, const REAL *y) 
{
    // Set boundary conditions.
    if (y[i*NX+j] == LY)
        {u[i*NX+j] = 100; unew[i*NX+j] = 100;}

    else if (x[i*NX+j] == 0 || x[i*NX+j] == LX || y[i*NX+j] == 0)
        {u[i*NX+j] = 0; unew[i*NX+j] = 0;}

    // Solve for the other values.
    else {
        unew[i*NY+j] = 
        u[i*NY+j] + (DT / (H * H)) *
        ( u[(i+1)*NY+j] + u[(i-1)*NY+j] +
        u[i*NY+(j+1)] + u[i*NY+(j-1)] - 4.0f*u[i*NY+j] );
    }
}

INT main() {

    REAL *uExact, *x, *y;
    REAL *u, *unew, *tmp;

    // Allocate memory.
    unew   = calloc(NX * NY, sizeof(*unew));
    u      = calloc(NX * NY, sizeof(*u));
    x      = calloc(NX * NY, sizeof(*x));
    y      = calloc(NX * NY, sizeof(*y));
    uExact = calloc(NX * NY, sizeof(*uExact));

    // Construct the mesh.
    meshGrid(x, y);

    // Compute the analytical solution.
    exactSolution(uExact, x, y);

    // Solving for every time step.
    for (INT t = 1; t < STEPS; t++) {

        for (INT i=1; i < NX; i++) {
            for (INT j=1; j < NY; j++) {

            // Predict the next value.
            heatEquation(unew, u, i, j, x, y);

            }
        }

        // Update the values.
        tmp = unew;
        unew = u;
        u = tmp;

    }

    // Write the results.
    writeOutput(x, y, u);
    //writeOutputComparison(x, y, u, uExact);
    //writeOutputLine(x, y, u, uExact);

    // Free memory.
    free(x);
    free(y);
    free(u);
    free(unew);
    free(uExact);

    return EXIT_SUCCESS;
}
