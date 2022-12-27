/*
 * MPI Program for Solving the Two-Dimensional Heat Equation.
 *
 * Author: Michael Paleos
 * Date: 12/14/2022
 * 
 * to compile: mpicc -O2 -std=c99 -o main.exe main.c functions.c -lm
 * to execute: mpiexec -n <number of processes> main.exe
 * 
 */

#include "functions.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>

#define NX 129
#define NY 129

#define LX 20.0f
#define LY 20.0f

#define DX LX / ((REAL)(NX - 1.0))
#define DY LY / ((REAL)(NY - 1.0))

#define H DX
#define K   1.0f
#define DT  0.20f * DX * DX / K

#define TIME 200
#define STEPS (INT) ((double) (TIME) / (double) (DT))

INT main(int argc, char *argv[]) {
    
    INT n_procs, rank;
    INT source, dest;
    INT start, end;

    REAL *u_total, *tmp;

    // Initializing MPI.
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&n_procs);

    INT n_dims = 1;
    INT dimension[n_dims];
    INT is_periodic[n_dims];
    INT reorder = 1;
    const INT root = 0;

    dimension[ 0 ]  = n_procs;
    is_periodic[ 0 ] = 0;

    MPI_Comm comm1D;
    MPI_Cart_create(MPI_COMM_WORLD, n_dims, dimension, is_periodic, reorder, &comm1D);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Cart_shift(comm1D, 0, 1, &source, &dest);

    MPI_Barrier(MPI_COMM_WORLD);

    // Define the number of ghost layers to account for.
    INT n_ghost_layers;
    if (rank == root || rank == n_procs - 1)
    {n_ghost_layers = 1;}

    else {n_ghost_layers = 2;}

    // Decompose the mesh.
    decomposeMesh(NY, n_procs, rank, &start, &end, n_ghost_layers);

    // Define start and end parameters.
    INT nrow = end - start + 1;
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank != root)
    {start += 1;}
    else if (rank != n_procs - 1)
    {end += 1;}

    //printf("These are the rows: %d \n", nrow);

    // Memory allocation.
    REAL *u = calloc(NX * (nrow + n_ghost_layers), sizeof(*u));
    REAL *unew = calloc(NX * (nrow + n_ghost_layers), sizeof(*unew));
    REAL *u_exact = calloc(NX * NY, sizeof(*u_exact));

    REAL *x = calloc(NX * NY, sizeof(*x));
    REAL *y = calloc(NX * NY, sizeof(*y));

    INT * recv_displs = (INT *) calloc(n_procs, sizeof(*recv_displs));
    INT * recv_counts = (INT *) calloc (n_procs, sizeof(*recv_counts));

    if (rank == root) 
    {u_total = (REAL *) calloc(NX * NY, sizeof(*u_total));} // only on the root.

    // Construct the mesh.
    meshGrid(x, y);

    MPI_Barrier(MPI_COMM_WORLD);
    REAL start_time = MPI_Wtime();

    // Starting the solution loop.
    for (INT t = 1; t < STEPS; t++) {

        // Predict the next value.
        heatEquation(unew, u, start, end, x, y);

        // Exchange information.
        exchangeSendRecv(u, start, end, source, dest, rank, n_procs);

        // Update.
        tmp = unew;
        unew = u;
        u = tmp;
    }

    // Construct the displacements and count arrays.
    for (INT i=0; i<n_procs; i++) {
        recv_displs[i] = (NY * (nrow - n_ghost_layers)) * i;
        recv_counts[i] = NY * (nrow - n_ghost_layers);
    }
    
    double *begin;
    if (rank == root)
        {begin = u;}
    else {begin = u + NY;}

    // Gather all data to an array on root.
    MPI_Gatherv(begin, NY * (nrow - n_ghost_layers), MPI_DOUBLE, 
                u_total, recv_counts, recv_displs, 
                MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    REAL end_time = MPI_Wtime();

    REAL elapsed_time = end_time - start_time;
    REAL wall_time;

    MPI_Reduce(&elapsed_time, &wall_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == root) {

        //Construct the total mesh.
        meshGrid(x, y);

        // Compute the analytical solution.
        exactSolution(u_exact, x, y);

        // Write the results.
        writeOutput(x, y, u_total);
        //writeOutputLine(x, y, u_total, u_exact);
        //writeOutputComparison(x, y, u_total, u_exact);

        free(u_total);

        // Print timing result. 
        printf("Duration = %.6f (ms) \n", wall_time * 1e3);
    }

    // Free memory.
    free(x);
    free(y);
    free(u_exact);
    free(u);
    free(unew);
    free(recv_counts);
    free(recv_displs);

    MPI_Finalize( );
    return EXIT_SUCCESS;
}