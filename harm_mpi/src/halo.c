/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

#include "decs.h"

/* Halo exchange of 2 layers for double matrix of size [N1+4][N2+4][NPR] */
// Pointer is the displaced one with [0][0] corresponding to [2][2]
void halo_npr(double (** pv)[NPR]) 
{
    MPI_Request req[4];
    MPI_Status sts[4];
    int num_req;
    // First exchange rows
    // This ensures diagonal consistency without explicit communication
    num_req = 0;
    if(RowRank > 0) {
        // Receive rows -2 & -1 from (RowRank - 1) on CommCol
        MPI_Irecv(pv[-2][-2], 1, d_row_type, (RowRank - 1),
                halo_count, CommCol, &req[num_req++]);
    }
    if(RowRank < (NumRows - 1)) {
        // Receive rows N1 & N1+1 from (RowRank + 1) on CommCol
        MPI_Irecv(pv[N1][-2], 1, d_row_type, (RowRank + 1),
                halo_count, CommCol, &req[num_req++]);
    }

    if(RowRank > 0) {
        // Send rows 0 & 1 to (RowRank - 1) on CommCol
        MPI_Isend(pv[0][-2], 1, d_row_type, (RowRank - 1),
                halo_count, CommCol, &req[num_req++]);
    }
    if(RowRank < (NumRows - 1)) {
        // Send rows 0 & 1 to (RowRank - 1) on CommCol
        MPI_Isend(pv[N1-2][-2], 1, d_row_type, (RowRank + 1),
                halo_count, CommCol, &req[num_req++]);
    }
    MPI_Waitall(num_req, req, sts);

    // Then exchange columns
    // This is slightly tricky and needs MPI data types
    // This ensures diagonal consistency without explicit communication
    num_req = 0;
    if(ColRank > 0) {
        // Receive cols -2 & -1 from (ColRank - 1) on CommRow
        MPI_Irecv(pv[-2][-2], 1, d_col_type, (ColRank - 1),
                halo_count, CommRow, &req[num_req++]);
    }
    if(ColRank < (NumCols - 1)) {
        // Receive cols N2 & N2+1 from (ColRank + 1) on CommRow
        MPI_Irecv(pv[-2][N2], 1, d_col_type, (ColRank + 1),
                halo_count, CommRow, &req[num_req++]);
    }
    
    if(ColRank > 0) {
        // Send cols 0 & 1 to (ColRank - 1) on CommRow
        MPI_Isend(pv[-2][0], 1, d_col_type, (ColRank - 1),
                halo_count, CommRow, &req[num_req++]);
    }
    if(ColRank < (NumCols - 1)) {
        // Send cols N2-2 & N2-1 to (ColRank + 1) on CommRow
        MPI_Isend(pv[-2][N2-2], 1, d_col_type, (ColRank + 1),
                halo_count, CommRow, &req[num_req++]);
    }
    MPI_Waitall(num_req, req, sts);

    halo_count++;
}

/* Halo exchange of 2 layers of pflag */
// Pointer is the displaced one with [0][0] corresponding to [2][2]
void halo_pflag() 
{
    MPI_Request req[4];
    MPI_Status sts[4];
    int num_req;
    // First exchange rows
    // This ensures diagonal consistency without explicit communication
    num_req = 0;
    if(RowRank > 0) {
        // Receive rows -2 & -1 from (RowRank - 1) on CommCol
        MPI_Irecv(&pflag[-2][-2], 1, d_row_type, (RowRank - 1),
                halo_count, CommCol, &req[num_req++]);
    }
    if(RowRank < (NumRows - 1)) {
        // Receive rows N1 & N1+1 from (RowRank + 1) on CommCol
        MPI_Irecv(&pflag[N1][-2], 1, d_row_type, (RowRank + 1),
                halo_count, CommCol, &req[num_req++]);
    }

    if(RowRank > 0) {
        // Send rows 0 & 1 to (RowRank - 1) on CommCol
        MPI_Isend(&pflag[0][-2], 1, d_row_type, (RowRank - 1),
                halo_count, CommCol, &req[num_req++]);
    }
    if(RowRank < (NumRows - 1)) {
        // Send rows 0 & 1 to (RowRank - 1) on CommCol
        MPI_Isend(&pflag[N1-2][-2], 1, d_row_type, (RowRank + 1),
                halo_count, CommCol, &req[num_req++]);
    }
    MPI_Waitall(num_req, req, sts);

    // Then exchange columns
    // This is slightly tricky and needs MPI data types
    // This ensures diagonal consistency without explicit communication
    num_req = 0;
    if(ColRank > 0) {
        // Receive cols -2 & -1 from (ColRank - 1) on CommRow
        MPI_Irecv(&pflag[-2][-2], 1, d_col_type, (ColRank - 1),
                halo_count, CommRow, &req[num_req++]);
    }
    if(ColRank < (NumCols - 1)) {
        // Receive cols N2 & N2+1 from (ColRank + 1) on CommRow
        MPI_Irecv(&pflag[-2][N2], 1, d_col_type, (ColRank + 1),
                halo_count, CommRow, &req[num_req++]);
    }
    
    if(ColRank > 0) {
        // Send cols 0 & 1 to (ColRank - 1) on CommRow
        MPI_Isend(&pflag[-2][0], 1, d_col_type, (ColRank - 1),
                halo_count, CommRow, &req[num_req++]);
    }
    if(ColRank < (NumCols - 1)) {
        // Send cols N2-2 & N2-1 to (ColRank + 1) on CommRow
        MPI_Isend(&pflag[-2][N2-2], 1, d_col_type, (ColRank + 1),
                halo_count, CommRow, &req[num_req++]);
    }
    MPI_Waitall(num_req, req, sts);

    halo_count++;
}
/*
void print_npr(double (** pv)[NPR])
{
    int i, j, k;
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    for(k = 0; k < NumProcs; k++) {
        if(WorldRank == k) {
            ZSLOOP(-2,N1+1,-2,N2+1) {
                fflush(stdout);
                printf("%d:pv[%d][%d] = %g\n",k,N1*RowRank + i, N2*ColRank + j, pv[i][j][RHO]);
                fflush(stdout);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
*/
/*void print_npr(double (** pv)[NPR], int val)
{
    int i, j, r, c;
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    for(r = 0; r < NumRows; r++) {
        if(RowRank == r) {
            for(i = -2; i < N1+2; i++) {
                for(c = 0; c < NumCols; c++) {
                    if(ColRank == c) {
                        for(j = -2; j < N2+2; j++) {
                            fflush(stdout);
                            printf("%0.4e  ", pv[i][j][val]);
                            fflush(stdout);
                        }
                        fflush(stdout);
                        printf("\t");
                        fflush(stdout);
                    }
                    fflush(stdout);
                    MPI_Barrier(CommRow);
                    fflush(stdout);
                }
                if(ColRank == 0) {
                    fflush(stdout);
                    printf("\n");
                    fflush(stdout);
                }
                fflush(stdout);
                MPI_Barrier(CommRow);
                fflush(stdout);
            }
            if(ColRank == 0) {
                fflush(stdout);
                printf("\n");
                fflush(stdout);
            }
            fflush(stdout);
            MPI_Barrier(CommRow);
            fflush(stdout);
        }
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
        fflush(stdout);
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    fflush(stdout);
}
*/

void print_npr(double (** pv)[NPR], int val) {
    double (*** tmp)[NPR];
    int i, j, r, c;

    MPI_Status sts;

    if(WorldRank == 0) {
        tmp = malloc (sizeof(double *) * NumCols);
        for(j = 0; j < NumCols; j++) {
            D_Alloc_2D_1(tmp[j], N1+4, N2+4, NPR);
        }
    }
    fflush(stderr);
    fflush(stdout);

    for(r= 0; r < NumRows; r++) {
        for(c = 0; c < NumCols; c++) {
            if((RowRank == r) && (ColRank == c)) {
                if(WorldRank == 0) {
                    for(i = 0; i < N1+4; i++) {
                        for(j = 0; j < N2+4; j++) {
                            tmp[c][i][j][val] = pv[i-2][j-2][val];
                        }
                    }

                } else {
                    MPI_Send(pv[-2][-2], (N1+4)*(N2+4)*(NPR), MPI_DOUBLE, 0, r*NumCols + c, MPI_COMM_WORLD);
                }
            }
            else if(WorldRank == 0) {
                MPI_Recv(tmp[c][0][0], (N1+4)*(N2+4)*(NPR), MPI_DOUBLE, r*NumCols + c, r*NumCols + c, MPI_COMM_WORLD, &sts);
            }
        }
        if(WorldRank == 0) {
            fflush(stdout);
            for(i = 0; i < N1+4; i++) {
                for(c = 0; c < NumCols; c++) {
                    for(j = 0; j < N2+4; j++) {
                            printf("%e\t", tmp[c][i][j][val]);
                    }
                    printf("\t");
                }
                printf("\n");
            }
            printf("\n");
            fflush(stdout);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(WorldRank == 0) {
        for(i = 0; i < NumCols; i++) {
            Free_2D(tmp[i]);
        }
        free(tmp);
    }
}

void print_npr2(double (** pv)[NPR], int val) {
    double (*** tmp)[NPR];
    int i, j, r, c;

    MPI_Status sts;

    if(WorldRank == 0) {
        tmp = malloc (sizeof(double *) * NumCols);
        for(j = 0; j < NumCols; j++) {
            D_Alloc_2D_1(tmp[j], N1+4, N2+4, NPR);
        }
    }
    fflush(stderr);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    for(r= 0; r < NumRows; r++) {
        for(c = 0; c < NumCols; c++) {
            if((RowRank == r) && (ColRank == c)) {
                if(WorldRank == 0) {
                    for(i = 0; i < N1+4; i++) {
                        for(j = 0; j < N2+4; j++) {
                            tmp[c][i][j][val] = pv[i-2][j-2][val];
                        }
                    }

                } else {
                    MPI_Send(pv[-2][-2], (N1+4)*(N2+4)*(NPR), MPI_DOUBLE, 0, r*NumCols + c, MPI_COMM_WORLD);
                }
            }
            else if(WorldRank == 0) {
                MPI_Recv(tmp[c][0][0], (N1+4)*(N2+4)*(NPR), MPI_DOUBLE, r*NumCols + c, r*NumCols + c, MPI_COMM_WORLD, &sts);
            }
        }
        if(WorldRank == 0) {
            fflush(stdout);
            if(r == 0) {
                // Print first 2 rows
                for(i = 0; i < 2; i++) {
                    // Print first 2 columns
                    printf("%e\t%e\t", tmp[0][i][0][val], tmp[0][i][1][val]);
                    // Print N2 columns
                    for(c = 0; c < NumCols; c++) {
                        for(j = 2; j < N2+2; j++) {
                            printf("%e\t", tmp[c][i][j][val]);
                        }
                    }
                    // Print last 2 columns
                    printf("%e\t%e\n", tmp[NumCols-1][i][N2+2][val], tmp[NumCols-1][i][N2+3][val]);
                }
            }
            // Print N1 rows
            for(i = 2; i < N1+2; i++) {
                // Print first 2 columns
                printf("%e\t%e\t", tmp[0][i][0][val], tmp[0][i][1][val]);
                // Print N2 columns
                for(c = 0; c < NumCols; c++) {
                    for(j = 2; j < N2+2; j++) {
                        printf("%e\t", tmp[c][i][j][val]);
                    }
                }
                // Print last 2 columns
                printf("%e\t%e\n", tmp[NumCols-1][i][N2+2][val], tmp[NumCols-1][i][N2+3][val]);
            }

            if(r == (NumRows -1)){
                // Print last 2 rows
                for(i = N1+2; i < N1+4; i++) {
                    // Print first 2 columns
                    printf("%e\t%e\t", tmp[0][i][0][val], tmp[0][i][1][val]);
                    // Print N2 columns
                    for(c = 0; c < NumCols; c++) {
                        for(j = 2; j < N2+2; j++) {
                            printf("%e\t", tmp[c][i][j][val]);
                        }
                    }
                    // Print last 2 columns
                    printf("%e\t%e\n", tmp[NumCols-1][i][N2+2][val], tmp[NumCols-1][i][N2+3][val]);
                }
                printf("\n");
            }
            fflush(stdout);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(WorldRank == 0) {
        for(i = 0; i < NumCols; i++) {
            Free_2D(tmp[i]);
        }
        free(tmp);
    }
}

