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
#include "defs.h"

/*****************************************************************/

/*****************************************************************
   main():
   ------

     -- Initializes, time-steps, and concludes the simulation. 
     -- Handles timing of output routines;
     -- Main is main, what more can you say.  

-*****************************************************************/
int main(int argc,char *argv[])
{
	double tdump,timage,tlog ;
	int nfailed = 0 ;

	if(argc < 3) {
		fprintf(stderr, "Insufficient arguments\nUsage: %s N1 N2\n", argv[0]);
		exit(1);
	}

	GlobalN1 = atoi(argv[1]);
	GlobalN2 = atoi(argv[2]);

	init_mpi(&argc, &argv);
	
	create_arrays();

	nstep = 0 ;

	/* Perform Initializations, either directly or via checkpoint */
	if(WorldRank == 0)
		system("mkdir -p dumps images");  
    
	if(!restart_init()) {
		init() ;
	}

	/* do initial diagnostics */
	diag(INIT_OUT) ;

	tdump = t+DTd ;
	timage = t+DTi ;
	tlog = t+DTl ;

	defcon = 1. ;
	while(t < tf) {


		/* step variables forward in time */
		nstroke = 0 ;
		step_ch() ;

        if(WorldRank == 0)
            fprintf(stderr,"%10.5g %10.5g %8d %10.5g\n",t,dt,nstep,
                    nstroke/(2.*GlobalN1*GlobalN2)) ;

		/* Handle output frequencies: */
		if(t >= tdump) {
			diag(DUMP_OUT) ;
			tdump += DTd ;
		}
		if(t >= timage) {
			diag(IMAGE_OUT) ;
			timage += DTi ;
		}
		if(t >= tlog) {
			diag(LOG_OUT) ;
			tlog += DTl ;
		}

		/* restart dump */
		nstep++ ;
		if(nstep%DTr == 0)
			restart_write() ; // this is a barrier


		/* deal with failed timestep, though we usually exit upon failure */
		if (failed) {
			restart_init() ;
			failed = 0 ;
			nfailed = nstep ;
			defcon = 0.3 ;
		}
		if(nstep > nfailed + DTr*4.*(1 + 1./defcon)) defcon = 1. ;

	}

    if(WorldRank == 0)
	    fprintf(stderr,"ns,ts: %d %d\n",nstep,nstep*N1*N2) ;

	/* do final diagnostics */
	diag(FINAL_OUT) ;

  clean_up();
	return(0) ;
}


/*****************************************************************/
/*****************************************************************
  set_arrays():
  ----------

       -- sets to zero all arrays
 *****************************************************************/
void set_arrays()
{
	int i,j,k ;

	/* everything must be initialized to zero */
	ZSLOOP(-2,N1+1,-2,N2+1) {
		PLOOP {
			p[i][j][k]   = 0. ;
			ph[i][j][k]  = 0. ;
			dq[i][j][k]  = 0. ;
			F1[i][j][k]  = 0. ;
			F2[i][j][k]  = 0. ;
		}
		pflag[i][j] = 0 ;
	}

	k = 0;
	IMAGELOOP {  
	    failimage[0][k] = failimage[1][k] = failimage[2][k] = failimage[3][k] = failimage[4][k++] = 0 ; 
	}


}


/*****************************************************************/
/*****************************************************************
  set_grid():
  ----------

       -- calculates all grid functions that remain constant 
          over time, such as the metric (gcov), inverse metric 
          (gcon), connection coefficients (conn), and sqrt of 
          the metric's determinant (gdet).

 *****************************************************************/
void set_grid()
{
	int i,j,k ;
	double X[NDIM] ;
	struct of_geom geom ;

	/* set up boundaries, steps in coordinate grid */
	set_points() ;
	dV = dx[1]*dx[2] ;

	DLOOPA X[j] = 0. ;

	ZSLOOP(-2,N1+1,-2,N2+1) {
		
		/* zone-centered */
		coord(i,j,CENT,X) ;
		gcov_func(X, gcov[i][j][CENT]) ;
		gdet[i][j][CENT] = gdet_func(gcov[i][j][CENT]) ;
		gcon_func(gcov[i][j][CENT], gcon[i][j][CENT]) ;

		get_geometry(i,j,CENT,&geom) ;
		conn_func(X, &geom, conn[i][j]) ;

		/* corner-centered */
		coord(i,j,CORN,X) ;
		gcov_func(X,gcov[i][j][CORN]) ;
		gdet[i][j][CORN] = gdet_func(gcov[i][j][CORN]) ;

		gcon_func(gcov[i][j][CORN],gcon[i][j][CORN]) ;

		/* r-face-centered */
		coord(i,j,FACE1,X) ;
		gcov_func(X,gcov[i][j][FACE1]) ;
		gdet[i][j][FACE1] = gdet_func(gcov[i][j][FACE1]) ;
		gcon_func(gcov[i][j][FACE1],gcon[i][j][FACE1]) ;

		/* theta-face-centered */
		coord(i,j,FACE2,X) ;
		gcov_func(X,gcov[i][j][FACE2]) ;
		gdet[i][j][FACE2] = gdet_func(gcov[i][j][FACE2]) ;
		gcon_func(gcov[i][j][FACE2],gcon[i][j][FACE2]) ;
	}

    // Halo exchange: gcov, gcon, conn, gdet
    // No need here as the values are rightly initialized here for all nodes

	/* done! */
}

/*****************************************************************
  init_mpi():
  ----------

       -- initializes MPI, creates 2D grid topology,
            row and column communicators, and divides N1 and N2.

 *****************************************************************/
void init_mpi(int *argc, char*** argv)
{
	int obtained, dimsize[2], periodic[2], coords[2];

	// Initialize MPI
	MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &obtained);
	assert(obtained == MPI_THREAD_FUNNELED);

	// Obtain number of processes and rank
	MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &WorldRank);

	// Create 2D grid topology
	dimsize[0] = dimsize[1] = 0;
	MPI_Dims_create(NumProcs, 2, dimsize);
	periodic[0] = periodic[1] = 0;	// Change if needed later on
	MPI_Cart_create(MPI_COMM_WORLD, 2, dimsize, periodic, 0, &Comm2D);
	MPI_Cart_coords(Comm2D, WorldRank, 2, coords);
  NumRows = dimsize[0];
  NumCols = dimsize[1];
	RowRank = coords[0];
	ColRank = coords[1];

	// Create separate column and row communicators	
	MPI_Comm_split(Comm2D, RowRank, ColRank, &CommRow);
	MPI_Comm_split(Comm2D, ColRank, RowRank, &CommCol);

	// Divide N1 and N2 for each process
	N1 = GlobalN1/NumRows;
	N2 = GlobalN2/NumCols;

  MPI_Type_vector(N1+4, 2*NPR, (N2+4)*NPR, MPI_DOUBLE, &d_col_type);
  MPI_Type_commit(&d_col_type);

  MPI_Type_contiguous(2*(N2+4)*NPR, MPI_DOUBLE, &d_row_type);
  MPI_Type_commit(&d_row_type);

  MPI_Type_vector(N1+4, 2, N2+4, MPI_INT, &i_col_type);
  MPI_Type_commit(&i_col_type);

  MPI_Type_contiguous(2*(N2+4), MPI_INT, &i_row_type);
  MPI_Type_commit(&i_row_type);

  halo_count = 0;

  // Custom types for printing restart files
  int globalsizes[2] = {(N1+4)*NumRows, (N2+4)*NumCols};
  int localsizes[2]  = {N1+4, N2+4};
  int starts[2]      = {(N1+4)*RowRank, (N2+4)*ColRank};
  int order          = MPI_ORDER_C;

  MPI_Type_contiguous(NPR*CHARSPERNUM, MPI_CHAR, &array_as_string);
  MPI_Type_commit(&array_as_string);
  MPI_Type_create_subarray(2, globalsizes, localsizes, starts, order,
                           array_as_string, &local_array);
  MPI_Type_commit(&local_array);

  // Custom types for creating images
  globalsizes[0] = GlobalN1; globalsizes[1] = GlobalN2;
  localsizes[0] = N1; localsizes[1] = N2;
  starts[0] = N1*RowRank; starts[1] = N2*ColRank;
  order = MPI_ORDER_FORTRAN;

  MPI_Type_contiguous(3, MPI_CHAR, &rgb);
  MPI_Type_commit(&rgb);
  MPI_Type_create_subarray(2, globalsizes, localsizes, starts, order, rgb,
      &color_array);
  MPI_Type_commit(&color_array);
}

/*****************************************************************
  create_arrays():
  ----------

    -- allocates all the necessary arrays,
          plus performs pointer trick 
          so that grid arrays can legitimately refer to ghost 
          zone quantities at positions  i = -2, -1, N1, N1+1 and 
          j = -2, -1, N2, N2+1
 *****************************************************************/
void create_arrays()
{
	int i, j;
	D_Alloc_2D_1(a_p, N1+4, N2+4, NPR);
	D_Alloc_2D_1(a_dq, N1+4, N2+4, NPR);
	D_Alloc_2D_1(a_F1, N1+4, N2+4, NPR);
	D_Alloc_2D_1(a_F2, N1+4, N2+4, NPR);
	D_Alloc_2D_1(a_ph, N1+4, N2+4, NPR);
	I_Alloc_2D(a_pflag, N1+4, N2+4);

	D_Alloc_2D_1(psave, N1, N2, NPR);
	D_Alloc_2D(fimage, NIMG, N1*N2);
	I_Alloc_2D(failimage, 5, N1*N2);

	D_Alloc_2D_3(a_conn, N1+4, N2+4, NDIM, NDIM, NDIM);
	D_Alloc_2D_3(a_gcon, N1+4, N2+4, NPG, NDIM, NDIM);
	D_Alloc_2D_3(a_gcov, N1+4, N2+4, NPG, NDIM, NDIM);
	D_Alloc_2D_1(a_gdet, N1+4, N2+4, NPG);

#if(DO_FONT_FIX)
	D_Alloc_1D(Katm, N1);
#endif

	D_Alloc_2D(emf, N1+1, N2+1);
    
    // Specially allocate for p, dq, F1, F2, ph, pflag
    AD_Alloc_2D(p, a_p, N1+4, 2, 2, p_tmp);
    AD_Alloc_2D(dq, a_dq, N1+4, 2, 2, dq_tmp);
    AD_Alloc_2D(F1, a_F1, N1+4, 2, 2, F1_tmp);
    AD_Alloc_2D(F2, a_F2, N1+4, 2, 2, F2_tmp);
    AD_Alloc_2D(ph, a_ph, N1+4, 2, 2, ph_tmp);
    AI_Alloc_2D(pflag, a_pflag, N1+4, 2, 2, pflag_tmp);
	
    /* grid functions */
    AD_Alloc_2D(conn, a_conn, N1+4, 2, 2, conn_tmp);
    AD_Alloc_2D(gcon, a_gcon, N1+4, 2, 2, gcon_tmp);
    AD_Alloc_2D(gcov, a_gcov, N1+4, 2, 2, gcov_tmp);
    AD_Alloc_2D(gdet, a_gdet, N1+4, 2, 2, gdet_tmp);

}

/*****************************************************************
  clean_up():
  ----------

       -- finalize MPI, free arrays.

 *****************************************************************/
void clean_up()
{
    Free_Offset(p, -2);
    Free_Offset(dq, -2);
    Free_Offset(F1, -2);
    Free_Offset(F2, -2);
    Free_Offset(ph, -2);
    Free_Offset(pflag, -2);
    
    Free_Offset(conn, -2);
    Free_Offset(gcon, -2);
    Free_Offset(gcov, -2);
    Free_Offset(gdet, -2);
    
    Free_2D(a_p);
    Free_2D(a_dq);
    Free_2D(a_F1);
    Free_2D(a_F2);
    Free_2D(a_ph);
    Free_2D(a_pflag);
    
    Free_2D(psave);
    Free_2D(fimage);
    Free_2D(failimage);
    Free_2D(a_conn);
    Free_2D(a_gcon);
    Free_2D(a_gcov);
    Free_2D(a_gdet);

    MPI_Type_free(&array_as_string);
    MPI_Type_free(&local_array);
    MPI_Type_free(&rgb);
    MPI_Type_free(&color_array);

    MPI_Finalize();
}
