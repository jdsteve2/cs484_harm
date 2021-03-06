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

/**
 *
 * this contains the generic piece of code for advancing
 * the primitive variables 
 *
**/


#include "decs.h"
#include <omp.h>
#include <stdbool.h>
#include "support.h"


/** algorithmic choices **/

/* use local lax-friedrichs or HLL flux:  these are relative weights on each numerical flux */
#define HLLF  (0.0)
#define LAXF  (1.0)

/** end algorithmic choices **/


double advance( double pi[][N2+4][NPR], double pb[][N2+4][NPR], 
	double Dt, double pf[][N2+4][NPR]) ;
double fluxcalc( double pr[][N2+4][NPR], double F[][N2+4][NPR], 
	int dir ) ;
void   flux_cd(double F1[][N2+4][NPR], double F2[][N2+4][NPR]) ;


/***********************************************************************************************/
/***********************************************************************************************
  step_ch():
  ---------
     -- handles the sequence of making the time step, the fixup of unphysical values, 
        and the setting of boundary conditions;

     -- also sets the dynamically changing time step size;

***********************************************************************************************/
void step_ch()
{
	double ndt ;
	int i,j,k ;

	fprintf(stderr,"h") ;
	ndt = advance(p, p, 0.5*dt, ph) ;   /* time step primitive variables to the half step */

	fixup(ph) ;         /* Set updated densities to floor, set limit for gamma */ // added omp calls
	bound_prim(ph) ;    /* Set boundary conditions for primitive variables, flag bad ghost zones */  // added omp calls
	fixup_utoprim(ph);  /* Fix the failure points using interpolation and updated ghost zone values */ // added omp calls
	bound_prim(ph) ;    /* Reset boundary conditions with fixed up points */ // added omp calls

	/* Repeat and rinse for the full time (aka corrector) step:  */
	fprintf(stderr,"f") ;
	#pragma omp parallel for \
		default(shared) \
		private(i,j,k)
	ZLOOP PLOOP psave[i][j][k] = p[i][j][k] ;
	ndt = advance(p, ph, dt,    p) ;

	fixup(p) ;  // added omp calls
    bound_prim(p) ; // added omp calls
	fixup_utoprim(p); // added omp calls
	bound_prim(p) ; // added omp calls


	/* Determine next time increment based on current characteristic speeds: */
        if(dt < 1.e-9) {
                fprintf(stderr,"timestep too small\n") ;
                exit(11) ;
        }

        /* increment time */
        t += dt ;

        /* set next timestep */
        if(ndt > SAFE*dt) ndt = SAFE*dt ;
        dt = ndt ;
        if(t + dt > tf) dt = tf - t ;  /* but don't step beyond end of run */

        /* done! */
}

/***********************************************************************************************/
/***********************************************************************************************
  advance():
  ---------
     -- responsible for what happens during a time step update, including the flux calculation, 
         the constrained transport calculation (aka flux_ct()), the finite difference 
         form of the time integral, and the calculation of the primitive variables from the 
         update conserved variables;
     -- also handles the "fix_flux()" call that sets the boundary condition on the fluxes;

***********************************************************************************************/
double advance(
	double pi[][N2+4][NPR], 
	double pb[][N2+4][NPR], 
	double Dt,
	double pf[][N2+4][NPR]
	)
{
  int i,j,k;
  double ndt,ndt1,ndt2,U[NPR],dU[NPR] ;
  struct of_geom geom ;
  struct of_state q ;
  
  double Uo[NPR],po[NPR] ;

  #pragma omp parallel for \
    default(shared) \
    private(i,j,k)
  ZLOOP PLOOP pf[i][j][k] = pi[i][j][k] ;        /* needed for Utoprim */
  
  fprintf(stderr,"0") ;
  ndt1 = fluxcalc(pb, F1, 1) ;
  ndt2 = fluxcalc(pb, F2, 2) ;

  fix_flux(F1,F2) ;  // added omp calls

  flux_ct(F1,F2) ;  // added omp calls

  /* evaluate diagnostics based on fluxes */
  diag_flux(F1,F2) ; 

  fprintf(stderr,"1") ;
  /** now update pi to pf **/
  #pragma omp parallel for \
    default(shared) \
    private(i,j,k,geom,q,U,dU)
  ZLOOP {

    get_geometry(i,j,CENT,&geom) ;

    source(pb[i][j],&geom,i,j,dU,Dt) ;

    get_state(pi[i][j],&geom,&q) ;
    primtoU(pi[i][j],&q,&geom,U) ;

    PLOOP {
      U[k] += Dt*(
		  - (F1[i+1][j][k] - F1[i][j][k])/dx[1]
		  - (F2[i][j+1][k] - F2[i][j][k])/dx[2]
		  + dU[k]
		  ) ;
    }

    pflag[i][j] = Utoprim_2d(U, geom.gcov, geom.gcon, geom.g, pf[i][j]);
    if( pflag[i][j] ) failimage[0][i+j*N1]++ ;
#if( DO_FONT_FIX ) 
    if( pflag[i][j] ) {
      pflag[i][j] = Utoprim_1dvsq2fix1(U, geom.gcov, geom.gcon, geom.g, pf[i][j], Katm[i] );
      if( pflag[i][j] ) { 
	failimage[1][i+j*N1]++ ;
	pflag[i][j] = Utoprim_1dfix1(U, geom.gcov, geom.gcon, geom.g, pf[i][j], Katm[i] );
	if( pflag[i][j] ) failimage[2][i+j*N1]++ ;
      }
    }
#endif

  }

  ndt = defcon * 1./(1./ndt1 + 1./ndt2) ;
  fprintf(stderr,"2") ;

  return(ndt) ;
}


/***********************************************************************************************/
/***********************************************************************************************
  fluxcalc():
  ---------
     -- sets the numerical fluxes, avaluated at the cell boundaries using the slope limiter
        slope_lim();

     -- only has HLL and Lax-Friedrichs  approximate Riemann solvers implemented;
        
***********************************************************************************************/
double fluxcalc(
	double pr[][N2+4][NPR], 
	double F[][N2+4][NPR], 
	int dir 
	)
{
	int i,j,k,idel,jdel,face ;
	double p_l[NPR],p_r[NPR],F_l[NPR],F_r[NPR],U_l[NPR],U_r[NPR] ;
	double cmax_l,cmax_r,cmin_l,cmin_r,cmax,cmin,ndt,dtij ;
	double ctop ;
	struct of_geom geom ;
	struct of_state state_l,state_r ;
	void rescale(double *pr, int which, int dir, int ii, int jj, int face, struct of_geom *geom) ;
	double bsq ;

        if     (dir == 1) {idel = 1; jdel = 0; face = FACE1;}
	else if(dir == 2) {idel = 0; jdel = 1; face = FACE2;}
	else { exit(10); }

	Timer time,tOMP;
	double tot_time,tot_timeOMP;
	int tid, threadstart, threadfinish, rows_per_thread;
	double ndt_private;


#if(TIMING_COMPARISON)
	//TESTING
	int size_N1p4_N2p4_NPR_d = (N1+4)*(N2+4)*NPR*sizeof(double);
	double * test_F = (double*) malloc(size_N1p4_N2p4_NPR_d);
	memcpy(test_F, &F[-2][-2][0], size_N1p4_N2p4_NPR_d);
	double * test_pr = (double*) malloc(size_N1p4_N2p4_NPR_d);
	memcpy(test_pr, &pr[-2][-2][0], size_N1p4_N2p4_NPR_d);
	double * test_dq = (double*) malloc(size_N1p4_N2p4_NPR_d);
	memcpy(test_dq, &dq[-2][-2][0], size_N1p4_N2p4_NPR_d);
	double test_ndt;
	startTime(&tOMP);
#endif //TIMING_COMPARISON

	ndt = 1.e9;

#if(RESCALE)
	#pragma omp parallel for \
		default(shared) \
		private(i,j,geom)
	ZSLOOP(-2,N1+1,-2,N2+1) {
		get_geometry(i,j,CENT,&geom);
		rescale(pr[i][j],FORWARD, dir, i,j,CENT,&geom);
	}
#endif //RESCALE

	/* then evaluate slopes */
	#pragma omp parallel for \
		default(shared) \
		private(i,j,k)
	ZSLOOP(-1,N1,-1,N2) {
		PLOOP {
			//-new get_geometry(i,j,CENT,&geom) ;
			//-new bsq = bsq_calc(pr[i][j],&geom) ;
			//-new if(bsq/pr[i][j][RHO] > 10. ||
			//-new    bsq/pr[i][j][UU]  > 1.e3) lim = MINM ;
			//-new else lim = MC ;
			dq[i][j][k] = slope_lim(
				pr[i-idel][j-jdel][k],
				pr[i][j][k],
				pr[i+idel][j+jdel][k]
				) ;
		}
	}
/*
	#pragma omp parallel for \
		default(shared) \
		private(i,j,k, \
				cmax_l,cmax_r,cmin_l,cmin_r,cmax,cmin,ctop,dtij, \
				geom,state_l,state_r, \
				p_l,p_r,F_l,F_r,U_l,U_r)
	ZSLOOP(-jdel,N1,-idel,N2) {

		// this avoids problems on the pole 
		if(dir == 2 && (j == 0 || j == N2)) {
			PLOOP F[i][j][k] = 0. ;
		} else {

            PLOOP {
                    p_l[k] = pr[i-idel][j-jdel][k] 
				+ 0.5*dq[i-idel][j-jdel][k] ;
                    p_r[k] = pr[i][j][k]   
				- 0.5*dq[i][j][k] ;
            }

			get_geometry(i,j,face,&geom) ;

#if(RESCALE)
			rescale(p_l,REVERSE,dir,i,j,face,&geom) ;
			rescale(p_r,REVERSE,dir,i,j,face,&geom) ;
#endif //RESCALE
			get_state(p_l,&geom,&state_l) ;
			get_state(p_r,&geom,&state_r) ;
		
			primtoflux(p_l,&state_l,dir,&geom,F_l) ;
			primtoflux(p_r,&state_r,dir,&geom,F_r) ;

			primtoflux(p_l,&state_l,TT, &geom,U_l) ;
			primtoflux(p_r,&state_r,TT, &geom,U_r) ;

			vchar(p_l,&state_l,&geom,dir,&cmax_l,&cmin_l) ;
			vchar(p_r,&state_r,&geom,dir,&cmax_r,&cmin_r) ;

			cmax = fabs(MY_MAX(MY_MAX(0., cmax_l),  cmax_r)) ;
			cmin = fabs(MY_MAX(MY_MAX(0.,-cmin_l), -cmin_r)) ;
			ctop = MY_MAX(cmax,cmin) ;

			PLOOP F[i][j][k] = 
				HLLF*(
				(cmax*F_l[k] + cmin*F_r[k] 
					- cmax*cmin*(U_r[k] - U_l[k]))/
						(cmax + cmin + SMALL) 
				) +
				LAXF*(
				0.5*(F_l[k] + F_r[k] 
					- ctop*(U_r[k] - U_l[k])) 
				) ;

		            // evaluate restriction on timestep 
		            cmax = MY_MAX(cmax,cmin) ;
		            dtij = cour*dx[dir]/cmax ;
			#pragma omp critical
			{ //TODO why isn't omp min reduction working?, make each thread keep local min and merge at end?
				if(dtij < ndt) ndt = dtij ;
			}
		}
	}
*/

	rows_per_thread = ceil(((float)N1+1-(-jdel))/numthreads);
	#pragma omp parallel \
		private(i,j,k,tid,threadstart,threadfinish,cmax_l,cmax_r,cmin_l,cmin_r,cmax,cmin,ctop,dtij,ndt_private, \
				geom, state_l, state_r, p_l, p_r, F_l, F_r, U_l, U_r )
	{
    tid = omp_get_thread_num();
	threadstart = -jdel+tid*rows_per_thread;
	threadfinish = (threadstart+rows_per_thread < N1+1 ? threadstart+rows_per_thread : N1+1);
	ndt_private = ndt;
	//ZSLOOP(-jdel,N1,-idel,N2) {
	//for(i=-jdel;i<=N1;i++) {
	//  for(j=-idel;j<=N2;j++) {
	for(i=threadstart;i<threadfinish;i++) {
	  for(j=-idel;j<=N2;j++) {

		/* this avoids problems on the pole */
		if(dir == 2 && (j == 0 || j == N2)) {
        	PLOOP F[i][j][k] = 0. ;
		}
		else {

            PLOOP {
                p_l[k] = pr[i-idel][j-jdel][k] 
					+ 0.5*dq[i-idel][j-jdel][k] ;
                p_r[k] = pr[i][j][k]   
					- 0.5*dq[i][j][k] ;
            }

			get_geometry(i,j,face,&geom) ;

#if(RESCALE)
			rescale(p_l,REVERSE,dir,i,j,face,&geom) ;
			rescale(p_r,REVERSE,dir,i,j,face,&geom) ;
#endif //RESCALE
			get_state(p_l,&geom,&state_l) ;
			get_state(p_r,&geom,&state_r) ;
		
			primtoflux(p_l,&state_l,dir,&geom,F_l) ;
			primtoflux(p_r,&state_r,dir,&geom,F_r) ;

			primtoflux(p_l,&state_l,TT, &geom,U_l) ;
			primtoflux(p_r,&state_r,TT, &geom,U_r) ;

			vchar(p_l,&state_l,&geom,dir,&cmax_l,&cmin_l) ;
			vchar(p_r,&state_r,&geom,dir,&cmax_r,&cmin_r) ;

			cmax = fabs(MY_MAX(MY_MAX(0., cmax_l),  cmax_r)) ;
			cmin = fabs(MY_MAX(MY_MAX(0.,-cmin_l), -cmin_r)) ;
			ctop = MY_MAX(cmax,cmin) ;


			PLOOP F[i][j][k] = 
				HLLF*(
						(cmax*F_l[k] + cmin*F_r[k] 
						- cmax*cmin*(U_r[k] - U_l[k]))/
						(cmax + cmin + SMALL) 
					 ) +
				LAXF*(
						0.5*(F_l[k] + F_r[k] 
						- ctop*(U_r[k] - U_l[k])) 
					 ) ;

		            /* evaluate restriction on timestep */
		    cmax = MY_MAX(cmax,cmin) ;
		    dtij = cour*dx[dir]/cmax ;
			if(dtij < ndt_private) ndt_private = dtij ;
		}
	  }
	}
	#pragma omp critical
	{
	if(ndt_private < ndt) ndt = ndt_private ;
	}// end omp critical

	} //end omp parallel




#if(RESCALE)
	#pragma omp parallel for \
		default(shared) \
		private(i,j,geom)
	ZSLOOP(-2,N1+1,-2,N2+1) {
		get_geometry(i,j,CENT,&geom);
		rescale(pr[i][j],REVERSE, dir, i,j,CENT,&geom);
	}
#endif //RESCALE

#if(TIMING_COMPARISON)
	stopTime(&tOMP);
#endif //TIMING_COMPARISON



// NON-OMP VERSION //

#if(TIMING_COMPARISON)
	startTime(&time);
#if(RESCALE)
	/** evaluate slopes of primitive variables **/
	/* first rescale */
	ZSLOOP(-2,N1+1,-2,N2+1) {
		get_geometry(i,j,CENT,&geom) ;
		rescale(&test_pr[((i+2)*(N2+4)+(j+2))*NPR],FORWARD, dir, i,j,CENT,&geom) ;
	}
#endif //RESCALE

	ZSLOOP(-1,N1,-1,N2) PLOOP {
   	        //-new get_geometry(i,j,CENT,&geom) ;
		//-new bsq = bsq_calc(pr[i][j],&geom) ;
		//-new if(bsq/pr[i][j][RHO] > 10. ||
		//-new    bsq/pr[i][j][UU]  > 1.e3) lim = MINM ;
		//-new else lim = MC ;
		test_dq[((i+2)*(N2+4)+(j+2))*NPR+k] = slope_lim(
				test_pr[((i-idel+2)*(N2+4)+(j-jdel+2))*NPR+k],
				test_pr[((i+2)*(N2+4)+(j+2))*NPR+k],
				test_pr[((i+idel+2)*(N2+4)+(j+jdel+2))*NPR+k]
				) ;
	}

	test_ndt = 1.e9 ;
	ZSLOOP(-jdel,N1,-idel,N2) {

		/* this avoids problems on the pole */
		if(dir == 2 && (j == 0 || j == N2)) {
                        //PLOOP test_F[i][j][k] = 0. ;
                        PLOOP test_F[((i+2)*(N2+4)+(j+2))*NPR+k] = 0. ;
		}
		else {

                PLOOP { 
                        p_l[k] = test_pr[((i-idel+2)*(N2+4)+(j-jdel+2))*NPR+k] 
					+ 0.5*test_dq[((i-idel+2)*(N2+4)+(j-jdel+2))*NPR+k] ;
                        p_r[k] = test_pr[((i+2)*(N2+4)+(j+2))*NPR+k]   
					- 0.5*test_dq[((i+2)*(N2+4)+(j+2))*NPR+k] ;
                }

		get_geometry(i,j,face,&geom) ;

#if(RESCALE)
		rescale(p_l,REVERSE,dir,i,j,face,&geom) ;
		rescale(p_r,REVERSE,dir,i,j,face,&geom) ;
#endif //RESCALE
		get_state(p_l,&geom,&state_l) ;
		get_state(p_r,&geom,&state_r) ;
		
		primtoflux(p_l,&state_l,dir,&geom,F_l) ;
		primtoflux(p_r,&state_r,dir,&geom,F_r) ;

		primtoflux(p_l,&state_l,TT, &geom,U_l) ;
		primtoflux(p_r,&state_r,TT, &geom,U_r) ;

		vchar(p_l,&state_l,&geom,dir,&cmax_l,&cmin_l) ;
		vchar(p_r,&state_r,&geom,dir,&cmax_r,&cmin_r) ;

		cmax = fabs(MY_MAX(MY_MAX(0., cmax_l),  cmax_r)) ;
		cmin = fabs(MY_MAX(MY_MAX(0.,-cmin_l), -cmin_r)) ;
		ctop = MY_MAX(cmax,cmin) ;


		//PLOOP test_F[i][j][k] = 
		PLOOP test_F[((i+2)*(N2+4)+(j+2))*NPR+k] = 
			HLLF*(
			(cmax*F_l[k] + cmin*F_r[k] 
				- cmax*cmin*(U_r[k] - U_l[k]))/
					(cmax + cmin + SMALL) 
			) +
			LAXF*(
			0.5*(F_l[k] + F_r[k] 
				- ctop*(U_r[k] - U_l[k])) 
			) ;

                /* evaluate restriction on timestep */
                cmax = MY_MAX(cmax,cmin) ;
                dtij = cour*dx[dir]/cmax ;
		if(dtij < test_ndt) test_ndt = dtij ;

		}
	}

#if(RESCALE)

	/** evaluate slopes of primitive variables **/
	/* first rescale */
	ZSLOOP(-2,N1+1,-2,N2+1) {
		get_geometry(i,j,CENT,&geom) ;
		rescale(&test_pr[((i+2)*(N2+4)+(j+2))*NPR],REVERSE, dir, i,j,CENT,&geom) ;
	}

#endif //RESCALE

	stopTime(&time);
	tot_time = elapsedTime(time);
	tot_timeOMP = elapsedTime(tOMP);
	tot_time_flux += tot_time;
	tot_time_flux_omp += tot_timeOMP;

	//printf("TESTING pr... \n"); fflush(stdout);
	verify_d(&pr[-2][-2][0],test_pr,(N1+4)*(N2+4)*NPR); // for testing dq,pr,F 
	//printf("FINISHED TESING pr... \n"); fflush(stdout);
	//printf("TESTING dq... \n"); fflush(stdout);
	verify_d(&dq[-2][-2][0],test_dq,(N1+4)*(N2+4)*NPR); // for testing dq,pr,F 
	//printf("FINISHED TESING dq... \n"); fflush(stdout);
	//printf("TESTING F... \n"); fflush(stdout);
	verify_d(&F[-2][-2][0],test_F,(N1+4)*(N2+4)*NPR); // for testing dq,pr,F 
	//printf("FINISHED TESING F... \n"); fflush(stdout);
	//printf("TESTING ndt... \n"); fflush(stdout);
	verify_d(&ndt,&test_ndt,1); // for testing dq,pr,F 
	//printf("FINISHED TESING ndt... \n"); fflush(stdout);
	free(test_F);
	free(test_dq);
	free(test_pr);


#endif //TIMING_COMPARISON

	return(ndt) ;
}


/***********************************************************************************************/
/***********************************************************************************************
  flux_ct():
  ---------
     -- performs the flux-averaging used to preserve the del.B = 0 constraint (see Toth 2000);
        
***********************************************************************************************/
void flux_ct(double F1[][N2+4][NPR], double F2[][N2+4][NPR])
{
	int i,j ;
	static double emf[N1+1][N2+1] ;

	/* calculate EMFs */
	/* Toth approach: just average */
	#pragma omp parallel for \
		default(shared) \
		private(i,j)
	ZSLOOP(0,N1,0,N2) emf[i][j] = 0.25*(F1[i][j][B2] + F1[i][j-1][B2]
					  - F2[i][j][B1] - F2[i-1][j][B1]) ;

	/* rewrite EMFs as fluxes, after Toth */
	#pragma omp parallel for \
		default(shared) \
		private(i,j)
    ZSLOOP(0,N1,0,N2-1) {
            F1[i][j][B1] = 0. ;
            F1[i][j][B2] =  0.5*(emf[i][j] + emf[i][j+1]) ;
    }

	#pragma omp parallel for \
		default(shared) \
		private(i,j)
    ZSLOOP(0,N1-1,0,N2) {
            F2[i][j][B1] = -0.5*(emf[i][j] + emf[i+1][j]) ;
            F2[i][j][B2] = 0. ;
	}

}


void init_mpi_omp(int argc, char *argv[])
{
	//MPI_Init(&argc, &argv);
	//MPI_Comm_size(MPI_COMM_WORLD,&numranks);

	if (argc > 1) {
		numthreads = atoi(argv[1]);
	} else {
		numthreads = MAXTHREADS;
	}

	printf("Threads requested: %d\n", numthreads);
	omp_set_num_threads(numthreads);


#if(TIMING_COMPARISON)
	tot_time_flux_dq = 0;
	tot_time_flux_dq_omp = 0;
	tot_time_flux_F = 0;
	tot_time_flux_F_omp = 0;
	tot_time_flux_pr = 0;
	tot_time_flux_pr_omp = 0;
	tot_time_flux = 0;
	tot_time_flux_omp = 0;
#endif
}

void finalize_mpi_omp()
{

#if(TIMING_COMPARISON)
	//printf("Total time Fluxcalc dq loop (orig,omp): (%f,%f), speedup: %f\n", tot_time_flux_dq, tot_time_flux_dq_omp, tot_time_flux_dq/tot_time_flux_dq_omp);
	//printf("Total time Fluxcalc F loop (orig,omp): (%f,%f), speedup: %f\n", tot_time_flux_F, tot_time_flux_F_omp, tot_time_flux_F/tot_time_flux_F_omp);
	//printf("Total time Fluxcalc rescale1 loop (orig,omp): (%f,%f), speedup: %f\n", tot_time_flux_pr, tot_time_flux_pr_omp, tot_time_flux_pr/tot_time_flux_pr_omp);
	printf("Total time Fluxcalc (orig,omp): (%f,%f), speedup: %f\n", tot_time_flux, tot_time_flux_omp, tot_time_flux/tot_time_flux_omp);
#endif
}





