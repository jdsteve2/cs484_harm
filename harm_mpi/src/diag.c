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

/* all diagnostics subroutine */

void diag(call_code)
int call_code ;
{
	char dfnam[100],ifnam[100] ;
	int i,j,k ;
	FILE *dump_file;
	double U[NPR],pp,e,rmed,divb,divbmax,e_fin,m_fin,gamma ;
	struct of_geom geom ;
	struct of_state q ;
	int imax,jmax ;
    double lred[3], gred[3];
    int tmp[2];
    struct {
        double d;
        int i;
    } ldi, gdi;

	static double e_init,m_init ;
	static FILE *ener_file ;

	if(call_code==INIT_OUT) {
		/* set things up */
        if(WorldRank == 0) {
            ener_file = fopen("ener.out","a") ;
            if(ener_file==NULL) {
                fprintf(stderr,"error opening energy output file\n") ;
                exit(1) ;
            }
        }
	}

	/* calculate conserved quantities */
	if(call_code==INIT_OUT || 
	   call_code==LOG_OUT ||
	   call_code==FINAL_OUT &&
	   !failed) {
		pp = 0. ;
		e = 0. ;
		rmed = 0. ;
		divbmax = 0. ;
		imax = 0 ; 
		jmax = 0 ;
		ZSLOOP(0,N1-1,0,N2-1) {
			get_geometry(i,j,CENT,&geom) ;
			get_state(p[i][j],&geom,&q) ;
			primtoU(p[i][j],&q,&geom,U) ;

			rmed += U[RHO]*dV ;
			pp += U[U3]*dV ;
			e += U[UU]*dV ;

			/* flux-ct defn */
			divb = fabs( 0.5*(
				p[i][j][B1]*gdet[i][j][CENT] 
				+ p[i][j-1][B1]*gdet[i][j-1][CENT]
				- p[i-1][j][B1]*gdet[i-1][j][CENT] 
				- p[i-1][j-1][B1]*gdet[i-1][j-1][CENT]
				)/dx[1] +
				0.5*(
				p[i][j][B2]*gdet[i][j][CENT] 
				+ p[i-1][j][B2]*gdet[i-1][j][CENT]
				- p[i][j-1][B2]*gdet[i][j-1][CENT] 
				- p[i-1][j-1][B2]*gdet[i-1][j-1][CENT]
				)/dx[2]) ;

			if(divb > divbmax && ((RowRank*N1 + i) > 0) && ((ColRank*N2 + j) > 0)) {
				imax = i ;
				jmax = j ;
				divbmax = divb ;
			}
		}

        // Reduce: imax, jmax, divbmax such that i and j correspond to divbmax (to zero)
        ldi.d = divbmax;
        ldi.i = WorldRank;
        MPI_Allreduce(&ldi, &gdi, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        if(gdi.i != 0) {
            if(WorldRank == 0) {
                MPI_Recv(tmp, 2, MPI_INT, gdi.i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                divbmax = gdi.d;
                imax = tmp[0];
                jmax = tmp[1];
            }
            if(WorldRank == gdi.i) {
                tmp[0] = RowRank*N1 + imax;
                tmp[1] = ColRank*N2 + jmax;
                MPI_Send(tmp, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
        }
        
        // Reduce: rmed, pp, e : sum (to zero)
        lred[0] = rmed; lred[1] = pp; lred[2] = e;
        MPI_Reduce(lred, gred, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        rmed = gred[0]; pp = gred[1]; e = gred[2];
	}

	if(call_code == INIT_OUT) {
		e_init = e ;
		m_init = rmed ;
	}

	if(call_code == FINAL_OUT) {
		e_fin = e ;
		m_fin = rmed ;
        if(WorldRank == 0) {
            fprintf(stderr,"\n\nEnergy: ini,fin,del: %g %g %g\n",
                    e_init,e_fin,(e_fin-e_init)/e_init) ;
            fprintf(stderr,"mass: ini,fin,del: %g %g %g\n",
                    m_init,m_fin,(m_fin-m_init)/m_init) ;
        }
	}

	if(call_code == INIT_OUT || 
	   call_code == LOG_OUT ||
	   call_code == FINAL_OUT) {
        if(WorldRank == 0) {
            fprintf(stderr,"LOG      t=%g \t divbmax: %d %d %g\n",
                    t,imax,jmax,divbmax) ;
            // MPI Get the last value from where required to print
            fprintf(ener_file,"%10.5g %10.5g %10.5g %10.5g %15.8g %15.8g ",
                    t,rmed,pp,e,p[N1/2][N2/2][UU]*pow(p[N1/2][N2/2][RHO],-gam),
                    p[N1/2][N2/2][UU]) ;
            fprintf(ener_file,"%15.8g %15.8g %15.8g ",mdot,edot,ldot) ;
            fprintf(ener_file,"\n") ;
            fflush(ener_file) ;
        }
	}


	/* dump at regular intervals */
	if(call_code == INIT_OUT || 
	   call_code == DUMP_OUT ||
	   call_code == FINAL_OUT) {

		/* make regular dump file */
    if (WorldRank == 0) {
        sprintf(dfnam, "dumps/dump%03d", dump_cnt) ;
        fprintf(stderr, "DUMP     file=%s\n", dfnam) ;
    }

		dump(dfnam);
		dump_cnt++;
	}
	
	/* image dump at regular intervals */
	if(call_code == IMAGE_OUT || 
	   call_code == INIT_OUT ||
	   call_code == FINAL_OUT) {

        // TODO: MPI Concurrent Write Image
//		image_all( image_cnt );

		image_cnt++ ;
	}
}


/** some diagnostic routines **/


void fail(int fail_type)
{

	failed = 1 ;

	fprintf(stderr,"\n\nfail: %d %d %d\n",icurr,jcurr,fail_type) ;

	area_map(icurr,jcurr,p) ;
	
	fprintf(stderr,"fail\n") ;

	diag(FINAL_OUT) ;

	/* for diagnostic purposes */
	exit(0) ;
}



/* map out region around failure point */
void area_map(int i, int j, double (** prim)[NPR])
{
	int k ;

	fprintf(stderr,"area map\n") ;

	PLOOP {
		fprintf(stderr,"variable %d \n",k) ;
		fprintf(stderr,"i = \t %12d %12d %12d\n",i-1,i,i+1) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j+1,
				prim[i-1][j+1][k],
				prim[i][j+1][k],
				prim[i+1][j+1][k]) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j,
				prim[i-1][j][k],
				prim[i][j][k],
				prim[i+1][j][k]) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j-1,
				prim[i-1][j-1][k],
				prim[i][j-1][k],
				prim[i+1][j-1][k]) ;
	}

	/* print out other diagnostics here */

}

/* evaluate fluxed based diagnostics; put results in
 * global variables */

void diag_flux(double (** F1)[NPR], double (** F2)[NPR])
{
	int j ;
    double lsum[3], gsum[3];

        mdot = edot = ldot = 0. ;
        if(RowRank == 0) {
            for(j=0;j<N2;j++) {
                mdot += F1[0][j][RHO]*2.*M_PI*dx[2] ;
                edot -= (F1[0][j][UU] - F1[0][j][RHO])*2.*M_PI*dx[2] ;
                ldot += F1[0][j][U3] *2.*M_PI*dx[2] ;
            }
            // Reduce mdot, edot, ldot : sum (To zero, only used for print above)
            // Only in CommRow
            lsum[0] = mdot;
            lsum[1] = edot;
            lsum[2] = ldot;
            MPI_Reduce(lsum, gsum, 3, MPI_DOUBLE, MPI_SUM, 0, CommRow);
            // FIXME: Assign gsum to mdot, edot and ldot, right?
        }
}
