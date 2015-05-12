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

void dump(char *filename)
{
	int i,j,k ;
	double divb ;
	double X[NDIM] ;
	double r,th,vmin,vmax ;
	struct of_geom geom ;
	struct of_state q ;
	FILE *fp;

	MPI_File file;
	MPI_Status status;
	MPI_Datatype data_as_string;
	MPI_Datatype dump_array;

	/***************************************************************
	  Rank 0 writes header information :
	***************************************************************/

	if (WorldRank == 0) {
	  fp = fopen(filename, "w");

    fprintf(fp, FMT_DBL_OUT, t        );
    fprintf(fp, FMT_INT_OUT, N1       );
    fprintf(fp, FMT_INT_OUT, N2       );
    fprintf(fp, FMT_DBL_OUT, startx[1]);
    fprintf(fp, FMT_DBL_OUT, startx[2]);
    fprintf(fp, FMT_DBL_OUT, dx[1]    );
    fprintf(fp, FMT_DBL_OUT, dx[2]    );
    fprintf(fp, FMT_DBL_OUT, tf       );
    fprintf(fp, FMT_INT_OUT, nstep    );
    fprintf(fp, FMT_DBL_OUT, a        );
    fprintf(fp, FMT_DBL_OUT, gam      );
    fprintf(fp, FMT_DBL_OUT, cour     );
    fprintf(fp, FMT_DBL_OUT, DTd      );
    fprintf(fp, FMT_DBL_OUT, DTl      );
    fprintf(fp, FMT_DBL_OUT, DTi      );
    fprintf(fp, FMT_INT_OUT, DTr      );
    fprintf(fp, FMT_INT_OUT, dump_cnt );
    fprintf(fp, FMT_INT_OUT, image_cnt);
    fprintf(fp, FMT_INT_OUT, rdump_cnt);
    fprintf(fp, FMT_DBL_OUT, dt       );
    fprintf(fp, FMT_INT_OUT, lim      );
    fprintf(fp, FMT_INT_OUT, failed   );
    fprintf(fp, FMT_DBL_OUT, Rin      );
    fprintf(fp, FMT_DBL_OUT, Rout     );
    fprintf(fp, FMT_DBL_OUT, hslope   );
    fprintf(fp, FMT_DBL_OUT, R0       );

    fprintf(fp,"\n") ;
    fclose(fp);
	}
		
	/***************************************************************
	  Write simulation information :
	***************************************************************/

	int num;
	if (!failed)
	  num = 34;
	else
	  num = 13;

	int count = 0;
	char *data_as_text = (char *) malloc(CHARSPERNUM*num*N1*N2*sizeof(char));
	ZSLOOP(0, N1-1, 0, N2-1) {
		coord(i, j, CENT, X);
		bl_coord(X, &r, &th);

		sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, X[1]);
		sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, X[2]);
		sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, r);
		sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, th);
		PLOOP sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, p[i][j][k]);

//		fprintf(fp, FMT_DBL_OUT, X[1]       );
//		fprintf(fp, FMT_DBL_OUT, X[2]       );
//		fprintf(fp, FMT_DBL_OUT, r          );
//		fprintf(fp, FMT_DBL_OUT, th         );
//		PLOOP fprintf(fp, FMT_DBL_OUT, p[i][j][k] );
		
                /* divb flux-ct defn; corner-centered.  Use
		   only interior corners */
		if(((RowRank*N1 + i) > 0) && ((ColRank*N2 + j) > 0) &&
		   ((RowRank*N1 + i) < GlobalN1) && ((ColRank*N2 + j) < GlobalN2)) {
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
		}
		else {
		  divb = 0. ;
		}

		if (!failed)
		  sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, divb);
		else
		  sprintf(&data_as_text[count++*CHARSPERNUM], ENDFMT_RST, divb);

//		fprintf(fp, FMT_DBL_OUT, divb     );

		if(!failed) {
			get_geometry(i, j, CENT, &geom) ;
			get_state(p[i][j], &geom, &q) ;

			for(k=0;k<NDIM;k++) sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, q.ucon[k]);
			for(k=0;k<NDIM;k++) sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, q.ucov[k]);
			for(k=0;k<NDIM;k++) sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, q.bcon[k]);
			for(k=0;k<NDIM;k++) sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, q.bcov[k]);
//			for(k=0;k<NDIM;k++) fprintf(fp,FMT_DBL_OUT,q.ucon[k]) ;
//			for(k=0;k<NDIM;k++) fprintf(fp,FMT_DBL_OUT,q.ucov[k]) ;
//			for(k=0;k<NDIM;k++) fprintf(fp,FMT_DBL_OUT,q.bcon[k]) ;
//			for(k=0;k<NDIM;k++) fprintf(fp,FMT_DBL_OUT,q.bcov[k]) ;

			vchar(p[i][j], &q, &geom, 1, &vmax, &vmin);
			sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, vmin);
			sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, vmax);
//			fprintf(fp, FMT_DBL_OUT, vmin );
//			fprintf(fp, FMT_DBL_OUT, vmax );

			vchar(p[i][j], &q, &geom, 2, &vmax, &vmin);
			sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, vmin);
			sprintf(&data_as_text[count++*CHARSPERNUM], FMT_RST, vmax);
//			fprintf(fp, FMT_DBL_OUT, vmin );
//			fprintf(fp, FMT_DBL_OUT, vmax );

			sprintf(&data_as_text[count++*CHARSPERNUM], ENDFMT_RST, geom.g);
//			fprintf(fp, FMT_DBL_OUT, geom.g );
		}

//		fprintf(fp,"\n") ;
	}

	/*
	 * NOTE: The following type creation has to be done every time because we
	 * don't know what the value of 'num' will be. In most cases, it will be 34,
	 * but when a failure occurs, it'll be 13. Since there is no way of knowing
	 * this information during the initialization phase, we must, unfortunately,
	 * perform this tedious calculation every single time. Hopefully, it doesn't
	 * kill performance by too much.
	 */

	/* need a custom type for printing a whole row in one go */
	MPI_Type_contiguous(num*CHARSPERNUM, MPI_CHAR, &data_as_string);
	MPI_Type_commit(&data_as_string);

	/* need a new type to set file view since we're looking at N1xN2 elements */
	int globalsizes[2] = {GlobalN1, GlobalN2};
  int localsizes[2]  = {N1, N2};
  int starts[2]      = {N1*RowRank, N2*ColRank};
  int order          = MPI_ORDER_C;
  MPI_Type_create_subarray(2, globalsizes, localsizes, starts, order,
                           data_as_string, &dump_array);
  MPI_Type_commit(&dump_array);

	/* open the file, set the view, write data and close the file */
  MPI_File_open(MPI_COMM_WORLD, filename,
                MPI_MODE_APPEND,
                MPI_INFO_NULL, &file);

  MPI_File_set_view(file, 0, MPI_CHAR, dump_array, "native", MPI_INFO_NULL);
  MPI_File_write_all(file, data_as_text, N1*N2, data_as_string, &status);
  MPI_File_close(&file);

  /* clean up */
  free(data_as_text);
}
