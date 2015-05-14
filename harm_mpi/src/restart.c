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


/* restart functions; restart_init and restart_dump */

#include "decs.h"



/***********************************************************************/
/***********************************************************************
  restart_write():
     -- writes current state of primitive variables to the
        checkpointing/restart file.
     -- uses ASCII text format ;
     -- when changing this routine, be sure to make analogous changes
        in restart_read();
************************************************************************/
void restart_write()
{
  FILE *fp;
  char filename[100];
  int idum, i, j, k;

  MPI_File file;
  MPI_Status status;

  /*************************************************************
    Rank 0 writes the header of the restart file:
  *************************************************************/
  if (WorldRank == 0) {
    if (rdump_cnt%2 == 0) {
      fp = fopen("dumps/rdump0","wt") ;
      sprintf(filename, "dumps/rdump0");
      fprintf(stderr, "RESTART  file=dumps/rdump0\n");
    }
    else {
      fp = fopen("dumps/rdump1","wt") ;
      sprintf(filename, "dumps/rdump1");
      fprintf(stderr, "RESTART file=dumps/rdump1\n");
    }

    if (fp == NULL) {
      fprintf(stderr, "Cannot open restart file\n");
      sprintf(filename, "");
    } else {
      fprintf(fp, FMT_INT_OUT, GlobalN1 );
      fprintf(fp, FMT_INT_OUT, GlobalN2 );
      fprintf(fp, FMT_DBL_OUT, t        );
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

      fprintf(fp, "\n");
      fclose(fp) ;
    }
  }

  MPI_Bcast(filename, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
  if (!strcmp(filename, "")) {
      exit(2);
  }

  /* convert data to a string */
  char *data_as_text = (char *) malloc(
      CHARSPERNUM * NPR * (N1 + 4) * (N2 + 4) * sizeof(char));
  int count = 0;
  ZSLOOP (-2, N1+1, -2, N2+1) {
    for (k = 0; k < NPR - 1; k++) {
      sprintf(&data_as_text[count*CHARSPERNUM], FMT_FILE, p[i][j][k]);
      count++;
    }

    sprintf(&data_as_text[count*CHARSPERNUM], ENDFMT_FILE, p[i][j][NPR-1]);
    count++;
  }

  /* open the file, set the view, write data and close the file */
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY, MPI_INFO_NULL,
      &file);
  MPI_File_set_view(file, RST_HEADER_DISP, array_as_string, local_array,
      "native", MPI_INFO_NULL);
  MPI_File_write_all(file, data_as_text, (N1 + 4) * (N2 + 4), array_as_string,
      &status);
  MPI_File_close(&file);

  /* clean up */
  free(data_as_text);
  rdump_cnt++;
  return;
}

/***********************************************************************/
/***********************************************************************
  restart_init():
     -- main driver for setting initial conditions from a checkpoint
        or restart file.
     -- determines if there are any restart files to use and then
         lets the  user choose if there are more than one file.
     -- then calls initializes run with restart data;
************************************************************************/
int restart_init()
{
  FILE *fp0, *fp1;
  char filename[100];
  char ans[100];
  int i, j, k;

  /* set up global arrays */
  set_arrays();

  /********************************************************************
   Rank 0 checks to see which restart files exist.
   Use the only one that exists, else prompt user to decide
     which one to use if we have a choice : 
  ********************************************************************/
  if (WorldRank == 0) {
    fp0 = fopen("dumps/rdump0","r") ;
    fp1 = fopen("dumps/rdump1","r") ;

    if ( (fp0 == NULL) && (fp1 == NULL) ) {
      fprintf(stderr,"No restart file\n") ;
      sprintf(filename, "");
    }
    else {
      fprintf(stderr,"\nRestart file exists! \n") ;
      if ( fp0 == NULL ) {
        fprintf(stderr,"Using dumps/rdump1 ... \n");
        sprintf(filename, "dumps/rdump1");
      }
      else if ( fp1 == NULL ) {
        fprintf(stderr,"Using dumps/rdump0 ... \n");
        sprintf(filename, "dumps/rdump0");
      }
      else {
        fprintf(stderr,"Use dumps/rdump0 (0) or dumps/rdump1 (1)?   [0|1]  \n");
        fscanf(stdin,"%s",ans) ;
        if (strncmp(ans,"0",1) == 0) {
          sprintf(filename, "dumps/rdump0");
        }
        else {
          sprintf(filename, "dumps/rdump1");
        }
      }
    }
  }

  MPI_Bcast(filename, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
  if (!strcmp(filename, "")) {
      return(0) ;
  }

  /********************************************************************
   Now that we know we are restarting from a checkpoint file, then
   we need to read in data, assign grid functions and define the grid:
  ********************************************************************/
  restart_read(filename);

  /* set half-step primitives everywhere */
  ZSLOOP(-2,N1+1,-2,N2+1) PLOOP ph[i][j][k] = p[i][j][k] ;

  /* set metric functions */
  set_grid() ;

  /* bound */
  bound_prim(p) ;
  bound_prim(ph) ;


#if( DO_FONT_FIX ) 
  set_Katm();
#endif 

  /***********************************************************************
    Make any changes to parameters in restart file  here:
      e.g., cour = 0.4 , change in limiter...
  ************************************************************************/
  //lim = MC ;
  //cour = 0.9 ;
  //lim = VANL ;
  //tf = 4000. ;

  if (WorldRank == 0)
    fprintf(stderr,"done with restart init.\n") ;

  /* done! */
  return(1) ;
}

/***********************************************************************/
/***********************************************************************
  restart_read():
     -- reads in data from the restart file, which is specified in
         restart_init() but is usually named "dumps/rdump[0,1]"
************************************************************************/
void restart_read(char *filename)
{
  FILE *fp;
  int idum, i, j, k ;

  MPI_File file;
  MPI_Status status;

  /*************************************************************
	  Everyone can simultaneously read the restart file:
  *************************************************************/
  fp = fopen(filename, "r");

  fscanf(fp, "%d", &idum);
  if (idum != GlobalN1) {
    fprintf(stderr, "error reading restart file; N1 differs\n");
    exit(3);
  }
  fscanf(fp, "%d", &idum);
  if (idum != GlobalN2) {
    fprintf(stderr, "error reading restart file; N2 differs\n");
    exit(4);
  }

  fscanf(fp, "%lf", &t        );
  fscanf(fp, "%lf", &tf       );
  fscanf(fp, "%d",  &nstep    );
  fscanf(fp, "%lf", &a        );
  fscanf(fp, "%lf", &gam      );
  fscanf(fp, "%lf", &cour     );
  fscanf(fp, "%lf", &DTd      );
  fscanf(fp, "%lf", &DTl      );
  fscanf(fp, "%lf", &DTi      );
  fscanf(fp, "%d",  &DTr      );
  fscanf(fp, "%d",  &dump_cnt );
  fscanf(fp, "%d",  &image_cnt);
  fscanf(fp, "%d",  &rdump_cnt);
  fscanf(fp, "%lf", &dt       );
  fscanf(fp, "%d",  &lim      );
  fscanf(fp, "%d",  &failed   );
  fscanf(fp, "%lf", &Rin      );
  fscanf(fp, "%lf", &Rout     );
  fscanf(fp, "%lf", &hslope   );
  fscanf(fp, "%lf", &R0       );

  fclose(fp);

  /* open the file, set the view, read the contents and close the file */
  char *data_as_text = (char *) malloc(
      CHARSPERNUM * NPR * (N1 + 4) * (N2 + 4) * sizeof(char));
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL,
      &file);
  MPI_File_set_view(file, RST_HEADER_DISP, array_as_string, local_array,
      "native", MPI_INFO_NULL);
  MPI_File_read_all(file, data_as_text, (N1 + 4) * (N2 + 4), array_as_string,
      &status);
  MPI_File_close(&file);

  /* convert strings to doubles */
  int count = 0;
  ZSLOOP(-2, N1+1, -2, N2+1) {
    PLOOP {
      sscanf(&data_as_text[count*CHARSPERNUM], "%lf", &(p[i][j][k]));
      count++;
    }
  }

  free(data_as_text);
  return ;
}

#undef FMT_DBL_OUT
#undef FMT_INT_OUT
