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
#include <ctype.h>


/* Decide whether to use raw R8 (0) or PPM (1) image formats */
/* NOTE: PLEASE DON'T USE R8 FOR NOW AS IT HASN'T BEEN PORTED TO MPI */
#define MAKE_PPM_IMAGE (1)

#if( MAKE_PPM_IMAGE ) 
#define IMGEXT "ppm"
#else 
#define IMGEXT "r8"
#endif


/* Number of colors used in the images: */
#define NCOLORS (256)


/* Local Mnemonics */
#define IRHO  (0) 
#define IUU   (1)
#define IBSQ  (2)
#define IGAM  (3)

#define RED   (0)
#define GREEN (1)
#define BLUE  (2)


/* Global variables for PPM images */
static unsigned int color_map[3][NCOLORS];   
static int ppm_height, ppm_ncolors, ppm_width;


/* Local functions */
void get_color_map(void); 
void image(double *f, char *fname);
void image_r8(double *f, char *fname);
void image_ppm(double *f, char *fname);


/******************************************************************************/
/******************************************************************************
  image_all(): 
  -----------
   -- Main driver for generating "image" or snapshots of any desired quantity;
   -- This is responsible for calculating the outputted quantities and setting
       the names of the images
   -- This is the only routine in this file that a user should modify, all 
      other routines merely control how the images are created. 
   -- The image generating routines use a linear scale, so be sure to 
      take the log of any quantity you want to see in log color scale.

******************************************************************************/
void image_all( int image_count ) 
{ 

  int i,j,k, i_img;
  static int first_call = 1;
  static char ifnam[3*NIMG+1][100];
  double gamma;
  struct of_geom geom ;
  static const double fimage_logmin = 1.e-15;

  
  if( (IGAM+1) != NIMG ) {
    fprintf(stderr, "image_all(): Index problem with fimage[] \n");
    fflush(stderr);
    exit(1);
  }

#if( MAKE_PPM_IMAGE ) 
  if( first_call ) { 
    get_color_map();
    first_call = 0 ;
  }
#endif

  /************************************************************************
    Set the names of the image files to be generated now : 
  ************************************************************************/
  i_img = 0 ;
  sprintf(ifnam[i_img++], "images/im_rho_%04d.%s",image_count, IMGEXT) ;
  sprintf(ifnam[i_img++], "images/im_u_%04d.%s"  ,image_count, IMGEXT) ;
  sprintf(ifnam[i_img++], "images/im_bsq_%04d.%s",image_count, IMGEXT) ;
  sprintf(ifnam[i_img++], "images/im_gam_%04d.%s",image_count, IMGEXT) ;

  sprintf(ifnam[i_img++], "images/im_lrho_%04d.%s",image_count, IMGEXT) ;
  sprintf(ifnam[i_img++], "images/im_lu_%04d.%s"  ,image_count, IMGEXT) ;
  sprintf(ifnam[i_img++], "images/im_lbsq_%04d.%s",image_count, IMGEXT) ;
  sprintf(ifnam[i_img++], "images/im_lgam_%04d.%s",image_count, IMGEXT) ;

  sprintf(ifnam[i_img++], "images/failu2p1_%04d.%s",image_count, IMGEXT) ;
  sprintf(ifnam[i_img++], "images/failu2p2_%04d.%s",image_count, IMGEXT) ;
  sprintf(ifnam[i_img++], "images/failu2p3_%04d.%s",image_count, IMGEXT) ;
  sprintf(ifnam[i_img++], "images/failgamc_%04d.%s",image_count, IMGEXT) ;
  sprintf(ifnam[i_img++], "images/failfint_%04d.%s",image_count, IMGEXT) ;


  /************************************************************************
    Calculate the functions to be imaged : 
        -- the log versions overwrite the non-log versions;
        -- calculate only the non-log version here
  ************************************************************************/
  k = 0 ;
  IMAGELOOP {
    get_geometry(i,j,CENT,&geom) ;
    if( gamma_calc(p[i][j],&geom,&gamma) ) { gamma = 1.; }
    
    fimage[IRHO][k] = p[i][j][RHO] ; 
    fimage[IUU ][k] = p[i][j][UU ] ;
    fimage[IBSQ][k] = bsq_calc( p[i][j], &geom ) ; 
    fimage[IGAM][k] = gamma ;
    k++;
  }


  /************************************************************************
    Output non-log versions:
  ************************************************************************/
  for( i_img = 0 ; i_img < NIMG; i_img++ ) {  
    image( fimage[i_img], ifnam[i_img] );  
  }

  /************************************************************************
    Make log version of the image functions 
  ************************************************************************/
  for( i = 0 ; i < NIMG*N1*N2; i++ ) { 
    fimage[0][i] = log( fabs(fimage[0][i]) + fimage_logmin );
  }
  for( i_img = 0 ; i_img < NIMG; i_img++ ) { 
    image( fimage[i_img], ifnam[i_img+NIMG] );
  }

  /************************************************************************
    Output log versions:
  ************************************************************************/
  for( j = 0 ; j < 5; j++ ) { 
    for( i = 0 ; i < N1*N2; i++ ) { 
      fimage[0][i] = (double) failimage[j][i];
    }
    image( fimage[0], ifnam[j+2*NIMG] );
  }

  /* Reset array after every image dump: */
  for( i = 0 ; i < 5*N1*N2; i++ ) { 
    failimage[0][i] = 0;
  }

  return;
}

/******************************************************************************/
/******************************************************************************
  image(): 
  -----------
     -- wrapper for image_r8() and image_ppm()

******************************************************************************/
void image(double *f, char *fname)
{
#if( MAKE_PPM_IMAGE )
  image_ppm( f, fname );
#else 
  image_r8( f, fname );
#endif
  return;
}

/******************************************************************************/
/******************************************************************************
  image_r8(): 
  -----------
     -- generates a raw, unmapped "image" or pixelated output file. 

******************************************************************************/
void image_r8(double *f, char *fname)
{
  int i,j ;
  double iq,liq,scale,max,min,f_ncolors;
  FILE *fp;

  f_ncolors = (double) (NCOLORS-1);


  if( (fp = fopen(fname,"w")) == NULL ) {
    fflush(stderr);
    fprintf(stderr,"image(): Cannot open %s !! \n",fname);
    fflush(stderr);
    return;
  }

  /*  mapping is in NCOLORS steps lmax and lmin */
  max = min = f[0];
  for( i = 0 ; i < N1*N2; i++ ) { 
    if(f[i] > max) max = f[i] ;
    if(f[i] < min) min = f[i] ;
  }

  scale = f_ncolors/(max - min) ;

  /* Header information: */
  fprintf(fp, "RAW\n#  min=%g  , max=%g \n%d %d\n%d\n", min, max, N1, N2, (NCOLORS-1) );
  fflush(fp);

  for( i = 0 ; i < N1*N2; i++ ) { 
    iq = f[i];
    liq = scale*(iq - min) ;
    if(liq > f_ncolors) liq = f_ncolors ;
    if(liq < 0.) liq = 0. ;
    fprintf(fp,"%c",(unsigned char)((int)liq)) ;
  }

  fclose(fp);

  return;
}


/******************************************************************************/
/******************************************************************************
  image_ppm(): 
  -----------
     -- generates a color mapped "image" or pixelated output file following 
        the "raw" PPM  format ("man ppm" for more details). 

     -- color map determined by get_color_map()

******************************************************************************/
void image_ppm(double *f, char *fname)
{
  int i,j, icolor ;
  double iq, liq, scale, f_ncolors;
  double lmax, lmin;
  double gmax, gmin;
  int offset;
  int fail = 0;
  FILE *fp;

  MPI_File file;
  MPI_Status status;

  f_ncolors = (double) ppm_ncolors;

  if (WorldRank == 0) {
    if ((fp = fopen(fname, "w")) == NULL ) {
      fprintf(stderr, "image(): Cannot open %s !!\n", fname);
      fail = 1;
    }
  }

  MPI_Bcast(&fail, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (fail)
    return;

  /* mapping is in ppm_ncolors steps lmax and lmin */
  lmax = lmin = f[0];
  for( i = 0 ; i < N1*N2; i++ ) { 
    if(f[i] > lmax) lmax = f[i] ;
    if(f[i] < lmin) lmin = f[i] ;
  }

  /* need global minimum and maximum */
  MPI_Allreduce(&lmin, &gmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&lmax, &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  scale = f_ncolors/(gmax - gmin);

  /* Rank 0 writes header information: */
  if (WorldRank == 0) {
    fprintf(fp, "P6\n#  min=%g  , max=%g \n%d %d\n%d\n%n", gmin, gmax,
        GlobalN1, GlobalN2, ppm_ncolors, &offset);
    fclose(fp);
  }

  /* file offset information is required to set the view */
  MPI_Bcast(&offset, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* convert data to a string */
  char *data_as_text = (char *) malloc(RGB_SIZE * N1 * N2 * sizeof(char));
  int count = 0;
  for (i = 0 ; i < N1*N2; i++) {
    iq = f[i];
    liq = scale*(iq - gmin) ;
    icolor = (int) liq;
    if( icolor > ppm_ncolors ) { icolor = ppm_ncolors; }
    if( icolor < 0           ) { icolor = 0;           }

    sprintf(&data_as_text[count++], "%c", color_map[RED  ][icolor]);
    sprintf(&data_as_text[count++], "%c", color_map[GREEN][icolor]);
    sprintf(&data_as_text[count++], "%c", color_map[BLUE ][icolor]);
  }

  /* open the file, set the view, write data and close the file */
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
  MPI_File_set_view(file, offset, rgb, color_array, "native", MPI_INFO_NULL);
  MPI_File_write_all(file, data_as_text, N1 * N2, rgb, &status);
  MPI_File_close(&file);

  /* clean up */
  free(data_as_text);
  return;
}



/*******************************************************************************/
/*******************************************************************************
   get_color_map():
   ---------------

     -- reads in a plain-type ppm file to generate the color map 
     -- see  http://netpbm.sourceforge.net/doc/ppm.html  or "man ppm" 
          for PPM format specifications 

*******************************************************************************/
void get_color_map(void)
{

  int i, c0,c1,c2, iw, ibody;
  char c, word[100] ;
  int ipixel, icolor;
  FILE *fp; 
  const char magic_raw[]   = "P6";
  const char magic_plain[] = "P3";
  enum GET_type { get_magic, get_width, get_height, get_ncolors, get_pixels, get_newline };
  enum GET_type read_type, old_read_type;
  enum PPM_type { ppm_plain, ppm_raw };
  enum PPM_type ppm_type, ppm_body;

  if (WorldRank == 0) {
    fprintf(stdout,"Starting get_color_map().... \n");
    fflush(stdout);
  }
  
  if( (fp = fopen("map.ppm","r")) == NULL ) {
    fflush(stderr);
    fprintf(stderr,"get_color_map(): Cannot open map.ppm !! \n");
    fflush(stderr);
    exit(1);
  }

  old_read_type = read_type = get_magic;  /* Need to read in the magic number or ppm type first */
  ppm_type  = ppm_plain;                  /* Always read in the header as plain text */ 
  iw = 0;
  icolor = ipixel = 0;

  /* Read in the entire stream as a character stream */
  while( (c = getc(fp)) != EOF ) { 

    //fprintf(stderr,"%c",c);fflush(stderr);

    switch (read_type) { 

    case get_magic :  
      word[iw++] = c;
      if( iw == 2 ) { 
	word[iw] = '\0';
	if(      strcmp(word,magic_raw  ) == 0 ) {  ppm_body = ppm_raw;  }
	else if( strcmp(word,magic_plain) == 0 ) {  ppm_body = ppm_plain;}
	else { 
	  fprintf(stderr,"get_color_map(): error getting magic number, type = %s \n", word);
	  fflush(stderr);	    fclose(fp);	    exit(1);
	}
	old_read_type = read_type;
	read_type     = get_width;
	iw = 0;
      }
      break;
	
    case get_width :
      if( c == '#' ) { old_read_type = get_width;  read_type = get_newline; break; }
      if( isdigit(c) ) { word[iw++] = c; }
      else if( iw > 0 ) { 
	word[iw] = '\0';
	ppm_width = atoi(word);
	old_read_type = read_type;
	read_type     = get_height;
	iw = 0;
      }
      break;

    case get_height :
      if( c == '#' ) { old_read_type = get_height;  read_type = get_newline; break; }
      if( isdigit(c) ) { word[iw++] = c; }
      else if( iw > 0 ) { 
	word[iw] = '\0';
	ppm_height    = atoi(word);
	old_read_type = read_type;
	read_type     = get_ncolors;
	iw = 0;
      }
      break;

    case get_ncolors :
      if( c == '#' ) { old_read_type = get_ncolors;  read_type = get_newline; break; }
      if( isdigit(c) ) { word[iw++] = c; }
      else if( iw > 0 ) { 
	word[iw] = '\0';
	ppm_ncolors   = atoi(word);
	old_read_type = read_type;
	read_type     = get_pixels;
	iw = 0;
	if( ppm_ncolors > NCOLORS ) { 
	  fprintf(stderr,
		  "get_color_map(): too many colors: ncolors = %d,  maxcolors = %d \n", 
		  ppm_ncolors, NCOLORS);
	  fflush(stderr);	    fclose(fp);	    exit(1);
	}
      }
      break;

    case get_pixels :
      /* We assume that pixels come in threes and there are no comments within a pixel */
      if( c == '#' ) { old_read_type = get_pixels;  read_type = get_newline; break; }
      if( isdigit(c) ) { word[iw++] = c; }
      else if( iw > 0 ) { 
	word[iw] = '\0';
	color_map[icolor++][ipixel] = atoi(word);   
	iw = 0;
	if( icolor == 3 ) { 
	  icolor = 0 ; 
	  ipixel++;
	}
      }
      break;

    case get_newline :
      if( c == '\n' ) { 
	read_type = old_read_type;
      }
      break;

    }  /* switch off */

  }  /* end while */

  if( ipixel != (ppm_width*ppm_height) ) { 
    fprintf(stderr,
	    "get_color_map(): missing pixels: npixel = %d,  width = %d,  height = %d \n", 
	    ipixel, ppm_width, ppm_height);
    fflush(stderr);	    fclose(fp);	    exit(1);
  }

  if( (ppm_ncolors+1) != ipixel  ) { 
    fprintf(stderr,"get_color_map(): missing colors: npixel = %d,  ppm_ncolors = %d \n",
	    ipixel, ppm_ncolors);
    fflush(stderr);	    fclose(fp);	    exit(1);
  }

  fclose(fp);


  /********************************************************************************
    test get_color_map():  Make sure that testmap.ppm looks identical to your map.ppm file. 
  ********************************************************************************/
  // Can probably comment this out - no need for testing... I think
  /*if( (fp = fopen("testmap.ppm","w")) == NULL ) {
    exit(1);
  }
  fprintf(fp,"P6\n%d %d\n%d\n",ppm_width, ppm_height,ppm_ncolors);  fflush(fp);
  for( i = 0 ; i <= ppm_ncolors; i++ ) { 
    fputc(color_map[RED  ][i],fp);
    fputc(color_map[GREEN][i],fp);
    fputc(color_map[BLUE ][i],fp);
  }
  fclose(fp);*/

  if (WorldRank == 0) {
    fprintf(stdout, "Ending get_color_map().... \n");
    fflush(stdout);
  }

  return;
}


#undef MAKE_PPM_IMAGE
#undef NCOLORS 
#undef IRHO  
#undef IUU   
#undef IBSQ  
#undef IGAM  
#undef RED   
#undef GREEN 
#undef BLUE  
#undef IMGEXT
