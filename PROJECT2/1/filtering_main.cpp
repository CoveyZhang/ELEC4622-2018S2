/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#include <math.h>
#include <time.h>
#include <cmath>

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension()
{
  int r, c;

  // First extend upwards
  float *first_line = buf;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      first_line[-r*stride+c] = first_line[c];

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      last_line[r*stride+c] = last_line[c];

  // Now extend all rows to the left and to the right
  float *left_edge = buf-border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
        left_edge[-c] = left_edge[0];
        right_edge[c] = right_edge[0];
      }
}


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/

void apply_filter(my_image_comp *in, my_image_comp *out, float sigma, float alpha)
{
#define PI 3.1415926F
	int FILTER_EXTENT = ceil(3 * sigma);
	int FILTER_DIM = (2 * FILTER_EXTENT + 1);
	int FILTER_TAPS = (FILTER_DIM*FILTER_DIM);

  // Create the filter kernel as a local array on the stack, which can accept
  // row and column indices in the range -FILTER_EXTENT to +FILTER_EXTENT.
  float *filter_buf = new float[FILTER_TAPS];
  float *mirror_psf = filter_buf+(FILTER_DIM*FILTER_EXTENT)+FILTER_EXTENT;
          // `mirror_psf' points to the central tap in the filter
  int r, c;
  for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++)
	  for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++)
		  mirror_psf[r*FILTER_DIM + c] = ((r*r + c*c) / (sigma*sigma) - 2)*exp(-(r*r + c*c) / (2*sigma*sigma))/(2*PI*pow(sigma,4));

  // Check for consistent dimensions
  assert(in->border >= FILTER_EXTENT);
  assert((out->height <= in->height) && (out->width <= in->width));

  // Perform the convolution
  for (r=0; r < out->height; r++)
    for (c=0; c < out->width; c++)
      {
        float *ip = in->buf + r*in->stride + c;
        float *op = out->buf + r*out->stride + c;
        float sum = 0.0F;
        for (int y=-FILTER_EXTENT; y <= FILTER_EXTENT; y++)
          for (int x=-FILTER_EXTENT; x <= FILTER_EXTENT; x++)
            sum += ip[y*in->stride+x] * mirror_psf[y*FILTER_DIM+x];
        *op = 128 + sum*alpha;
      }
}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if (argc != 5)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> <sigma> <alpha>\n",argv[0]);
      return -1;
    }

  int err_code=0;
  
  clock_t start_time = clock();

  try {
      // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

	  float sigma = atof(argv[3]);
	  float alpha = atof(argv[4]);

      int width = in.cols, height = in.rows;
      int n, num_comps = in.num_components;
      my_image_comp *input_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        input_comps[n].init(height,width,ceil(3*sigma)); // Leave a border of 4
      
      int r; // Declare row index
      io_byte *line = new io_byte[width*num_comps];
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&in,line)) != 0)
            throw err_code;
          for (n=0; n < num_comps; n++)
            {
              io_byte *src = line+n; // Points to first sample of component n
              float *dst = input_comps[n].buf + r * input_comps[n].stride;
              for (int c=0; c < width; c++, src+=num_comps)
                dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
            }
        }
      bmp_in__close(&in);

      // Allocate storage for the filtered output
      my_image_comp *output_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        output_comps[n].init(height,width,0); // Don't need a border for output

      // Process the image, all in floating point (easy)
      for (n=0; n < num_comps; n++)
        input_comps[n].perform_boundary_extension();

	  printf("Start!\r\n");

      for (n=0; n < num_comps; n++)
        apply_filter(input_comps+n,output_comps+n,sigma,alpha);

	  printf("Finish!\r\n");

      // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2],width,height,num_comps)) != 0)
        throw err_code;
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = line+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
              for (int c=0; c < width; c++, dst+=num_comps)
			  {
				  if (src[c] > 255) {
					  src[c] = 255;
				  }
				  else if (src[c] < 0) {
					  src[c] = 0;
				  }
				  *dst = (io_byte)src[c];
			  }
            }
          bmp_out__put_line(&out,line);
        }
      bmp_out__close(&out);
      delete[] line;
      delete[] input_comps;
	  delete[] output_comps;

	  clock_t end_time = clock();
	  float elaps = ((float)(end_time - start_time)) / CLOCKS_PER_SEC;
	  printf_s("The runtime is %f seconds!\n\r", elaps);

    }
  catch (int exc) {
      if (exc == IO_ERR_NO_FILE)
        fprintf(stderr,"Cannot open supplied input or output file.\n");
      else if (exc == IO_ERR_FILE_HEADER)
        fprintf(stderr,"Error encountered while parsing BMP file header.\n");
      else if (exc == IO_ERR_UNSUPPORTED)
        fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
      else if (exc == IO_ERR_FILE_TRUNC)
        fprintf(stderr,"Input or output file truncated unexpectedly.\n");
      else if (exc == IO_ERR_FILE_NOT_OPEN)
        fprintf(stderr,"Trying to access a file which is not open!(?)\n");
      return -1;
    }
  return 0;
}
