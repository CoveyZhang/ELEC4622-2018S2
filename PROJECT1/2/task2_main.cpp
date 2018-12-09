/*****************************************************************************/
// File: vertical_filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include<time.h>

#include "io_bmp.h"
#include "aligned_image_comps.h"
#include <math.h>
#include <time.h>


/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*              my_aligned_image_comp::perform_boundary_extension            */
/*****************************************************************************/

void my_aligned_image_comp::perform_boundary_extension()
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

void my_aligned_image_comp::bilinear(my_aligned_image_comp *in)
{
	for (int i = 0;i < height+1;i++) {
		for (int j = 0;j < width+1;j++) {
#if 0
			if (i % 2 == 0 && j % 2 == 0)
			{
				float *ip = in->buf + (i*3/5)*in->stride + j*3/5;
				float *op = this->buf + (i)*stride + j;
				*op = ip[0];
			}
			if (i % 2 == 1 && j % 2 == 0)
			{
				float *ip = in->buf + (i*3/5)*in->stride + j*3/5;
				float *op = this->buf + (i)*stride + j;
				*op = (ip[0] + ip[in->stride]) / 2;
			}
			if (i % 2 == 0 && j % 2 == 1)
			{
				float *ip = in->buf + (i*3/5)*in->stride + j*3/5;
				float *op = this->buf + (i)*stride + j;
				*op = (ip[0] + ip[1]) / 2;
			}
			if (i % 2 == 1 && j % 2 == 1)
			{
				float *ip = in->buf + (i*3/5)*in->stride + j*3/5;
				float *op = this->buf + (i)*stride + j;
				*op = (ip[0] + ip[in->stride] + ip[1] + ip[in->stride + 1]) / 4;
			}
#else			
			int ii = (i) * 0.6F;
			int jj = (j) * 0.6F;
			float *ip = in->buf + (ii)*in->stride + jj;
			float *op = this->buf + (i)*stride + j;
			float y = float(i) * 0.6F - ii;
			float x = float(j) * 0.6F - jj;
			*op = ip[0] * (1 - x)*(1 - y) + ip[in->stride] * y*(1 - x) + ip[1] * (1 - y)*x + ip[in->stride + 1] * x*y;
#endif
		}
	}
}

/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
	clock_t start_time = clock();

  if (argc != 3)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file>\n",argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

      int width = in.cols, height = in.rows;
      int n, num_comps = in.num_components;

      my_aligned_image_comp *input_comps =
        new my_aligned_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        input_comps[n].init(height,width,2); 
      
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
	  delete[] line;
      bmp_in__close(&in);

	  my_aligned_image_comp *output_comps =
		  new my_aligned_image_comp[num_comps];
	  

	  int out_width, out_height;
	  out_width = (int)ceilf(width * 1.6667F);
	  out_height = (int)ceilf(height * 1.6667F);

	  for (n = 0; n < num_comps; n++) {
		  output_comps[n].init(out_height, out_width, 1); // Don't need a border for output
	  }

	  for (n = 0; n < num_comps; n++) {
		  input_comps[n].perform_boundary_extension();
	  }
	  printf("Start!\r\n");

	  for (n = 0; n < num_comps; n++) {
		  output_comps[n].bilinear(input_comps + n);
	  }
	  
	  printf("Finish!\r\n");


      // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2],out_width,out_height,num_comps)) != 0)
        throw err_code;

	  io_byte *oline = new io_byte[out_width*num_comps];

      for (r=out_height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = oline+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
			  for (int c = 0; c < out_width; c++, dst += num_comps) {
				  if (src[c] > 255) {
					  src[c] = 255;
				  }
				  else if (src[c] < 0) {
					  src[c] = 0;
				  }
				  *dst = (io_byte)src[c]; // The cast to type "io_byte" is
						// required here, since floats cannot generally be
						// converted to bytes without loss of information.  The
						// compiler will warn you of this if you remove the cast.
						// There is in fact not the best way to do the
						// conversion.  You should fix it up in the lab.
			  }
            }
          bmp_out__put_line(&out,oline);
        }
      bmp_out__close(&out);
      delete[] oline;
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
