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

void my_aligned_image_comp::difference(my_aligned_image_comp *in1, my_aligned_image_comp *in2)
{
	float sum = 0,sum1 = 0,temp = 0;
	for (int i = 0;i < height;i++) {
		for (int j = 0;j < width;j++) {
			float *ip1 = in1->buf + (i)*in1->stride + j;
			float *ip2 = in2->buf + (i)*in2->stride + j;
			float *op = this->buf + (i)*stride + j;
			temp = (ip1[0] - ip2[0]);
			*op = 128 + 0.5 * temp;
			sum += temp;
			sum1 += temp*temp;
		}
	}
	float me = sum / (height*width);
	float mse = sum1 / (height*width);
	float psnr = 10 * log10(255*255/mse);
	printf("ME:%f\r\n",me);
	printf("MSE:%f\r\n",mse);
	printf("PSNR:%f\r\n",psnr);
	printf("\r\n");
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

  if (argc != 4)
    {
      fprintf(stderr,"Usage: %s <in bmp file 1> <in bmp file 2><out bmp file>\n",argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      // Read the input image
      bmp_in in1,in2;
      if ((err_code = bmp_in__open(&in1,argv[1])) != 0)
        throw err_code;

      int width1 = in1.cols, height1 = in1.rows;
      int n, num_comps = in1.num_components;

      my_aligned_image_comp *input_comps1 =
        new my_aligned_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        input_comps1[n].init(height1,width1,0); 
      
      int r1; // Declare row index
      io_byte *line = new io_byte[width1*num_comps];
      for (r1=height1-1; r1 >= 0; r1--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&in1,line)) != 0)
            throw err_code;
          for (n=0; n < num_comps; n++)
            {
              io_byte *src = line+n; // Points to first sample of component n
              float *dst = input_comps1[n].buf + r1 * input_comps1[n].stride;
              for (int c=0; c < width1; c++, src+=num_comps)
                dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
            }
        }
	  delete[] line;
      bmp_in__close(&in1);

	  if ((err_code = bmp_in__open(&in2, argv[2])) != 0)
		  throw err_code;

	  int width2 = in2.cols, height2 = in2.rows;
	  //int n, num_comps = in2.num_components;

	  my_aligned_image_comp *input_comps2 =
		  new my_aligned_image_comp[num_comps];
	  for (n = 0; n < num_comps; n++)
		  input_comps2[n].init(height2, width2, 0);

	  int r2; // Declare row index
	  io_byte *line2 = new io_byte[width2*num_comps];
	  for (r2 = height1 - 1; r2 >= 0; r2--)
	  { // "r" holds the true row index we are reading, since the image is
		// stored upside down in the BMP file.
		  if ((err_code = bmp_in__get_line(&in2, line2)) != 0)
			  throw err_code;
		  for (n = 0; n < num_comps; n++)
		  {
			  io_byte *src = line2 + n; // Points to first sample of component n
			  float *dst = input_comps2[n].buf + r2 * input_comps2[n].stride;
			  for (int c = 0; c < width2; c++, src += num_comps)
				  dst[c] = (float)*src; // The cast to type "float" is not
										// strictly required here, since bytes can always be
										// converted to floats without any loss of information.
		  }
	  }
	  delete[] line2;
	  bmp_in__close(&in2);

	  my_aligned_image_comp *output_comps =
		  new my_aligned_image_comp[num_comps];
	  

	  int minw, minh;
	  if (width1 > width2) {
		  minw = width2;
	  }
	  else
	  {
		  minw = width1;
	  }
	  if (height1 > height2) {
		  minh = height2;
	  }
	  else
	  {
		  minh = height1;
	  }
	  int out_width, out_height;
	  out_width = minw;
	  out_height = minh;

	  for (n = 0; n < num_comps; n++) {
		  output_comps[n].init(out_height, out_width, 0); // Don't need a border for output
	  }



	  printf("Start!\r\n");

	  for (n = 0; n < num_comps; n++) {
		  printf("This is plane %d.\r\n",n+1);
		  output_comps[n].difference(input_comps1 + n, input_comps2 + n);
	  }
	  
	  printf("Finish!\r\n");


      // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[3],out_width,out_height,num_comps)) != 0)
        throw err_code;

	  io_byte *oline = new io_byte[out_width*num_comps];
	  int r;
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
	  delete[] input_comps1;
	  delete[] input_comps2;
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
