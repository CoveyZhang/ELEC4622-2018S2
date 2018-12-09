/*****************************************************************************/
// File: motion_main.cpp
// Author: David Taubman
// Last Revised: 30 September, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#include "motion.h"
#include <math.h>

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension()
{
#if 0
  int r, c;
  // First extend upwards
  int *first_line = buf;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      first_line[-r*stride+c] = first_line[c];

  // Now extend downwards
  int *last_line = buf+(height-1)*stride;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      last_line[r*stride+c] = last_line[c];

  // Now extend all rows to the left and to the right
  int *left_edge = buf-border*stride;
  int *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
        left_edge[-c] = left_edge[0];
        right_edge[c] = right_edge[0];
      }
#else
	int r, c;

	// First extend upwards
	int *first_line = buf;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			first_line[-r*stride + c] = 2*first_line[c]-first_line[(r-1)*stride+c];

	// Now extend downwards
	int *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			last_line[r*stride + c] = 2*last_line[c]-last_line[(1-r)*stride+c];

	// Now extend all rows to the left and to the right
	int *left_edge = buf - border*stride;
	int *right_edge = left_edge + width - 1;
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
		for (c = 1; c <= border; c++)
		{
			left_edge[-c] = 2*left_edge[0]-left_edge[c-1];
			right_edge[c] = 2*right_edge[0]-right_edge[-c+1];
		}
#endif
}

/*****************************************************************************/
/* STATIC                         find_motion                                */
/*****************************************************************************/

static mvector
  find_motion(my_image_comp *ref, my_image_comp *tgt,
              int start_row, int start_col, int block_width, int block_height)
  /* This function finds the motion vector which best describes the motion
     between the `ref' and `tgt' frames, over a specified block in the
     `tgt' frame.  Specifically, the block in the `tgt' frame commences
     at the coordinates given by `start_row' and `start_col' and extends
     over `block_width' columns and `block_height' rows.  The function finds
     the translational offset (the returned vector) which describes the
     best matching block of the same size in the `ref' frame, where
     the "best match" is interpreted as the one which minimizes the sum of
     absolute differences (SAD) metric. */
{
  mvector vec, best_vec;
  int mse, best_mse=256*256;
  for (vec.y=-8; vec.y <= 8; vec.y++)
    for (vec.x=-8; vec.x <= 8; vec.x++)
      {
        int ref_row = start_row-vec.y;
        int ref_col = start_col-vec.x;
        if ((ref_row < 0) || (ref_col < 0) ||
            ((ref_row+block_height) > ref->height) ||
            ((ref_col+block_width) > ref->width))
          continue; // Translated block not containe within reference frame
        int r, c;
        int *rp = ref->buf + ref_row*ref->stride + ref_col;
        int *tp = tgt->buf + start_row*tgt->stride + start_col;
        for (mse=0, r=block_height; r > 0; r--,
             rp+=ref->stride, tp+=tgt->stride)
          for (c=0; c < block_width; c++)
            {
              int diff = tp[c] - rp[c];
			  mse += diff*diff;
            }
		mse /= (block_height*block_width);
        if (mse < best_mse)
          {
            best_mse = mse;
            best_vec = vec;
          }
      }

  return best_vec;
}

/*****************************************************************************/
/* STATIC                         motion_comp                                */
/*****************************************************************************/

static void
  motion_comp(my_image_comp *ref, my_image_comp *tgt, mvector vec,
              int start_row, int start_col, int block_width, int block_height)
  /* This function transfers data from the `ref' frame to a block within the
     `tgt' frame, thereby realizing motion compensation.  The motion in
     question has already been found by `find_motion' and is captured by
     the `vec' argument.  The block in the `tgt' frame commences
     at the coordinates given by `start_row' and `start_col' and extends
     for `block_width' columns and `block_height' rows. */
{
  int r, c;
  int ref_row = start_row - vec.y;
  int ref_col = start_col - vec.x;
  int *rp = ref->buf + ref_row*ref->stride + ref_col;
  int *tp = tgt->buf + start_row*tgt->stride + start_col;
  for (r=block_height; r > 0; r--,
       rp+=ref->stride, tp+=tgt->stride)
    for (c=0; c < block_width; c++)
      tp[c] = rp[c];
}

static void
motion_comp_global_step1(my_image_comp *ref, my_image_comp *tgt, global_vector vec,
	int start_row, int start_col, int block_width, int block_height)
{
	int r, c;
	int ref_row = start_row - vec.y_int;
	int ref_col = start_col - vec.x_int;
	int *rp = ref->buf + ref_row*ref->stride + ref_col;
	int *tp = tgt->buf + start_row*tgt->stride + start_col;
	for (r = block_height; r > 0; r--,
		rp += ref->stride, tp += tgt->stride)
		for (c = 0; c < block_width; c++)
			tp[c] = rp[c];
}

static void
motion_comp_global_step2(my_image_comp *in, my_image_comp *out, global_vector vec)
{
	int r, c;
	float x = vec.x_float;
	float y = vec.y_float;
	float temp;
	if (x > 0 && y > 0)
	{
		for (r = 0; r < out->height; r++)
			for (c = 0; c < out->width; c++)
			{
				int *ip = in->buf + r*in->stride + c;
				int *op = out->buf + r*out->stride + c;
				temp = (1 - x)*(1 - y)*ip[0] + x*(1 - y)*ip[-1] + (1 - x)*y*ip[-in->stride] + x*y*ip[-in->stride - 1];
				*op = int(temp + 0.5);
			}
	}
	if (x > 0 && y < 0)
	{
		for (r = 0; r < out->height; r++)
			for (c = 0; c < out->width; c++)
			{
				int *ip = in->buf + r*in->stride + c;
				int *op = out->buf + r*out->stride + c;
				temp = (1 - x)*(1 + y)*ip[0] + x*(1 + y)*ip[-1] + (1 - x)*(-y)*ip[in->stride] + x*(-y)*ip[in->stride - 1];
				*op = int(temp + 0.5);
			}
	}
	if (x < 0 && y > 0)
	{
		for (r = 0; r < out->height; r++)
			for (c = 0; c < out->width; c++)
			{
				int *ip = in->buf + r*in->stride + c;
				int *op = out->buf + r*out->stride + c;
				temp = (1 + x)*(1 - y)*ip[0] + (-x)*(1 - y)*ip[1] + (1 + x)*y*ip[-in->stride] + (-x)*y*ip[-in->stride + 1];
				*op = int(temp + 0.5);
			}
	}
	if (x < 0 && y < 0)
	{
		for (r = 0; r < out->height; r++)
			for (c = 0; c < out->width; c++)
			{
				int *ip = in->buf + r*in->stride + c;
				int *op = out->buf + r*out->stride + c;
				temp = (1 + x)*(1 + y)*ip[0] + (-x)*(1 + y)*ip[1] + (-y)*(1 + x)*ip[in->stride] + (-x)*(-y)*ip[in->stride + 1];
				*op = int(temp + 0.5);
			}
	}
}

float display_mse(my_image_comp *tgt, my_image_comp *out) {
	float sum = 0, temp = 0;
	for (int i = 0;i < tgt->height;i++) {
		for (int j = 0;j < tgt->width;j++) {
			int *ip1 = tgt->buf + i*tgt->stride + j;
			int *ip2 = out->buf + i*out->stride + j;
			temp = (ip1[0] - ip2[0]);
			sum += temp*temp;
		}
	}
	float mse = sum / (tgt->height*tgt->width);
	return mse;
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
  if (argc != 4)
    {
      fprintf(stderr,
              "Usage: %s <bmp frame 1> <bmp frame 2> <bmp MC out>\n",
              argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      // Read the input image
      bmp_in in[2];
      if ((err_code = bmp_in__open(&in[0],argv[1])) != 0)
        throw err_code;
      if ((err_code = bmp_in__open(&in[1],argv[2])) != 0)
        throw err_code;

      int width = in[0].cols, height = in[0].rows;
      if ((width != in[1].cols) || (height != in[1].rows))
        {
          fprintf(stderr,"The two input frames have different dimensions.\n");
          return -1;
        }
      my_image_comp mono[2];
      mono[0].init(height,width,8); // Leave a border of 4 (in case needed)
      mono[1].init(height,width,8); // Leave a border of 4 (in case needed)
	  
      int n, r, c;
      int num_comps = in[0].num_components;
      io_byte *line = new io_byte[width*num_comps];
      for (n=0; n < 2; n++)
        {
          for (r=height-1; r >= 0; r--)
            { // "r" holds the true row index we are reading, since the image
              // is stored upside down in the BMP file.
              if ((err_code = bmp_in__get_line(&(in[n]),line)) != 0)
                throw err_code;
              io_byte *src = line; // Points to first sample of component n
              int *dst = mono[n].buf + r * mono[n].stride;
              for (c=0; c < width; c++, src+=num_comps)
                dst[c] = *src;
            }
          bmp_in__close(&(in[n]));
        }

	  mono[0].perform_boundary_extension();
	  mono[1].perform_boundary_extension();

      // Allocate storage for the motion compensated output
      my_image_comp output;
      output.init(height,width,0); 

	  my_image_comp output_new;
	  output_new.init(height, width, 1); 

	  my_image_comp output_newest;
	  output_newest.init(height, width, 0);

	  int *buf_temp = output_new.buf;
	  int *buf_tgt = mono[1].buf;
	  for (r = 0; r < height; r++)
	  {
		  for (c = 0;c < width; c++)
		  {
			  int *op = buf_temp + r*output_new.stride + c;
			  int *ip = buf_tgt + r*mono[1].stride + c;
			  *op = *ip;
		  }
	  }

      // Now perform simple motion estimation and compensation
      int nominal_block_width = 8;
      int nominal_block_height = 8;
      int block_width, block_height;
	  global_vector g_vec;
	  int num=0;
      for (r=0; r < height; r+=block_height)
        {
          block_height = nominal_block_height;
          if ((r+block_height) > height)
            block_height = height-r;
          for (c=0; c < width; c+=block_width)
            {
              block_width = nominal_block_width;
              if ((c+block_width) > width)
                block_width = width-c;
              mvector vec = find_motion(&(mono[0]),&(mono[1]),
                                        r,c,block_width,block_height);
              motion_comp(&(mono[0]),&output,vec,
                          r,c,block_width,block_height);
			  g_vec.x += vec.x;
			  g_vec.y += vec.y;
			  num++;
			  //printf("x:%d ; y:%d ; cnt:%d\n",vec.x,vec.y,counter );
            }
        }
	  g_vec.x /= num;
	  g_vec.y /= num;
	  g_vec.calculate();
	  float MSE_O = display_mse(&mono[1], &output);
	  printf("Block MSE: %f\n", MSE_O);

	  printf("Global motion estimation vector (x,y): (%f,%f)\n",g_vec.x,g_vec.y);
	  printf("integer part:(%d,%d)\n", g_vec.x_int, g_vec.y_int);
	  printf("float part:(%f,%f)\n", g_vec.x_float, g_vec.y_float);

	  for (r = 0; r < height; r += block_height)
	  {
		  block_height = nominal_block_height;
		  if ((r + block_height) > height)
			  block_height = height - r;
		  for (c = 0; c < width; c += block_width)
		  {
			  block_width = nominal_block_width;
			  if ((c + block_width) > width)
				  block_width = width - c;
			  motion_comp_global_step1(&(mono[0]), &output_new, g_vec,
				  r, c, block_width, block_height);
		  }
	  }

	  output_new.perform_boundary_extension();

	  motion_comp_global_step2(&output_new, &output_newest, g_vec);

	  float MSE_N = display_mse(&mono[1], &output_newest);
	  printf("Global MSE: %f\n", MSE_N);

	  float psnr = 10 * log10(255 * 255 / MSE_N);
	  printf("Global psnr: %f\n", psnr);



      // Write the motion compensated image out
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[3],width,height,1)) != 0)
        throw err_code;
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
		  
			io_byte *dst = line; // Points to first sample of component n
			int *src = output_newest.buf + r * output_newest.stride;
			for (int c = 0; c < width; c++, dst++) {
				if (src[c] > 255) {
					src[c] = 255;
				}
				if (src[c] < 0) {
					src[c] = 0;
				}
				*dst = (io_byte)src[c];
			 }
          bmp_out__put_line(&out,line);
        }
      bmp_out__close(&out);
      delete[] line;
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
