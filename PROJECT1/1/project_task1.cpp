/*****************************************************************************/
// File: vertical_filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "aligned_image_comps1.h"
#include <math.h>
#include <time.h>

#define v_mode 1
#define h_mode 2

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
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			first_line[-r*stride + c] = first_line[c];

	// Now extend downwards
	float *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			last_line[r*stride + c] = last_line[c];

	// Now extend all rows to the left and to the right
	float *left_edge = buf - border*stride;
	float *right_edge = left_edge + width - 1;
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
		for (c = 1; c <= border; c++)
		{
			left_edge[-c] = left_edge[0];
			right_edge[c] = right_edge[0];
		}
}

/*****************************************************************************/
/*              my_aligned_image_comp::filter (hanning window)               */
/*****************************************************************************/

void my_aligned_image_comp::filter(my_aligned_image_comp *in, int filter_length, int mode)
{
#define FILTER_EXTENT 14
#define FILTER_TAPS (2*FILTER_EXTENT+1)
#define PI 3.141592653589793F

	// Create the vertical filter PSF as a local array on the stack.
	float filter_buf0[FILTER_TAPS];
	float filter_buf1[FILTER_TAPS];
	float filter_buf2[FILTER_TAPS];

	for (int i = 0; i < FILTER_TAPS; ++i) {
		filter_buf0[i] = 0.0F;
		filter_buf1[i] = 0.0F;
		filter_buf2[i] = 0.0F;
	}

	float *mirror_psf0 = filter_buf0 + FILTER_EXTENT;
	float *mirror_psf1 = filter_buf1 + FILTER_EXTENT;
	float *mirror_psf2 = filter_buf2 + FILTER_EXTENT;

	// `mirror_psf' points to the central tap in the filter
	for (int t = -filter_length; t <= filter_length; t++)
	{
		mirror_psf1[t] = 0.6f*sinf(0.6f*PI*(t + 0.3333F)) / (0.6f*PI*(t + 0.3333F))*0.5F*(1 + cosf(PI*(t + 0.3333F) / (filter_length + 0.5F)));
		mirror_psf2[t] = 0.6f*sinf(0.6f*PI*(t - 0.3333F)) / (0.6f*PI*(t - 0.3333F))*0.5F*(1 + cosf(PI*(t - 0.3333F) / (filter_length + 0.5F)));

		if (t == 0)
			mirror_psf0[t] = 0.6f;
		/*
		else if (t > 0)
		mirror_psf0[t] = mirror_psf0[-t];
		*/
		else
			mirror_psf0[t] = 0.6f*sinf(0.6f*PI*t) / (0.6f*PI*t)*0.5F*(1 + cosf((PI*t) / (filter_length + 0.5F)));


	}
	float gain0 = 0, gain1 = 0, gain2 = 0;
	for (int t = -filter_length; t <= filter_length; t++)
	{
		gain0 += mirror_psf0[t];
		gain1 += mirror_psf1[t];
		gain2 += mirror_psf2[t];
	}

	for (int t = -filter_length; t <= filter_length; t++) {
		mirror_psf0[t] = mirror_psf0[t] / gain0;
		mirror_psf1[t] = mirror_psf1[t] / gain1;
		mirror_psf2[t] = mirror_psf2[t] / gain2;
	}

	// Check for consistent dimensions
	//assert(in->border >= FILTER_EXTENT);
	//assert((this->height <= in->height) && (this->width <= in->width));

	// Perform the convolution
	if (mode == 1) {
		for (int r = 0, rr = 0;(r < height) && (rr < in->height); r += 3, rr += 5)
			for (int c = 0; c < width; c++)
			{
				float *ip = in->buf + rr*in->stride + c;
				float *op = buf + r*stride + c;
				float sum = 0.0F;
				for (int y = -filter_length; y <= filter_length; y++)
					sum += ip[y*in->stride] * mirror_psf0[y];
				*op = sum;
			}
		for (int r = 0, rr = 0;(r < (height - 1)) && (rr < (in->height - 2)); r += 3, rr += 5)
			for (int c = 0; c < width; c++)
			{
				float *ip = in->buf + ((rr + 2)*in->stride) + c;
				float *op = buf + (r + 1)*stride + c;
				float sum = 0.0F;
				for (int y = -filter_length; y <= filter_length; y++)
					sum += ip[y*in->stride] * mirror_psf1[y];
				*op = sum;
			}
		for (int r = 0, rr = 0;(r < (height - 2)) && (rr < (in->height - 3)); r += 3, rr += 5)
			for (int c = 0; c < width; c++)
			{
				float *ip = in->buf + ((rr + 3)*in->stride) + c;
				float *op = buf + (r + 2)*stride + c;
				float sum = 0.0F;
				for (int y = -filter_length; y <= filter_length; y++)
					sum += ip[y*in->stride] * mirror_psf2[y];
				*op = sum;
			}
	}
	if (mode == 2) {
		for (int r = 0; r < height; r++) {
			for (int c = 0, cc = 0; (c < width) && (cc < in->stride); c += 3, cc += 5) {
				float *ip = in->buf + r*in->stride + cc;
				float *op = buf + r*stride + c;
				float sum = 0.0F;
				for (int y = -filter_length; y <= filter_length; y++)
					sum += ip[y] * mirror_psf0[y];
				*op = sum;
			}

			for (int c = 0, cc = 0; (c < (width - 1)) && (cc < (in->stride - 2)); c += 3, cc += 5) {
				float *ip = in->buf + r*in->stride + cc + 2;
				float *op = buf + r*stride + c + 1;
				float sum = 0.0F;
				for (int y = -filter_length; y <= filter_length; y++)
					sum += ip[y] * mirror_psf1[y];
				*op = sum;
			}

			for (int c = 0, cc = 0; (c < (width - 2)) && (cc < (in->stride - 3)); c += 3, cc += 5) {
				float *ip = in->buf + r*in->stride + cc + 3;
				float *op = buf + r*stride + c + 2;
				float sum = 0.0F;
				for (int y = -filter_length; y <= filter_length; y++)
					sum += ip[y] * mirror_psf2[y];
				*op = sum;
			}
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
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <in bmp file> <out bmp file> <filter length>\n", argv[0]);
		return -1;
	}

	clock_t start_time = clock();

	int err_code = 0;
	try {
		// Read the input image
		bmp_in in;
		if ((err_code = bmp_in__open(&in, argv[1])) != 0)
			throw err_code;

		int filter_len = atoi(argv[3]);
		int width = in.cols, height = in.rows;
		int n, num_comps = in.num_components;
		my_aligned_image_comp *input_comps =
			new my_aligned_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			input_comps[n].init(height, width, 15); // Leave a border of 4

		int r; // Declare row index
		io_byte *line = new io_byte[width*num_comps];
		for (r = height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are reading, since the image is
		  // stored upside down in the BMP file.
			if ((err_code = bmp_in__get_line(&in, line)) != 0)
				throw err_code;
			for (n = 0; n < num_comps; n++)
			{
				io_byte *src = line + n; // Points to first sample of component n
				float *dst = input_comps[n].buf + r * input_comps[n].stride;
				for (int c = 0; c < width; c++, src += num_comps)
					dst[c] = (float)*src; // The cast to type "float" is not
										  // strictly required here, since bytes can always be
										  // converted to floats without any loss of information.
			}
		}
		bmp_in__close(&in);
		// Change output size
		int owidth, oheight;
		owidth = (int)floorf(width * 0.6F);
		oheight = (int)floorf(height * 0.6F);
		// Allocate storage for the filtered output
		my_aligned_image_comp *output_comps =
			new my_aligned_image_comp[num_comps];
		my_aligned_image_comp * intermediea_comps =
			new my_aligned_image_comp[num_comps];
		for (n = 0; n < num_comps; n++) {
			intermediea_comps[n].init(oheight, width, 15);
			output_comps[n].init(oheight, owidth, 0); // Don't need a border for output
		}
		for (n = 0; n < num_comps; n++) {
			input_comps[n].perform_boundary_extension();
			intermediea_comps[n].perform_boundary_extension();
		}
		printf("Start filtering!\r\n");

		// Process the image, all in floating point (easy)

#if 1
		for (n = 0; n < num_comps; n++)
			intermediea_comps[n].filter(input_comps + n, filter_len, v_mode);

		for (n = 0; n < num_comps; n++) {
			intermediea_comps[n].perform_boundary_extension();
			output_comps[n].filter(intermediea_comps + n, filter_len, h_mode);
		}
#else
		for (n = 0; n < num_comps; n++)
			output_comps[n].vector_filter(input_comps + n);
#endif

		// Write the image back out again
		printf("Filtering end!\r\n");
		bmp_out out;
		if ((err_code = bmp_out__open(&out, argv[2], owidth, oheight, num_comps)) != 0)
			throw err_code;
		for (r = oheight - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			for (n = 0; n < num_comps; n++)
			{
				io_byte *dst = line + n; // Points to first sample of component n
				float *src = output_comps[n].buf + r * output_comps[n].stride;
				for (int c = 0; c < owidth; c++, dst += num_comps) {
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
			bmp_out__put_line(&out, line);
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
			fprintf(stderr, "Cannot open supplied input or output file.\n");
		else if (exc == IO_ERR_FILE_HEADER)
			fprintf(stderr, "Error encountered while parsing BMP file header.\n");
		else if (exc == IO_ERR_UNSUPPORTED)
			fprintf(stderr, "Input uses an unsupported BMP file format.\n  Current "
				"simple example supports only 8-bit and 24-bit data.\n");
		else if (exc == IO_ERR_FILE_TRUNC)
			fprintf(stderr, "Input or output file truncated unexpectedly.\n");
		else if (exc == IO_ERR_FILE_NOT_OPEN)
			fprintf(stderr, "Trying to access a file which is not open!(?)\n");
		return -1;
	}
	return 0;
}
