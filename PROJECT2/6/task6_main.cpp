/*****************************************************************************/
// File: vertical_filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#define PI 3.1415926F

#include "io_bmp.h"
#include "aligned_image_comps.h"
#include <math.h>
#include <time.h>

double time_1 = 0;

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
/*                        my_aligned_image_comp::filter                      */
/*****************************************************************************/


void my_aligned_image_comp::filter_1(my_aligned_image_comp *in, float sigma)
{
	int FILTER_EXTENT = ceil(3 * sigma);
	int FILTER_TAPS = 2 * FILTER_EXTENT + 1;

	// Create the vertical filter PSF as a local array on the stack.
	float *filter_buf_1 = new float[FILTER_TAPS];

	for (int i = 0; i < FILTER_TAPS; ++i) {
		filter_buf_1[i] = 0.0F;
	}

	float *mirror_kernal_1 = filter_buf_1 + FILTER_EXTENT;

	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT;t++)
		mirror_kernal_1[t] = ((t*t) / (sigma*sigma) - 1)*exp(-(t*t) / (2 * sigma*sigma)) / (2 * PI*pow(sigma, 4));

	// Perform the convolution

	for (int r = 0; r < height; r++) {
		for (int c = 0; c < width; c++)
		{
			float *ip = in->buf + (r*in->stride) + c;
			float *op = this->buf + r*stride + c;
			float sum = 0.0F;
			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; ++y)
				sum += ip[y*in->stride] * mirror_kernal_1[y];
			*op = sum;
		}
	}

}
void my_aligned_image_comp::filter_2(my_aligned_image_comp *in, float sigma)
{
	int FILTER_EXTENT = ceil(3 * sigma);
	int FILTER_TAPS = 2 * FILTER_EXTENT + 1;

	// Create the vertical filter PSF as a local array on the stack.
	float *filter_buf_2 = new float[FILTER_TAPS];

	for (int i = 0; i < FILTER_TAPS; ++i) {
		filter_buf_2[i] = 0.0F;
	}

	float *mirror_kernal_2 = filter_buf_2 + FILTER_EXTENT;

	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT;t++)
		mirror_kernal_2[t] = exp(-(t*t) / (2 * sigma*sigma));

	// Perform the convolution

	for (int r = 0; r < height; r++) {
		for (int c = 0; c < width; c++)
		{
			float *ip = in->buf + (r*in->stride) + c;
			float *op = this->buf + r*stride + c;
			float sum = 0.0F;
			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; ++y)
				sum += ip[y*in->stride] * mirror_kernal_2[y];
			*op = sum;
		}
	}

}
void my_aligned_image_comp::out_sum_filter(my_aligned_image_comp *in1, my_aligned_image_comp *in2, float sigma, float alpha)
{
	int FILTER_EXTENT = ceil(3 * sigma);
	int FILTER_TAPS = 2 * FILTER_EXTENT + 1;

	// Create the vertical filter PSF as a local array on the stack.
	float *filter_buf_1 = new float[FILTER_TAPS];
	float *filter_buf_2 = new float[FILTER_TAPS];

	for (int i = 0; i < FILTER_TAPS; ++i) {
		filter_buf_1[i] = 0.0F;
		filter_buf_2[i] = 0.0F;
	}

	float *mirror_kernal_1 = filter_buf_1 + FILTER_EXTENT;
	float *mirror_kernal_2 = filter_buf_2 + FILTER_EXTENT;

	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT;t++)
		mirror_kernal_1[t] = ((t*t) / (sigma*sigma) - 1)*exp(-(t*t) / (2 * sigma*sigma)) / (2 * PI*pow(sigma, 4));

	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT;t++)
		mirror_kernal_2[t] = exp(-(t*t) / (2 * sigma*sigma));

	// Perform the convolution

	for (int r = 0; r < height; ++r) {
		for (int c = 0; c < width; c++) {
			float *ip1 = in1->buf + r*in1->stride + c;
			float *ip2 = in2->buf + r*in2->stride + c;
			float *op = buf + r*stride + c;
			float sum = 0.0F;
			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++)
				sum += ip1[y] * mirror_kernal_1[y] + ip2[y] * mirror_kernal_2[y];
			*op = 128 + sum * alpha;
		}
	}
}
void my_aligned_image_comp::zerocrossing(my_aligned_image_comp *in)
{
	for (int r = 0; r < height; ++r) {
		for (int c = 0; c < width; c++) {
			float *ip = in->buf + r*in->stride + c;
			float *op = buf + r*stride + c;

			if (ip[0] > 126 && ip[0] < 130)
				if (ip[0] > 128)
				{
					if (ip[1] < 128)
						*op = 255;
					else if (ip[-1] < 128)
						*op = 255;
					else if (ip[stride] < 128)
						*op = 255;
					else if (ip[-stride] < 128)
						*op = 255;
					else if (ip[-stride - 1] < 128)
						*op = 255;
					else if (ip[-stride + 1] < 128)
						*op = 255;
					else if (ip[stride - 1] < 128)
						*op = 255;
					else if (ip[stride + 1] < 128)
						*op = 255;
					else
						*op = 0;
				}
				else
					*op = 0;
		}
	}
}
void my_aligned_image_comp::dilation(my_aligned_image_comp *in) {
	int r, c;
	for (r = 0; r < height; r++) {
		for (c = 0; c < width; c++) {
			float *ip = in->buf + r * in->stride + c;
			float *op = buf + r * stride + c;
			if (ip[0] == 0) {
				if (ip[-1] == 255) {
					*op = 255;
				}

				else if (ip[1] == 255) {
					*op = 255;
				}

				else if (ip[in->stride] == 255) {
					*op = 255;
				}

				else if (ip[-in->stride] == 255) {
					*op = 255;
				}

				else {
					*op = 0;
				}

			}
			else if (ip[0] == 255) {
				*op = 255;
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
	clock_t start_time = clock();

	if (argc != 6)
	{
		fprintf(stderr, "Usage: %s <in bmp file> <out bmp file> <out dilated bmp file> <sigma> <alpha>\n", argv[0]);
		return -1;
	}

	int err_code = 0;
	try {
		// Read the input image
		bmp_in in;
		if ((err_code = bmp_in__open(&in, argv[1])) != 0)
			throw err_code;

		float sigma = atof(argv[4]);
		float alpha = atof(argv[5]);

		int width = in.cols, height = in.rows;
		int n, num_comps = in.num_components;

		my_aligned_image_comp *input_comps =
			new my_aligned_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			input_comps[n].init(height, width, 9);

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
		delete[] line;
		bmp_in__close(&in);

		my_aligned_image_comp *output_comps_0 =
			new my_aligned_image_comp[num_comps];
		my_aligned_image_comp *output_comps_1 =
			new my_aligned_image_comp[num_comps];
		my_aligned_image_comp *output_comps =
			new my_aligned_image_comp[num_comps];
		my_aligned_image_comp * intermediea_comps_1 =
			new my_aligned_image_comp[num_comps];
		my_aligned_image_comp * intermediea_comps_2 =
			new my_aligned_image_comp[num_comps];

		int out_width, out_height;
		out_width = width;
		out_height = height;

		for (n = 0; n < num_comps; n++) {
			intermediea_comps_1[n].init(height, width, 9);
			intermediea_comps_2[n].init(height, width, 9);
			output_comps_0[n].init(out_height, out_width, 1);
			output_comps_1[n].init(out_height, out_width, 1);
			output_comps[n].init(out_height, out_width, 1);
		}

		// Process the image, all in floating point (easy)
		for (n = 0; n < num_comps; n++) {
			input_comps[n].perform_boundary_extension();
			intermediea_comps_1[n].perform_boundary_extension();
			intermediea_comps_2[n].perform_boundary_extension();
		}

		printf("Start filtering!\r\n");

		for (n = 0; n < num_comps; n++)
		{
			intermediea_comps_1[n].filter_2(input_comps + n, sigma);
			intermediea_comps_2[n].filter_1(input_comps + n, sigma);
		}
		for (n = 0; n < num_comps; n++) {
			intermediea_comps_1[n].perform_boundary_extension();
			intermediea_comps_2[n].perform_boundary_extension();
			output_comps_0[n].out_sum_filter(intermediea_comps_1 + n, intermediea_comps_2 + n,sigma,alpha);
		}
		for (n = 0; n < num_comps; n++) {
			output_comps_0[n].perform_boundary_extension();
			output_comps_1[n].zerocrossing(output_comps_0 + n);
		}
		for (n = 0; n < num_comps; n++) {
			output_comps_0[n].perform_boundary_extension();
			output_comps[n].dilation(output_comps_1 + n);
		}
		printf("Filtering end!\r\n");


		// Write the image back out again
		bmp_out out;
		if ((err_code = bmp_out__open(&out, argv[3], out_width, out_height, num_comps)) != 0)
			throw err_code;

		io_byte *oline = new io_byte[out_width*num_comps];

		for (r = out_height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			for (n = 0; n < num_comps; n++)
			{
				io_byte *dst = oline + n; // Points to first sample of component n
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
			bmp_out__put_line(&out, oline);
		}
		if ((err_code = bmp_out__open(&out, argv[2], out_width, out_height, num_comps)) != 0)
			throw err_code;

		//io_byte *oline = new io_byte[out_width*num_comps];

		for (r = out_height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			for (n = 0; n < num_comps; n++)
			{
				io_byte *dst = oline + n; // Points to first sample of component n
				float *src = output_comps_1[n].buf + r * output_comps_1[n].stride;
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
			bmp_out__put_line(&out, oline);
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
