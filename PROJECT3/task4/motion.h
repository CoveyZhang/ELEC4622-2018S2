// File: motion.h
// Author: David Taubman
// Last Revised: 30 September, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#ifndef MOTION_H

struct mvector {
    mvector() { x=y=0; }
    int x, y;
  };

struct global_vector
{
	global_vector() { 
		x = y = 0;
		x_int = int(x), y_int = int(y);
		x_float = x - x_int, y_float = y - y_int;
	}
	void calculate()
	{
		x_int = int(x), y_int = int(y);
		x_float = x - x_int, y_float = y - y_int;
	}
	float x, y;
	int x_int, y_int;
	float x_float, y_float;
};
#endif // MOTION_H
