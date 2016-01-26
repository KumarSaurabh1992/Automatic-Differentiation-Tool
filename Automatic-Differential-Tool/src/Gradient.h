/*
 * Gradient.h
 *
 *  Created on: Mar 4, 2015
 *      Author: saurabh
 */

#ifndef GRADIENT_H_
#define GRADIENT_H_

#include "ADT.h"
#include "Vector.h"

class Gradient {
	int number_of_variables;

public:
	Vector<float> grad;
	void set_number_of_variables(int );
	void create_gradient();


};

#endif /* GRADIENT_H_ */
