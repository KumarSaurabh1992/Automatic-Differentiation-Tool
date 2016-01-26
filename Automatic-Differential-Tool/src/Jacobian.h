/*
 * Jacobian.h
 *
 *  Created on: Mar 4, 2015
 *      Author: saurabh
 */

#ifndef JACOBIAN_H_
#define JACOBIAN_H_
#include "array.h"
#include "ADT.h"
#include "Vector.h"
#include "SparseMatrix.h"
class Jacobian {

	int number_of_functions;
	int number_of_variables;
public:
	void set_size(int , int);
	void createJacob(SparseMatrix &, Vector<double> &, Vector<double> &);
	void set_number_of_functions(int );
	void set_number_of_variables(int );
	void create_grad(ADT<float> &,int );
};

#endif /* JACOBIAN_H_ */
