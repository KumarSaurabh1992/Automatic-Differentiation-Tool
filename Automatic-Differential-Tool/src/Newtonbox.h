/*
 * Newtonbox.h
 *
 *  Created on: Mar 7, 2015
 *      Author: saurabh
 */

#ifndef NEWTONBOX_H_
#define NEWTONBOX_H_

#include "SparseMatrix.h"
#include "ADT.h"
#include "Jacobian.h"
#include "IterativeSolver.h"

class Newton_box {
	SparseMatrix spm;
	Vector<double> x;
	Vector<double> b;
public:
	void compute();
	double cal_error(Vector<double> &, Vector<double> &);
};

#endif /* NEWTONBOX_H_ */
