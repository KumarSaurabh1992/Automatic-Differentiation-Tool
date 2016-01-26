/*
 * Newtonbox.cpp
 *
 *  Created on: Mar 7, 2015
 *      Author: saurabh
 */

#include "Newtonbox.h"

void Newton_box::compute(){
	ADT<double> adt;
	int number_of_functions = 2;
	int number_of_variables = 2;
	SparseMatrix spm;
	spm.initialize(number_of_functions);
	x.set_size(number_of_functions);
	cout << "Enter the initial Guess" << endl;
	x.read_values();
	b.set_size(number_of_functions);
	Jacobian jcb;
	jcb.set_number_of_functions(number_of_functions);
	jcb.set_number_of_variables(number_of_variables);
	double err = 199;
	while(err > pow(10.0,-6)){
	jcb.createJacob(spm,b,x);
//	b.display();
//	cout << "****" <<endl;
//	x.display();
//	cout << "****" <<endl;
//	spm.display();
	Vector<double> y;
	y.set_size(number_of_functions);
	Iterative_Solver its;
	its.BiCG(spm,x,b);
	Vector<double> x1;
	x1.set_size(number_of_functions);
	x1.sub(x,y);
//	x.display();
	err = cal_error(x,x1);
//	x1.display();
	x.copy_values(x1);
	x.display();
	cout<<endl;
	}


}
//! This calculates the error based on the relative difference between two points
/*! @param X Vector X */
/*! @param X1 Vector X1*/

double Newton_box::cal_error(Vector<double> &X, Vector<double> &X1){
	double err = 0.0;
	for (int i = 0; i < X.get_size(); i++){
		if (fabs(X[i] - 0.0) > pow(10.0,-6)){
			if (err <  fabs((X[i] - X1[i])/X[i])){
				err = fabs((X[i] - X1[i])/X[i]);
			}

		}

	}
	//	cout << err << endl;
	return err;
}
