//============================================================================
// Name        : ADT.cpp
// Author      : saurabh
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "ADT.h"
#include "Gradient.h"
#include "Jacobian.h"
#include "Newtonbox.h"
#include "ComputeMatrix.h"
using namespace std;


int main() {
//	ADT<float> adt1,adt2,adt3,adt4;
//	float y =2.0;
//	adt1.ADvar(y);
//	adt2.sine(adt1);
//	adt3.exponential(adt1);
//	adt4.division(adt2,adt3);
//	cout << adt4.get_first_derivative() << endl;
//	cout << adt4.get_second_derivative();
//	Gradient grad;
//	grad.set_number_of_variables(2);
//	grad.create_gradient();
//	Jacobian jcb;
//	jcb.set_number_of_functions(2);
//	jcb.set_number_of_variables(2);
//	jcb.createJacob();
//	Newton_box nb;
//	nb.compute();
//	Vector<double> val(9);
//	val.initialize();
	ComputeMatrix cpm;
	cpm.solveNonLinear();
	return 0;
}
