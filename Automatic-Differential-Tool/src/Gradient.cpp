/*
 * Gradient.cpp
 *
 *  Created on: Mar 4, 2015
 *      Author: saurabh
 */

#include "Gradient.h"


void Gradient::set_number_of_variables(int N){
	number_of_variables = N;
}

void Gradient::create_gradient(){
	ADT<float> *AD  = new ADT<float> [number_of_variables];
	grad.set_size(number_of_variables);
	Vector<float> pointval(number_of_variables);
	pointval.read_values();
	for(int i = 0; i < number_of_variables; i++){
		AD[i].ADconst(pointval[i]);
	}
	for (int i = 0; i < number_of_variables; i++){
		AD[i].ADvar(pointval[i]);
		ADT<float> adt1,adt2;
		adt2.sine(AD[0]);
		adt1.add(adt2,AD[1]);
		grad[i] = adt1.get_first_derivative();
		AD[i].ADconst(pointval[i]);
	}
	grad.display();
}

