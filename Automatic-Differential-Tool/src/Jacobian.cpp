/*
 * Jacobian.cpp
 *
 *  Created on: Mar 4, 2015
 *      Author: saurabh
 */

#include "Jacobian.h"
#include "SparseMatrix.h"
//void Jacobian::set_size(int m, int n){
//	jacob.set_size(m,n);
//}

void Jacobian::set_number_of_functions(int N){
	number_of_functions = N;
}

void Jacobian::set_number_of_variables(int N){
	number_of_variables = N;
}
void Jacobian::createJacob(SparseMatrix & spm, Vector<double> & b, Vector<double> & pointval){
	//	array<float> jacob(number_of_functions,number_of_variables);

//	spm.initialize(number_of_functions);

	ADT<float> *AD = new ADT<float> [number_of_variables];

	ADT <float> adt1,adt2,adt3,adt4;
//	Vector<float> pointval(number_of_variables);
//	pointval.read_values();

	for(int i = 0; i < number_of_variables; i++){
		AD[i].ADconst(pointval[i]);
	}

	for(int j = 0; j < number_of_functions;j++){
		for(int i = 0; i < number_of_variables; i++){
			AD[i].ADvar(pointval[i]);
			if(j == 0){
				adt2.multiply(AD[1],2);
				adt1.substract(AD[0],adt2);
			}
			else{
				adt2.multiply(AD[0],2);
				adt3.multiply(AD[1],2);
				double z = 1;
				adt4.add(adt2,adt3);
				adt1.substract(adt4,z);
			}
			b[j] = adt1.get_fun_value();
			spm.update_matrix(adt1.get_first_derivative(),j,i);
			AD[i].ADconst(pointval[i]);
		}
	}
//	spm.display();

}


