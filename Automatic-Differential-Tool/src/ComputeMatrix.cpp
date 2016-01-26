/*
 * ComputeMatrix.cpp
 *
 *  Created on: Jan 30, 2015
 *      Author: saurabh
 */

#include "ComputeMatrix.h"
#include <fstream>
#include <math.h>
#include <time.h>
#include <fstream>

//! @return number of x - divisions
int ComputeMatrix::get_number_of_division_x(){
	return number_of_divisions_x;
}

//! @return number of y - divisions
int ComputeMatrix::get_number_of_division_y(){
	return number_of_divisions_y;
}

//! this sets the number of x - divisions and x-spacing
/*! @param N number of X - divisions*/
void ComputeMatrix::set_number_of_division_x(int N){
	h = 1.0/(N);
	number_of_divisions_x = N + 1;
}

//! This sets the number of y - divisions and y-spacing
/*! @param N number of y - divisions*/
void ComputeMatrix::set_number_of_division_y(int N){
	k = 1.0/(N);
	number_of_divisions_y = N + 1;
	//	cout <<"division y";
}

//! This takes the nput about the boundary condition
/*! Boundary condition is given in form of alpha*u + beeta*u' = gama*/
/*! The user inputs the boundary condition in the following order: bottom, right, top, left*/
void ComputeMatrix::get_boundary_condition(){
	int N = 4;
	set_boundary_size(N);
	double coeff_alpha,coeff_beeta,coeff_gama;
	//	alpha.initialize();
	//	beeta.initialize();
	//	gama.initialize();
	for(int i = 0; i < N; i++){
		//		cout << "Start" << endl;
		//		cout << alpha.get_size() << beeta.get_size() << gama.get_size() << endl;
		cout << "Enter the value of alpha" << endl;
		cin >> coeff_alpha;
		alpha[i] = coeff_alpha;
		//		alpha.display();
		cout << "Enter the value of beeta" << endl;
		cin >> coeff_beeta;
		beeta[i] = coeff_beeta;
		cout << "Enter the value of gamma" << endl;
		cin >> coeff_gama;
		//		cout << "End" << endl;
		//		cout << alpha.get_size() << beeta.get_size() << gama.get_size() << endl;
		gama[i] = coeff_gama;
		//		cout << "Done";
	}

}
//! Allocates the space for the boundary vectors
/*! @param N number of boundaries for the domain*/
void ComputeMatrix::set_boundary_size(int N){
	alpha.set_size(N);
	beeta.set_size(N);
	gama.set_size(N);
}


//! It creates the vector containing the information about spacing in the x direction
void ComputeMatrix::compute_h(){

	cout << "Enter the number of spacing in the x direction" << endl;
	cin >> number_of_x;
	x_spacing = new double*[number_of_x];
	for(int i = 0; i < number_of_x;i++){
		x_spacing[i] = new double[5];
	}
	double temp;
	for(int i = 0; i < number_of_x; i++){

		cout << "Enter the starting coordinate point" << endl;
		cin >> temp;
		x_spacing[i][0] = temp;
		cout << "Enter the end coordinate point" << endl;
		cin >> temp;
		x_spacing[i][1] = temp;
		cout << "Enter the value for h" << endl;
		cin >> temp;
		x_spacing[i][2] = temp;
		temp = (fabs(x_spacing[i][0] - x_spacing[i][1]))/(x_spacing[i][2]);
		//		cout << round(temp);
		//		cout << temp << endl;;
		//		cout << (fabs(temp - round(temp)));
		if ((fabs(temp - round(temp)) < pow(10.0,-6))) {
			x_spacing[i][3] = round(temp);
			x_spacing[i][4] = x_spacing[i][2];
			//			cout << "Entered";
		}
		else{
			x_spacing[i][3] = ceil(temp);
			x_spacing[i][4] = x_spacing[i][1] - (x_spacing[i][0]+ floor(temp)*x_spacing[i][2]);
		}
	}
	int m = 0;
	for (int i = 0; i < number_of_x ;i++){

		m += x_spacing[i][3];
		x_spacing[i][3] = m;
	}
	x_spacing[number_of_x - 1][3] +=1;
	//


}

//! It creates the vector containing information about spacing in the y - direction
void ComputeMatrix::compute_k(){
	cout << "Enter the number of spacing in the y direction" << endl;
	cin >> number_of_y;
	y_spacing = new double*[number_of_y];
	for(int i = 0; i < number_of_y;i++){
		y_spacing[i] = new double[5];
	}
	double temp;
	for(int i = 0; i < number_of_y; i++){

		cout << "Enter the starting coordinate point" << endl;
		cin >> temp;
		y_spacing[i][0] = temp;
		cout << "Enter the end coordinate point" << endl;
		cin >> temp;
		y_spacing[i][1] = temp;
		cout << "Enter the value for k" << endl;
		cin >> temp;
		y_spacing[i][2] = temp;
		temp = (fabs(y_spacing[i][0] - y_spacing[i][1]))/(y_spacing[i][2]);
		//		cout << round(temp);
		//		cout << temp << endl;;
		//		cout << (fabs(temp - round(temp)));
		if ((fabs(temp - round(temp)) < pow(10.0,-6))) {
			y_spacing[i][3] = round(temp);
			y_spacing[i][4] = y_spacing[i][2];
			//			cout << "Entered";
		}
		else{
			y_spacing[i][3] = ceil(temp);
			y_spacing[i][4] = y_spacing[i][1] - (y_spacing[i][0]+ floor(temp)*y_spacing[i][2]);
		}
	}
	int m = 0;
	for (int i = 0; i < number_of_y ;i++){

		m += y_spacing[i][3];
		y_spacing[i][3] = m;
	}
	y_spacing[number_of_y - 1][3] +=1;
	//

}

//! At a particular point, it returns the spacing in the left direction
/*! @param x_division the point of our interest*/
/*! @return the left x-spacing of that point*/
double ComputeMatrix::get_h_before(int x_division){
	if (x_division == number_of_divisions_x - 1){
		return x_spacing[number_of_x - 1][4];
	}
	for(int i = 0; i < number_of_x;i++){
		if(x_division == x_spacing[i][3]){
			return x_spacing[i][4];
		}
		else if (x_division < x_spacing[i][3]){
			return x_spacing[i][2];
		}
	}
	return 0.0;
}
//! At a particular point, it returns the spacing in the right direction
/*! @param x_division the point of our interest*/
/*! @return the right x-spacing of that point*/
double ComputeMatrix::get_h_after(int x_division){

	if (x_division == number_of_divisions_x - 2){
		return x_spacing[number_of_x - 1][4];
	}
	for(int i = 0; i < number_of_x;i++){
		if(x_division == x_spacing[i][3] - 1){
			return x_spacing[i][4];
		}
		else if (x_division < x_spacing[i][3]){
			return x_spacing[i][2];
		}
	}
	return 0.0;
}
//! At a particular point, it returns the spacing in the top direction
/*! @param y_division the point of our interest*/
/*! @return the top y-spacing of that point*/
double ComputeMatrix::get_k_after(int y_division){

	if (y_division == number_of_divisions_y - 2){
		return y_spacing[number_of_y - 1][4];
	}
	for(int i = 0; i < number_of_y;i++){
		if(y_division == y_spacing[i][3] - 1){
			return y_spacing[i][4];
		}
		else if (y_division < y_spacing[i][3]){
			return y_spacing[i][2];
		}
	}
	return 0.0;
}
//! At a particular point, it returns the spacing in the bottom direction
/*! @param y_division the point of our interest*/
/*! @return the bottom y-spacing of that point*/
double ComputeMatrix::get_k_before(int y_division){
	if (y_division == number_of_divisions_y - 1){
		return y_spacing[number_of_y - 1][4];
	}
	for(int i = 0; i < number_of_y;i++){
		if(y_division == y_spacing[i][3]){
			return y_spacing[i][4];
		}
		else if (y_division < y_spacing[i][3]){
			return y_spacing[i][2];
		}
	}
	return 0.0;
}

void ComputeMatrix::make_variable_f1(double value){
	adt_f1.ADvar(1.0);
}

void ComputeMatrix::make_constant_f1(double value){
	adt_f1.ADconst(1.0);
}

void ComputeMatrix::make_variable_f2(double value){
	adt_f2.ADvar(1.0);
}

void ComputeMatrix::make_constant_f2(double value){
	adt_f2.ADconst(1.0);
}

void ComputeMatrix::make_variable_f3(double value){

	adt_f3.ADvar(2.0);
}

void ComputeMatrix::make_constant_f3(double value){
	adt_f3.ADconst(2.0);
}

void ComputeMatrix::make_variable_f4(double value){
	adt_f4.ADvar(1.0);
}

void ComputeMatrix::make_constant_f4(double value){
	adt_f4.ADconst(1.0);
}

void ComputeMatrix::make_variable_f5(double value){
	adt_f5.ADvar(1.0);
}

void ComputeMatrix::make_constant_f5(double value){
	adt_f5.ADconst(1.0);
}

//! This is used to compute the matrix A and vector B when the grid is non - uniform in both x and y coordinates
/*! Matrix A is stored in sparse form.*/
/*! Vector X is initialized to 0*/
void ComputeMatrix::compute(Vector<double> & val){
	double temp_center, temp_left, temp_right, temp_bottom, temp_top,h_before,h_after,k_before,k_after;
	//	ADT<double> adt_f1, adt_f2, adt_f3;
	//	adt_f1.ADconst(1.0);
//	//	adt_f2.ADconst(1.0);
//	//	adt_f3.ADconst(1.0);
	ADT <double> adt_bottom, adt_left, adt_center, adt_right, adt_top,adt_x_deriv,adt_y_deriv;
	ADT<double> adt_rhs;
	ADT<double> adt_temp1,adt_temp2,adt_temp3,adt_temp4;
	ADT<double> adt_fun;
//
	for(int i = 0; i < number_of_divisions_y; i++){
		//		cout << i << endl;

		for(int j = 0; j < number_of_divisions_x; j++){

			int row = (i)*(number_of_divisions_x) + (j);
			int col_left = (i)*(number_of_divisions_x) + (j - 1);
			int col_right = i*(number_of_divisions_x) + (j + 1);
			int col_top = (i + 1)*(number_of_divisions_x) + j;
			int col_bottom = (i - 1)*(number_of_divisions_x) + j;
			int col_center = i*(number_of_divisions_x) + j;
			//
			if ((j == 0) && (i == 0)){									// bottom left point
				if((fabs(beeta[0] - 0.0) < pow(10.0,-6)) || (fabs(beeta[3] - 0.0) < pow(10.0,-6))){
					if (fabs(beeta[0] - 0.0) < pow(10.0,-6)){
						temp_center = alpha[0];
						spm.update_matrix(temp_center,row, col_center);
						B[row] = alpha[0]*val[row] - gama[0];
					}
					else{
						temp_center = alpha[3];
						spm.update_matrix(temp_center,row, col_center);
						B[row] = alpha[3]*val[row]-gama[3];
					}
				}
				else{

					h_after = get_h_after(j);
					k_after = get_k_after(i);
					double temp_center1 = -((2.0/(h_after*h_after))*(1 - (alpha[3]*h_after/beeta[3])));
					double temp_center2 = -((2.0/(k_after*k_after))*(1 - (alpha[0]*k_after/beeta[0])));

					double temp_right = 2.0/(h_after*h_after);
					double temp_top = 2.0/(k_after*k_after);
					double temp_rhs1 = -2.0*gama[3]/(beeta[3]*h_after);
					double temp_rhs2 = -2.0*gama[0]/(beeta[0]*k_after);
					//					double temp_rhs = -(rhs + 2.0*gama[3]/(beeta[3]*h_after) + (2.0*gama[0]/(beeta[0]*k_after)));
					//					temp_center = -(((2.0/(h_after*h_after))*(1 - (alpha[3]*h_after/beeta[3]))) + ((2.0/(k_after*k_after))*(1 - (alpha[0]*k_after/beeta[0]))));

					double temp_deriv_right1 = -2*h_after*alpha[3]/beeta[3];
					double temp_deriv_right2 = 2*h_after*gama[3]/beeta[3];
					double temp_deriv_top1 = -2*k_after*alpha[0]/beeta[0];
					double temp_deriv_top2 = 2*k_after*gama[0]/beeta[0];

					surrounding = new ADT<double>[3];
					surrounding[0].ADconst(val[col_center]);
					surrounding[1].ADconst(val[col_right]);
					surrounding[2].ADconst(val[col_top]);

//					spm.display();

					for(int l = 0; l < 3; l++){
						make_constant_f1(val[col_center]);
						make_constant_f2(val[col_center]);
						make_constant_f3(val[col_center]);
						make_constant_f4(val[col_center]);
						make_constant_f5(val[col_center]);
						if(l == 0){
							make_variable_f1(val[col_center]);
							make_variable_f2(val[col_center]);
							make_variable_f3(val[col_center]);
							make_variable_f4(val[col_center]);
							make_variable_f5(val[col_center]);
						}
						surrounding[l].ADvar();
						adt_temp1.multiply(surrounding[0],temp_center1);
						adt_temp2.multiply(surrounding[0],temp_center2);
						adt_temp1.multiply(adt_f1);
						adt_temp2.multiply(adt_f2);
						adt_center.add(adt_temp1,adt_temp2);
						//						surrounding[0].add(adt_temp1,adt_temp2);
						adt_temp1.multiply(surrounding[1],temp_right);
						adt_right.multiply(adt_temp1,adt_f1);
						//						surrounding[1].multiply(temp_right);
						//						surrounding[1].multiply(adt_f1);
						adt_temp1.multiply(surrounding[2],temp_top);
						adt_top.multiply(adt_temp1,adt_f2);
						//						surrounding[2].multiply(temp_top);
						//						surrounding[2].multiply(adt_f2);

						adt_temp1.multiply(surrounding[0],temp_deriv_right1);
						adt_temp2.multiply(adt_temp1,adt_f4);
						adt_temp3.multiply(adt_f4,temp_deriv_right2);
						adt_x_deriv.add(adt_temp2,adt_temp3);

						adt_temp1.multiply(surrounding[0],temp_deriv_top1);
						adt_temp2.multiply(adt_temp1,adt_f5);
						adt_temp3.multiply(adt_f5,temp_deriv_top2);
						adt_y_deriv.add(adt_temp2,adt_temp3);

						adt_temp1.multiply(adt_f1,temp_rhs1);
						adt_temp2.multiply(adt_f2,temp_rhs2);
						adt_rhs.add(adt_temp1,adt_temp2);
						adt_rhs.add(adt_f3);

						adt_temp1.add(adt_center,adt_right);
						adt_temp2.add(adt_top,adt_rhs);
						adt_temp3.add(adt_x_deriv,adt_y_deriv);

						adt_fun.add(adt_temp1,adt_temp2);
						adt_fun.add(adt_temp3);

						double val = adt_fun.get_first_derivative();

						if(l == 0){
							spm.update_matrix(val,row,col_center);
						}
						else if(l == 1){
							spm.update_matrix(val,row,col_right);
						}
						else if(l == 2){
							spm.update_matrix(val,row,col_top);
						}
						surrounding[l].ADConst();
					}
					B[row] = adt_fun.get_fun_value();
					delete surrounding;
				}

			}
			else if ((i == number_of_divisions_y - 1) && (j == 0)){			// top left point
				if((fabs(beeta[2] - 0.0) < pow(10.0,-6)) || (fabs(beeta[3] - 0.0) < pow(10.0,-6))){
					if (fabs(beeta[2] - 0.0) < pow(10.0,-6)){
						temp_center = alpha[2];
						spm.update_matrix(temp_center,row, col_center);
						B[row] = alpha[2]*val[row]-gama[2];
					}
					else{
						temp_center = alpha[3];
						spm.update_matrix(temp_center,row, col_center);
						B[row] = alpha[3]*val[row] - gama[3];
					}
				}
				else{
					h_after = get_h_after(j);
					k_before = get_k_before(i);

					double temp_center1 = -(((2.0/(h_after*h_after))*(1 - (alpha[3]*h_after/beeta[3]))));
					double temp_center2 = -(((2.0/(k_before*k_before))*(1 + (alpha[2]*k_before/beeta[2]))));
					//					temp_center = -(((2.0/(h_after*h_after))*(1 - (alpha[3]*h_after/beeta[3]))) + ((2.0/(k_before*k_before))*(1 + (alpha[2]*k_before/beeta[2]))));
					temp_right = 2.0/(h_after*h_after);
					temp_bottom = 2.0/(k_before*k_before);
					double temp_rhs1 = - 2.0*gama[3]/(beeta[3]*h_after);
					double temp_rhs2 = (2.0*gama[2]/(beeta[2]*k_before));

					double temp_deriv_right1 = -2*h_after*alpha[3]/beeta[3];
					double temp_deriv_right2 = 2*h_after*gama[3]/beeta[3];
					double temp_deriv_top1 = -2*k_before*alpha[2]/beeta[2];
					double temp_deriv_top2 = 2*k_before*gama[2]/beeta[2];

//					ADT<double> surrounding[3];
					surrounding = new ADT<double>[3];
					surrounding[0].ADconst(val[col_bottom]);
					surrounding[1].ADconst(val[col_center]);
					surrounding[2].ADconst(val[col_right]);

					for(int l = 0; l < 3; l++){
						make_constant_f1(val[col_center]);
						make_constant_f2(val[col_center]);
						make_constant_f3(val[col_center]);
						make_constant_f4(val[col_center]);
						make_constant_f5(val[col_center]);
						if(l == 1){
							make_variable_f1(val[col_center]);
							make_variable_f2(val[col_center]);
							make_variable_f3(val[col_center]);
							make_variable_f4(val[col_center]);
							make_variable_f5(val[col_center]);
						}
						surrounding[l].ADvar();
						adt_temp1.multiply(surrounding[0],temp_bottom);
						adt_bottom.multiply(adt_temp1,adt_f2);
						//						surrounding[0].multiply(temp_bottom);
						//						surrounding[0].multiply(adt_f2);
						adt_temp1.multiply(surrounding[1],temp_center1);
						adt_temp2.multiply(surrounding[1],temp_center2);
						adt_temp1.multiply(adt_f1);
						adt_temp2.multiply(adt_f2);
						adt_center.add(adt_temp1,adt_temp2);
						//						surrounding[1].add(adt_temp1,adt_temp2);
						adt_temp1.multiply(surrounding[2],temp_right);
						adt_right.multiply(adt_temp1,adt_f1);
						//						surrounding[2].multiply(temp_right);
						//						surrounding[2].multiply(adt_f1);
						adt_temp1.multiply(surrounding[1],temp_deriv_right1);
						adt_temp2.multiply(adt_temp1,adt_f4);
						adt_temp3.multiply(adt_f4,temp_deriv_right2);
						adt_x_deriv.add(adt_temp2,adt_temp3);

						adt_temp1.multiply(surrounding[1],temp_deriv_top1);
						adt_temp2.multiply(adt_temp1,adt_f5);
						adt_temp3.multiply(adt_f5,temp_deriv_top2);
						adt_y_deriv.add(adt_temp2,adt_temp3);

						adt_temp1.multiply(adt_f1,temp_rhs1);
						adt_temp2.multiply(adt_f2,temp_rhs2);
						adt_rhs.add(adt_temp1,adt_temp2);
						adt_rhs.add(adt_f3);

						adt_temp1.add(adt_bottom,adt_center);
						adt_temp2.add(adt_right,adt_rhs);
						adt_temp3.add(adt_x_deriv,adt_y_deriv);

						adt_fun.add(adt_temp1,adt_temp2);
						adt_fun.add(adt_temp3);

						double val = adt_fun.get_first_derivative();
						if(l == 0){
							spm.update_matrix(val,row,col_bottom);
						}
						else if(l == 1){
							spm.update_matrix(val,row,col_center);
						}
						else if(l == 2){
							spm.update_matrix(val,row,col_right);
						}
						surrounding[l].ADConst();
					}
					B[row] = adt_fun.get_fun_value();
					delete surrounding;
				}
			}
			else if ((i == 0) && (j == number_of_divisions_x - 1)){			// bottom right point
				if((fabs(beeta[0] - 0.0) < pow(10.0,-6)) || (fabs(beeta[1] - 0.0) < pow(10.0,-6))){
					if (fabs(beeta[0] - 0.0) < pow(10.0,-6)){
						temp_center = alpha[0];
						spm.update_matrix(temp_center,row, col_center);
						B[row] = alpha[0]*val[row] - gama[0];
					}
					else{
						temp_center = alpha[1];
						spm.update_matrix(temp_center,row, col_center);
						B[row] = alpha[1]*val[row] - gama[1];
					}
				}
				else{

					h_before = get_h_before(j);
					k_after = get_k_after(i);

					double temp_center1 = -(((2.0/(h_before*h_before))*(1 + (alpha[1]*h_before/beeta[1]))));
					double temp_center2 =  -(((2.0/(k_after*k_after))*(1 - (alpha[0]*k_after/beeta[0]))));
					//					temp_center = -(((2.0/(h_before*h_before))*(1 + (alpha[1]*h_before/beeta[1]))) + ((2.0/(k_after*k_after))*(1 - (alpha[0]*k_after/beeta[0]))));
					temp_left = 2.0/(h_before*h_before);
					temp_top = 2.0/(k_after*k_after);

					double temp_rhs1 = 2.0*gama[1]/(beeta[1]*h_before);
					double temp_rhs2 = -(2.0*gama[0]/(beeta[0]*k_after));

					double temp_deriv_right1 = -2*h_before*alpha[1]/beeta[1];
					double temp_deriv_right2 = 2*h_before*gama[1]/beeta[1];
					double temp_deriv_top1 = -2*k_after*alpha[0]/beeta[0];
					double temp_deriv_top2 = 2*k_after*gama[0]/beeta[0];

					surrounding = new ADT<double>[3];
//					ADT<double> surrounding[3];
					surrounding[0].ADconst(val[col_left]);
					surrounding[1].ADconst(val[col_center]);
					surrounding[2].ADconst(val[col_top]);

					for(int l = 0; l < 3; l++){
						make_constant_f1(val[col_center]);
						make_constant_f2(val[col_center]);
						make_constant_f3(val[col_center]);
						make_constant_f4(val[col_center]);
						make_constant_f5(val[col_center]);
						if(l == 1){
							make_variable_f1(val[col_center]);
							make_variable_f2(val[col_center]);
							make_variable_f3(val[col_center]);
							make_variable_f4(val[col_center]);
							make_variable_f5(val[col_center]);
						}
						surrounding[l].ADvar();
						adt_temp1.multiply(surrounding[0],temp_left);
						adt_left.multiply(adt_temp1,adt_f1);
						//						surrounding[0].multiply(temp_left);
						//						surrounding[0].multiply(adt_f1);
						adt_temp1.multiply(surrounding[1],temp_center1);
						adt_temp2.multiply(surrounding[1],temp_center2);
						adt_temp1.multiply(adt_f1);
						adt_temp2.multiply(adt_f2);
						adt_center.add(adt_temp1,adt_temp2);
						//						surrounding[1].add(adt_temp1,adt_temp2);


						adt_temp1.multiply(surrounding[2],temp_top);
						adt_top.multiply(adt_temp1,adt_f2);
						//						surrounding[2].multiply(temp_top);
						//						surrounding[2].multiply(adt_f2);

						adt_temp1.multiply(surrounding[1],temp_deriv_right1);
						adt_temp2.multiply(adt_temp1,adt_f4);
						adt_temp3.multiply(adt_f4,temp_deriv_right2);
						adt_x_deriv.add(adt_temp2,adt_temp3);

						adt_temp1.multiply(surrounding[1],temp_deriv_top1);
						adt_temp2.multiply(adt_temp1,adt_f5);
						adt_temp3.multiply(adt_f5,temp_deriv_top2);
						adt_y_deriv.add(adt_temp2,adt_temp3);

						adt_temp1.multiply(adt_f1,temp_rhs1);
						adt_temp2.multiply(adt_f2,temp_rhs2);
						adt_rhs.add(adt_temp1,adt_temp2);
						adt_rhs.add(adt_f3);

						adt_temp1.add(adt_left,adt_center);
						adt_temp2.add(adt_top,adt_rhs);
						adt_temp3.add(adt_x_deriv,adt_y_deriv);

						adt_fun.add(adt_temp1,adt_temp2);
						adt_fun.add(adt_temp3);

						double val = adt_fun.get_first_derivative();
						if(l == 0){
							spm.update_matrix(val,row,col_left);
						}
						else if(l == 1){
							spm.update_matrix(val,row,col_center);
						}
						else if(l == 2){
							spm.update_matrix(val,row,col_top);
						}
						surrounding[l].ADConst();
					}
					B[row] = adt_fun.get_fun_value();
					delete surrounding;


				}
			}
			else if ((i == number_of_divisions_y - 1) && (j == number_of_divisions_x - 1)){		// top right point
				if((fabs(beeta[2] - 0.0) < pow(10.0,-6)) || (fabs(beeta[1] - 0.0) < pow(10.0,-6))){
					if (fabs(beeta[2] - 0.0) < pow(10.0,-6)){
						temp_center = alpha[2];
						spm.update_matrix(temp_center,row, col_center);
						B[row] = alpha[2]*val[row] - gama[2];
					}
					else{
						temp_center = alpha[1];
						spm.update_matrix(temp_center,row, col_center);
						B[row] = alpha[1]*val[row] - gama[1];
					}
				}
				else{
					h_before = get_h_before(j);
					k_before = get_k_before(i);
					double temp_center1 = -(((2.0/(h_before*h_before))*(1 + (alpha[1]*h_before/beeta[1]))));
					double temp_center2 = - (((2.0/(k_before*k_before))*(1 + (alpha[2]*k_before/beeta[2]))));
					//					temp_center = -(((2.0/(h_before*h_before))*(1 + (alpha[1]*h_before/beeta[1]))) + ((2.0/(k_before*k_before))*(1 + (alpha[2]*k_before/beeta[2]))));
					temp_left = 2.0/(h_before*h_before);
					temp_bottom = 2.0/(k_before*k_before);

					double temp_rhs1 = 2.0*gama[1]/(beeta[1]*h_before);
					double temp_rhs2 = (2.0*gama[2]/(beeta[2]*k_before));

					double temp_deriv_right1 = -2*h_before*alpha[1]/beeta[1];
					double temp_deriv_right2 = 2*h_before*gama[1]/beeta[1];
					double temp_deriv_top1 = -2*k_before*alpha[2]/beeta[2];
					double temp_deriv_top2 = 2*k_before*gama[2]/beeta[2];

//					ADT<double> surrounding[3];
					surrounding = new ADT<double>[3];
					surrounding[0].ADconst(val[col_bottom]);
					surrounding[1].ADconst(val[col_left]);
					surrounding[2].ADconst(val[col_center]);

					for(int l = 0; l < 3; l++){
						make_constant_f1(val[col_center]);
						make_constant_f2(val[col_center]);
						make_constant_f3(val[col_center]);
						make_constant_f4(val[col_center]);
						make_constant_f5(val[col_center]);

						if(l == 2){
							make_variable_f1(val[col_center]);
							make_variable_f2(val[col_center]);
							make_variable_f3(val[col_center]);
							make_variable_f4(val[col_center]);
							make_variable_f5(val[col_center]);
						}
						surrounding[l].ADvar();
						adt_temp1.multiply(surrounding[0],temp_bottom);
						adt_bottom.multiply(adt_temp1,adt_f2);
						//						surrounding[0].multiply(temp_bottom);
						//						surrounding[0].multiply(adt_f2);
						adt_temp1.multiply(surrounding[1],temp_left);
						adt_left.multiply(adt_temp1,adt_f1);
						//						surrounding[1].multiply(temp_right);
						//						surrounding[1].multiply(adt_f1);
						adt_temp1.multiply(surrounding[2],temp_center1);
						adt_temp2.multiply(surrounding[2],temp_center2);
						adt_temp1.multiply(adt_f1);
						adt_temp2.multiply(adt_f2);
						adt_center.add(adt_temp1,adt_temp2);

						adt_temp1.multiply(surrounding[2],temp_deriv_right1);
						adt_temp2.multiply(adt_temp1,adt_f4);
						adt_temp3.multiply(adt_f4,temp_deriv_right2);
						adt_x_deriv.add(adt_temp2,adt_temp3);

						adt_temp1.multiply(surrounding[2],temp_deriv_top1);
						adt_temp2.multiply(adt_temp1,adt_f5);
						adt_temp3.multiply(adt_f5,temp_deriv_top2);
						adt_y_deriv.add(adt_temp2,adt_temp3);

						//						surrounding[2].add(adt_temp1,adt_temp2);
						adt_temp1.multiply(adt_f1,temp_rhs1);
						adt_temp2.multiply(adt_f2,temp_rhs2);
						adt_rhs.add(adt_temp1,adt_temp2);
						adt_rhs.add(adt_f3);

						adt_temp1.add(adt_bottom,adt_left);
						adt_temp2.add(adt_center,adt_rhs);
						adt_temp3.add(adt_x_deriv,adt_y_deriv);

						adt_fun.add(adt_temp1,adt_temp2);
						adt_fun.add(adt_temp3);

						double val = adt_fun.get_first_derivative();
						if(l == 0){
							spm.update_matrix(val,row,col_bottom);
						}
						else if(l == 1){
							spm.update_matrix(val,row,col_left);
						}
						else if(l == 2){
							spm.update_matrix(val,row,col_center);
						}
						surrounding[l].ADConst();
					}
					B[row] = adt_fun.get_fun_value();
					delete surrounding;
				}
			}
			else if ((i == 0)){												// bottom face
				if((fabs(beeta[0] - 0.0) < pow(10.0,-6))){
					temp_center = alpha[0];
					spm.update_matrix(temp_center,row, col_center);
					B[row] = alpha[0]*val[row] - gama[0];
				}
				else{
					h_before = get_h_before(j);
					h_after = get_h_after(j);
					k_after = get_k_after(i);

					double temp_center1 =  -((2.0/(h_after*h_before)));
					double temp_center2 = -(((2.0/(k_after*k_after))*(1 - (alpha[0]*k_after/beeta[0]))));
					//					temp_center = -((2.0/(h_after*h_before)) + ((2.0/(k_after*k_after))*(1 - (alpha[0]*k_after/beeta[0]))));

					double temp_deriv_top1 = -2*k_after*alpha[0]/beeta[0];
					double temp_deriv_top2 = 2*k_after*gama[0]/beeta[0];

					temp_left = 2.0/(h_before*(h_before + h_after));
					temp_top = 2.0/(k_after*k_after);
					temp_right = 2.0/(h_after*(h_before + h_after));


					double temp_rhs1 = - (2.0*gama[0]/(beeta[0]*k_after));
//					ADT<double> surrounding[4];
					surrounding = new ADT<double>[4];
					surrounding[0].ADconst(val[col_left]);
					surrounding[1].ADconst(val[col_center]);
					surrounding[2].ADconst(val[col_right]);
					surrounding[3].ADconst(val[col_top]);

					for(int l = 0; l < 4; l++){
						make_constant_f1(val[col_center]);
						make_constant_f2(val[col_center]);
						make_constant_f3(val[col_center]);
						make_constant_f4(val[col_center]);
						make_constant_f5(val[col_center]);
						if(l == 1){
							make_variable_f1(val[col_center]);
							make_variable_f2(val[col_center]);
							make_variable_f3(val[col_center]);
							make_variable_f4(val[col_center]);
							make_variable_f5(val[col_center]);
						}
						surrounding[l].ADvar();
						adt_temp1.multiply(surrounding[0],temp_left);
						adt_left.multiply(adt_temp1,adt_f1);
						//						surrounding[0].multiply(temp_left);
						//						surrounding[0].multiply(adt_f1);
						adt_temp1.multiply(surrounding[1],temp_center1);
						adt_temp2.multiply(surrounding[1],temp_center2);
						adt_temp1.multiply(adt_f1);
						adt_temp2.multiply(adt_f2);
						adt_center.add(adt_temp1,adt_temp2);

						//						surrounding[1].add(adt_temp1,adt_temp2);
						adt_temp1.multiply(surrounding[2],temp_right);
						adt_right.multiply(adt_temp1,adt_f1);
						//						surrounding[2].multiply(temp_right);
						//						surrounding[2].multiply(adt_f1);
						adt_temp1.multiply(surrounding[3],temp_top);
						adt_top.multiply(adt_temp1,adt_f2);

						adt_temp1.multiply(adt_f4,surrounding[2]);
						adt_temp2.multiply(adt_f4,surrounding[0]);
						adt_x_deriv.substract(adt_temp1,adt_temp2);

						adt_temp1.multiply(surrounding[1],temp_deriv_top1);
						adt_temp2.multiply(adt_temp1,adt_f5);
						adt_temp3.multiply(adt_f5,temp_deriv_top2);
						adt_y_deriv.add(adt_temp2,adt_temp3);
						//						surrounding[3].multiply(temp_top);
						//						surrounding[3].multiply(adt_f2);
						adt_temp1.multiply(adt_f2,temp_rhs1);
						adt_rhs.add(adt_f3,adt_temp1);
						adt_temp1.add(adt_left,adt_center);
						adt_temp2.add(adt_right,adt_top);
						adt_temp2.add(adt_rhs);
						adt_temp3.add(adt_x_deriv,adt_y_deriv);

						adt_fun.add(adt_temp1,adt_temp2);
						adt_fun.add(adt_temp3);

						double val = adt_fun.get_first_derivative();
						if(l == 0){
							spm.update_matrix(val,row,col_left);
						}
						else if(l == 1){
							//							cout << adt_rhs.get_first_derivative() << "*********" << endl;
							spm.update_matrix(val,row,col_center);

						}
						else if(l == 2){
							spm.update_matrix(val,row,col_right);
						}
						else if(l == 3){
							spm.update_matrix(val,row,col_top);
						}
						surrounding[l].ADConst();
					}

					B[row] = adt_fun.get_fun_value();
					delete surrounding;
				}
			}
			else if (i == number_of_divisions_y - 1){							// top face
				if((fabs(beeta[2] - 0.0) < pow(10.0,-6))){
					temp_center = alpha[2];
					spm.update_matrix(temp_center,row, col_center);
					B[row] = alpha[2]*val[row] - gama[2];
				}
				else{
					k_before = get_k_before(i);
					h_after = get_h_after(j);
					h_before = get_h_before(j);
					double temp_center1 = -((2.0/(h_before*h_after)));
					double temp_center2 = -((((2.0/(k_before*k_before))*(1 + (alpha[2]*k_before/beeta[2])))));
					//					temp_center = -((2.0/(h_before*h_after)) + ((2.0/(k_before*k_before))*(1 + (alpha[2]*k_before/beeta[2]))));

					double temp_deriv_top1 = -2*k_before*alpha[2]/beeta[2];
					double temp_deriv_top2 = 2*k_before*gama[2]/beeta[2];

					temp_left = 2.0/(h_before*(h_after + h_before));
					temp_bottom = 2.0/(k_after*k_after);
					temp_right = 2.0/(h_after*(h_after + h_before));
					double temp_rhs1 = (2.0*gama[2]/(beeta[2]*k_before));


//					ADT<double> surrounding[4];
					surrounding = new ADT<double>[4];
					surrounding[0].ADconst(val[col_bottom]);
					surrounding[1].ADconst(val[col_left]);
					surrounding[2].ADconst(val[col_center]);
					surrounding[3].ADconst(val[col_right]);

					for(int l = 0; l < 4; l++){
						make_constant_f1(val[col_center]);
						make_constant_f2(val[col_center]);
						make_constant_f3(val[col_center]);
						make_constant_f4(val[col_center]);
						make_constant_f5(val[col_center]);
						if(l == 2){
							make_variable_f1(val[col_center]);
							make_variable_f2(val[col_center]);
							make_variable_f3(val[col_center]);
							make_variable_f4(val[col_center]);
							make_variable_f5(val[col_center]);
						}
						surrounding[l].ADvar();
						adt_temp1.multiply(surrounding[0],temp_bottom);
						adt_bottom.multiply(adt_temp1,adt_f2);
						//						surrounding[0].multiply(temp_bottom);
						//						surrounding[0].multiply(adt_f2);
						adt_temp1.multiply(surrounding[1],temp_left);
						adt_left.multiply(adt_temp1,adt_f1);
						//						surrounding[1].multiply(temp_left);
						//						surrounding[1].multiply(adt_f1);
						adt_temp1.multiply(surrounding[2],temp_center1);
						adt_temp2.multiply(surrounding[2],temp_center2);
						adt_temp1.multiply(adt_f1);
						adt_temp2.multiply(adt_f2);
						adt_center.add(adt_temp1,adt_temp2);

						adt_temp1.multiply(adt_f4,surrounding[3]);
						adt_temp2.multiply(adt_f4,surrounding[1]);
						adt_x_deriv.substract(adt_temp1,adt_temp2);

						adt_temp1.multiply(surrounding[2],temp_deriv_top1);
						adt_temp2.multiply(adt_temp1,adt_f5);
						adt_temp3.multiply(adt_f5,temp_deriv_top2);
						adt_y_deriv.add(adt_temp2,adt_temp3);

						//						surrounding[2].add(adt_temp1,adt_temp2);
						adt_temp1.multiply(surrounding[3],temp_right);
						adt_right.multiply(adt_temp1,adt_f1);
						//						surrounding[3].multiply(temp_right);
						//						surrounding[3].multiply(adt_f1);
						adt_temp1.multiply(adt_f2,temp_rhs1);
						adt_rhs.add(adt_temp1,adt_f3);
						adt_temp1.add(adt_bottom,adt_left);
						adt_temp2.add(adt_center,adt_rhs);
						adt_temp2.add(adt_right);
						adt_temp3.add(adt_x_deriv,adt_y_deriv);

						adt_fun.add(adt_temp1,adt_temp2);
						adt_fun.add(adt_temp3);

						double val = adt_fun.get_first_derivative();
						if(l == 0){
							spm.update_matrix(val,row,col_bottom);
						}
						else if(l == 1){
							spm.update_matrix(val,row,col_left);
						}
						else if(l == 2){
							spm.update_matrix(val,row,col_center);

						}
						else if(l == 3){
							spm.update_matrix(val,row,col_right);
						}
						surrounding[l].ADConst();
					}
					B[row] = adt_fun.get_fun_value();
					delete surrounding;
				}
			}
			else if (j == number_of_divisions_x - 1){							// right face
				if((fabs(beeta[1] - 0.0) < pow(10.0,-6))){
					temp_center = alpha[1];
					spm.update_matrix(temp_center,row, col_center);
					B[row] = alpha[1]*val[row] - gama[1];
				}
				else{
					k_before = get_k_before(i);
					k_after = get_k_after(i);
					h_before = get_h_before(j);
					double temp_center1 =  -(((2.0/(h_before*h_before))*(1 + (alpha[1]*h_before/beeta[1]))));
					double temp_center2 = -(((2.0/(k_after*k_before))));
					//					temp_center = -(((2.0/(h_before*h_before))*(1 + (alpha[1]*h_before/beeta[1]))) + ((2.0/(k_after*k_before))));
					temp_left = 2.0/(h_before*h_before);
					temp_bottom = 2.0/(k_before*(k_before + k_after));
					temp_top = 2.0/(k_after*(k_before + k_after));

					double temp_deriv_right1 = -2*h_before*alpha[1]/beeta[1];
					double temp_deriv_right2 = 2*h_before*gama[1]/beeta[1];

					double temp_rhs1 = 2.0*gama[1]/(beeta[1]*h_before);
//					ADT<double> surrounding[4];
					surrounding = new ADT<double>[4];
					surrounding[0].ADconst(val[col_bottom]);
					surrounding[1].ADconst(val[col_left]);
					surrounding[2].ADconst(val[col_center]);
					surrounding[3].ADconst(val[col_top]);

					for(int l = 0; l < 4; l++){
						make_constant_f1(val[col_center]);
						make_constant_f2(val[col_center]);
						make_constant_f3(val[col_center]);
						make_constant_f4(val[col_center]);
						make_constant_f5(val[col_center]);
						if(l == 2){
							make_variable_f1(val[col_center]);
							make_variable_f2(val[col_center]);
							make_variable_f3(val[col_center]);
							make_variable_f4(val[col_center]);
							make_variable_f5(val[col_center]);
						}
						surrounding[l].ADvar();
						adt_temp1.multiply(surrounding[0],temp_bottom);
						adt_bottom.multiply(adt_temp1,adt_f2);
						//						surrounding[0].multiply(temp_bottom);
						//						surrounding[0].multiply(adt_f2);
						adt_temp1.multiply(surrounding[1],temp_left);
						adt_left.multiply(adt_temp1,adt_f1);
						//						surrounding[1].multiply(temp_left);
						//						surrounding[1].multiply(adt_f1);
						adt_temp1.multiply(surrounding[2],temp_center1);
						adt_temp2.multiply(surrounding[2],temp_center2);
						adt_temp1.multiply(adt_f1);
						adt_temp2.multiply(adt_f2);
						adt_center.add(adt_temp1,adt_temp2);
						//						surrounding[2].add(adt_temp1,adt_temp2);
						adt_temp1.multiply(surrounding[3],temp_top);
						adt_top.multiply(adt_temp1,adt_f2);
						//						surrounding[3].multiply(temp_top);
						//						surrounding[3].multiply(adt_f2);

						adt_temp1.multiply(surrounding[2],temp_deriv_right1);
						adt_temp2.multiply(adt_temp1,adt_f4);
						adt_temp3.multiply(adt_f4,temp_deriv_right2);
						adt_x_deriv.add(adt_temp2,adt_temp3);

						adt_temp1.multiply(surrounding[3],adt_f5);
						adt_temp2.multiply(surrounding[0],adt_f5);
						adt_y_deriv.substract(adt_temp1,adt_temp2);

						adt_temp1.multiply(adt_f1,temp_rhs1);
						adt_rhs.add(adt_temp1,adt_f3);
						adt_temp1.add(adt_bottom,adt_left);
						adt_temp2.add(adt_center,adt_top);
						adt_temp2.add(adt_rhs);

						adt_temp3.add(adt_x_deriv,adt_y_deriv);
						adt_fun.add(adt_temp1,adt_temp2);
						adt_fun.add(adt_temp3);

						double val = adt_fun.get_first_derivative();
						if(l == 0){
							spm.update_matrix(val,row,col_bottom);
						}
						else if(l == 1){
							spm.update_matrix(val,row,col_left);
						}
						else if(l == 2){
							spm.update_matrix(val,row,col_center);



						}
						else if(l == 3){
							spm.update_matrix(val,row,col_top);
						}
						surrounding[l].ADConst();
					}
					B[row] = adt_fun.get_fun_value();
					delete surrounding;
				}
			}
			else if (j == 0){													// left face
				if((fabs(beeta[3] - 0.0) < pow(10.0,-6))){
					temp_center = alpha[3];
					spm.update_matrix(temp_center,row, col_center);
					B[row] = alpha[3]*val[row] - gama[3];

				}
				else{
					h_after = get_h_after(j);
					k_before = get_k_before(i);
					k_after = get_k_after(i);
					double temp_center1 = -(((2.0/(h_after*h_after))*(1 - (alpha[3]*h_after/beeta[3]))));
					double temp_center2 = -(((2.0/(k_after*k_before))));
					//					temp_center = -(((2.0/(h_after*h_after))*(1 - (alpha[3]*h_after/beeta[3]))) + ((2.0/(k_after*k_before))));
					temp_right = 2.0/(h_after*h_after);
					temp_top = 2.0/(k_after*(k_before + k_after));
					temp_bottom = 2.0/(k_before*(k_after + k_before));

					double temp_deriv_right1 = -2*h_before*alpha[3]/beeta[3];
					double temp_deriv_right2 = 2*h_before*gama[3]/beeta[3];

					double temp_rhs1 = -(2.0*gama[3]/(beeta[3]*h_after));
//					ADT<double> surrounding[4];
					surrounding = new ADT<double>[4];
					surrounding[0].ADconst(val[col_bottom]);
					surrounding[1].ADconst(val[col_center]);
					surrounding[2].ADconst(val[col_right]);
					surrounding[3].ADconst(val[col_top]);

					for(int l = 0; l < 4; l++){
						make_constant_f1(val[col_center]);
						make_constant_f2(val[col_center]);
						make_constant_f3(val[col_center]);
						make_constant_f4(val[col_center]);
						make_constant_f5(val[col_center]);

						if(l == 1){
							make_variable_f1(val[col_center]);
							make_variable_f2(val[col_center]);
							make_variable_f3(val[col_center]);
							make_variable_f4(val[col_center]);
							make_variable_f5(val[col_center]);
						}
						surrounding[l].ADvar();
						adt_temp1.multiply(surrounding[0],temp_bottom);
						adt_bottom.multiply(adt_temp1,adt_f2);
						//						surrounding[0].multiply(temp_bottom);
						//						surrounding[0].multiply(adt_f2);
						adt_temp1.multiply(surrounding[1],temp_center1);
						adt_temp2.multiply(surrounding[1],temp_center2);
						adt_temp1.multiply(adt_f1);
						adt_temp2.multiply(adt_f2);
						adt_center.add(adt_temp1,adt_temp2);
						//						surrounding[1].add(adt_temp1,adt_temp2);
						adt_temp1.multiply(surrounding[2],temp_right);
						adt_right.multiply(adt_temp1,adt_f1);
						//						surrounding[2].multiply(temp_right);
						//						surrounding[2].multiply(adt_f1);
						adt_temp1.multiply(surrounding[3],temp_top);
						adt_top.multiply(adt_temp1,adt_f2);

						adt_temp1.multiply(surrounding[1],temp_deriv_right1);
						adt_temp2.multiply(adt_temp1,adt_f4);
						adt_temp3.multiply(adt_f4,temp_deriv_right2);
						adt_x_deriv.add(adt_temp2,adt_temp3);

						adt_temp1.multiply(adt_f5,surrounding[3]);
						adt_temp2.multiply(adt_f5,surrounding[0]);
						adt_y_deriv.substract(adt_temp1,adt_temp2);

						//						surrounding[3].multiply(temp_top);
						//						surrounding[3].multiply(adt_f2);
						adt_temp1.multiply(adt_f1,temp_rhs1);
						adt_rhs.add(adt_temp1,adt_f3);
						adt_temp1.add(adt_center,adt_right);
						adt_temp2.add(adt_top,adt_bottom);
						adt_temp2.add(adt_rhs);
						adt_temp3.add(adt_x_deriv,adt_y_deriv);

						adt_fun.add(adt_temp1,adt_temp2);
						adt_fun.add(adt_temp3);

						double val = adt_fun.get_first_derivative();
						if(l == 0){
							spm.update_matrix(val,row,col_bottom);
						}
						else if(l == 1){
							spm.update_matrix(val,row,col_center);

						}
						else if(l == 2){
							spm.update_matrix(val,row,col_right);
						}
						else if(l == 3){
							spm.update_matrix(val,row,col_top);
						}
						surrounding[l].ADConst();
					}
					B[row] = adt_fun.get_fun_value();
					delete surrounding;
				}

			}
			else{																					// center face
				h_after = get_h_after(j);
				h_before = get_h_before(j);
				k_after = get_k_after(i);
				k_before = get_k_before(i);


				double temp_center1 = -2.0*((1.0/(h_after*h_before)));
				double temp_center2 = -2.0* ((1.0/(k_after*k_before)));

				temp_bottom = 2.0/(k_before*(k_before + k_after));
				temp_right = 2.0/(h_after*(h_before + h_after));
				temp_top = 2.0/(k_after*(k_after + k_before));
				temp_left = 2.0/(h_before*(h_after + h_before));


				surrounding = new ADT<double>[5];
				surrounding[0].ADconst(val[col_bottom]);
				surrounding[1].ADconst(val[col_left]);
				surrounding[2].ADconst(val[col_center]);
				surrounding[3].ADconst(val[col_right]);
				surrounding[4].ADconst(val[col_top]);

				for(int l = 0 ;l < 5; l++){
					make_constant_f1(val[col_center]);
					make_constant_f2(val[col_center]);
					make_constant_f3(val[col_center]);
					if(l == 2){
						make_variable_f1(val[col_center]);
						make_variable_f2(val[col_center]);
						make_variable_f3(val[col_center]);
					}
					surrounding[l].ADvar();
					adt_temp1.multiply(surrounding[0],temp_bottom);
					adt_bottom.multiply(adt_temp1,adt_f2);
					adt_temp1.multiply(surrounding[1],temp_left);
					adt_left.multiply(adt_temp1,adt_f1);
					adt_temp1.multiply(surrounding[2],temp_center1);
					adt_temp2.multiply(surrounding[2],temp_center2);
					adt_temp1.multiply(adt_f1);
					adt_temp2.multiply(adt_f2);
					adt_center.add(adt_temp1,adt_temp2);
//					if(l == 2){
//						cout << adt_temp1.get_first_derivative() << endl;
//						cout << adt_temp2.get_first_derivative() << endl;
//						cout << adt_center.get_first_derivative() << endl;
//					}

					adt_temp1.multiply(adt_f4,surrounding[3]);
					adt_temp2.multiply(adt_f4,surrounding[1]);
					adt_x_deriv.substract(adt_temp1,adt_temp2);

					adt_temp1.multiply(adt_f5,surrounding[4]);
					adt_temp2.multiply(surrounding[0],adt_f5);
					adt_y_deriv.substract(adt_temp1,adt_temp2);

					adt_temp1.multiply(surrounding[3],temp_right);
					adt_right.multiply(adt_temp1,adt_f1);
					adt_temp1.multiply(surrounding[4],temp_top);
					adt_top.multiply(adt_temp1,adt_f2);
					adt_temp1.add(adt_left,adt_right);
					adt_temp2.add(adt_top,adt_bottom);
					adt_temp3.add(adt_temp1,adt_temp2);
//					if(l == 2){
//						cout << adt_temp1.get_first_derivative() << endl;
//						cout << adt_temp2.get_first_derivative() << endl;
//						cout << adt_temp3.get_first_derivative() << endl;
//					}


					adt_temp3.add(adt_f3);
					adt_temp4.add(adt_x_deriv,adt_y_deriv);


					adt_fun.add(adt_temp3,adt_center);
					adt_fun.add(adt_temp4);

					double val = adt_fun.get_first_derivative();
					if(l == 0){
						spm.update_matrix(val,row,col_bottom);
					}
					else if(l == 1){
						spm.update_matrix(val,row,col_left);
					}
					else if(l == 2){

						spm.update_matrix(val,row,col_center);

					}
					else if(l == 3){
						spm.update_matrix(val,row,col_right);
					}
					else{
						spm.update_matrix(val,row,col_top);
					}
					surrounding[l].ADConst();

				}
				B[row] = adt_fun.get_fun_value();
				delete surrounding;
			}

		}

	}

	spm.resize_matrix();
//	spm.display();
	//	spm.display();

//	B.display();
}


//! Displays the sparse matrix
void ComputeMatrix:: compute_matrix_display(){
	spm.display();
}

//! Solve the equation AX = B
void ComputeMatrix::solve_iterative(){
	Iterative_Solver its;
	cout << "Enter the solver you want to use" << endl;
	cout << "1. Gauss Siedel with or without relaxation" << endl;
	cout << "2. BiConjugate Gradient without preconditioner" << endl;
	cout << "3. Biconjugate Gradient with preconditioner" << endl;
	int choice;
	cin >> choice;
	if(choice == 1){
		clock_t t;
		t = clock();
		its.Gauss_Siedel(spm,X,B);
		t = clock() - t;
		cout << "Time taken by Gauss Seidel is  " << ((float)t)/CLOCKS_PER_SEC << "seconds" << endl;
	}
	else if (choice == 2){
		clock_t t;
		t = clock();
		its.BiCG(spm,X,B);
		t = clock() - t;
		cout << "Time taken by BiCG without preconditioner is  " << ((float)t)/CLOCKS_PER_SEC << "seconds" << endl;


	}
	else if (choice == 3){
		int choice2;
		cout << "0. Jacobi Preconditioner" << endl;
		cout << "1. ILU Preconditioner" << endl;
		cout << "2. Gauss Seidel Preconditioner" <<endl;
		cin >> choice2;
		clock_t t;
		t = clock();
		its.BiCG_preconditioner(spm,X,B,choice2);
		t = clock() - t;
		cout << "Time taken by BiCG with the choice of your preconditioner is  " << ((float)t)/CLOCKS_PER_SEC << "seconds" << endl;
	}

	//	cout <<"Reached Iter";
	//	X.display();


}


void ComputeMatrix::solveNonLinear(){
	compute_h();
	compute_k();
	number_of_divisions_x = x_spacing[number_of_x - 1][3];
	number_of_divisions_y = y_spacing[number_of_y - 1][3];
	cout << "*************" <<endl;
	cout << "Enter the boundary conditions" << endl;
	cout << "Boundary condition is of the form: alpha*u + beeta * u' = gama" << endl;
	cout << "Enter the boundary condition in order of bottom, right, top, left" << endl;
	get_boundary_condition();
	clock_t t;
	t = clock();
	X.set_size((number_of_divisions_x )*(number_of_divisions_y ));
	X.ones();
	Vector<double> y((number_of_divisions_x )*(number_of_divisions_y ));
	B.set_size((number_of_divisions_x )*(number_of_divisions_y));

	Vector<double> X1((number_of_divisions_x )*(number_of_divisions_y ));
//	spm.display();
//	cout << "*****";
//	B.display();
		ofstream fout("size_mat.txt");

		fout << number_of_divisions_x << endl;
		fout << number_of_divisions_y << endl;
		fout.close();
		cout << "File Written for size" << endl;
//		cout <<"*******" << endl;
		Iterative_Solver its;
		double err = 199;
		int number_of_iter = 0;
		while (err > pow(10.0,-5)){
////			cout <<  "Reached";
//
			spm.del();
			spm.initialize((number_of_divisions_x)*(number_of_divisions_y));
			compute(X);
//			X.display();
			cout <<endl;
//			spm.display();
//			cout << endl;
			y.copy_values(X);
//
			B.minus();
			its.BiCG_preconditioner(spm,y,B,1);
////			y.display();
//			cout << endl;
////			B.display();
//			cout << endl;
////			spm.display();
			X1.add(X,y);
			err = cal_error(B);
			X.copy_values(X1);
//
			number_of_iter++;
			cout << "Iteration of Newton box:" << number_of_iter << endl;
			cout << err << endl;

		}
//		//	spm.display();
//		//	X.display();
//
//
//
		t = clock() - t;
		cout << "Writing in file" << endl;
		ofstream fout1("solver.txt");

		for(int i = 0; i < X.get_size(); i++){
			fout1 << X[i] << endl;
		}
		fout1.close();
		cout << "File Written" << endl;

}
//! This calculates the error based on the relative difference between two points
/*! @param X Vector X */
/*! @param X1 Vector X1*/

double ComputeMatrix::cal_error(Vector<double> &X, Vector<double> &X1){
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

double ComputeMatrix::cal_error(Vector<double> &X){
	double err = 0.0;

	for (int i = 0; i < X.get_size(); i++){
		if (fabs(X[i]) > err){
			err = fabs(X[i]);
		}

	}
	//	cout << err << endl;
	return err;
}
