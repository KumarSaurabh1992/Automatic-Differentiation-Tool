//! This class computes the matrix for the poison process.
/*! @author Kumar Saurabh (MA14M004) */

#ifndef COMPUTEMATRIX_H_
#define COMPUTEMATRIX_H_

#include "SparseMatrix.h"
#include "IterativeSolver.h"
#include "ADT.h"

class ComputeMatrix {
	//! A in sparse Matrix Form
	SparseMatrix spm;
	//! Number of x divisions
	int number_of_divisions_x;
	//! Number of y divisions
	int number_of_divisions_y;
	//! x spacing
	double h;
	//! y spacing
	double k;
	//! Vector B in AX = B
	Vector<double> B;
	//! Vector X in AX = B
	Vector<double> X;
	//! Value of alpha in the boundary condition of all the faces
	Vector<double> alpha;
	//! Value of beta in the boundary condition of all the faces
	Vector<double> beeta;
	//! Value of gama in the boundary condition of all the faces
	Vector<double> gama;
	//! Vector containing different x spacing
	double **x_spacing;
	//! Vector containing various y spacing
	double **y_spacing;
	//! number of different spacing in x
	int number_of_x;
	//! number of different spacing in y
	int number_of_y;
	ADT<double> * surrounding;
	ADT<double> adt_f1, adt_f2, adt_f3, adt_f4, adt_f5;
	void set_boundary_size(int N);
	void compute_h();
	void compute_k();
	double get_h_before(int );
	double get_h_after(int );
	double get_k_before(int );
	double get_k_after(int );
	void get_boundary_condition();
	double cal_error(Vector<double> &, Vector<double> &);
	double cal_error(Vector<double> &);
	void make_variable_f1(double val);
	void make_constant_f1(double val);
	void make_variable_f2(double val);
	void make_constant_f2(double val);
	void make_variable_f3(double val);
	void make_constant_f3(double val);
	void make_variable_f4(double val);
	void make_constant_f4(double val);
	void make_variable_f5(double val);
	void make_constant_f5(double val);

public:
	void compute(Vector<double> &);
	void set_number_of_division_x(int N);
	void set_number_of_division_y(int N);
	int get_number_of_division_x();
	int get_number_of_division_y();
	float get_h();
	float get_k();
	void compute_matrix_display();
	void solve_iterative();
	void compute_general(Vector <double> &);
	void solveNonLinear();



};

#endif /* COMPUTEMATRIX_H_ */
