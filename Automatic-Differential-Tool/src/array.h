/*
 * array.h
 *
 *  Created on: Jan 16, 2015
 *      Author: saurabh
 */

#ifndef ARRAY_H_
#define ARRAY_H_

#include <iostream>
using namespace std;

template <class T> class array {
	int row,column;
	T **A;
	void getData();
	void initialize();
public:
	array<T>(int m,int n);
	array<T>();
	~array<T>();
	void display();
	void read_values();
	array operator + (array);
	array operator - (array);
	array operator * (array);
	array add(array);
	void set_size(int ,int);
	void set_element(T,int ,int);
};


template <class T> array<T>::array(){
	row = 0;
	column = 0;
	A = new T*[0];
	for (int i = 0; i < row; i++)
	{
		A[i] = new T[column];
	}
}
template <class T> array<T>::array(int m, int n){
	row = m;
	column = n;
	A = new T*[row];
	for (int i = 0; i < row; i++)
	{
		A[i] = new T[column];
	}
}

template <class T> void array<T>::set_element(T val, int ROW, int COL){
	A[ROW][COL] = val;
}
template<class T> void array<T>::set_size(int m, int n){
	row = m;
	column = n;
	A = new T*[row];
	for (int i = 0; i < row; i++)
	{
		A[i] = new T[column];
	}
}

template <class T> array<T>:: ~array<T>(){
}
template <class T> void array<T>::read_values(){
	getData();
}
template <class T> void array<T>::display()
		{
	for (int i = 0; i < row; i++)
	{
		for(int j = 0; j < column; j++)
		{
			cout<< A[i][j] << "\t";
		}
		cout << endl;
	}

		}

template <class T> void array<T>::getData()
		{
	cout << "Enter the value for the matrix" << endl;
	for (int i = 0; i < row; i++)
	{
		for(int j = 0; j < column; j++)
		{
			cin>> A[i][j];
		}
	}
		}

template <class T> array<T> array<T>::operator + (array matrix2){
	array<T> temp_matrix(row,column);
	for (int i = 0; i < row; i++)
	{
		for(int j = 0; j < column; j++)
		{
			temp_matrix.A[i][j] = A[i][j] + matrix2.A[i][j];

		}

	}
	return (temp_matrix);
}

template <class T> array<T> array<T>::operator - (array matrix2){
	array<T> temp_matrix(row,column);
	for (int i = 0; i < row; i++)
	{
		for(int j = 0; j < column; j++)
		{
			temp_matrix.A[i][j] = A[i][j] - matrix2.A[i][j];

		}
	}
	return (temp_matrix);
}

template <class T> array<T> array<T>::operator * (array matrix2){
	array<T> temp_matrix(row,matrix2.column);
	temp_matrix.initialize();
	for (int i = 0; i < row; i++)
	{
		for(int j = 0; j < matrix2.column; j++)
		{
			for(int k = 0; k < column; k++)
			{
				temp_matrix.A[i][k] += A[i][j]*matrix2.A[j][k];
			}
		}
	}
	return (temp_matrix);
}
//template <class T> void array<T>::add(array matrix2)
//{
//	matrix1
//}

template <class T> void array<T> :: initialize(){
	for (int i = 0; i < row; i++){
		for(int j = 0; j < column; j++){
			A[i][j] = 0;
		}
	}
}


#endif /* ARRAY_H_ */
