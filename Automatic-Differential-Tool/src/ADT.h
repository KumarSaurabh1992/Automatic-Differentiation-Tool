/*
 * ADT.h
 *
 *  Created on: Feb 27, 2015
 *      Author: saurabh
 */

#ifndef ADT_H_
#define ADT_H_

#include <iostream>
#include <math.h>
using namespace std;

template <class T> class ADT{
	T fun_val, diff1_val, diff2_val;
public:
	ADT<T>();
	T get_first_derivative();
	T get_fun_value();
	T get_second_derivative();
	//	void operator = (ADT &);
	void multiply(T);
	void multiply(ADT &, ADT &);
	void multiply(ADT &);
	void multiply(ADT &, T );
	void ADvar(T);
	void ADvar();
	void ADconst(T);
	void ADconst(ADT &);
	void ADConst();
	void set_fun_value(T);
	void set_first_derivative(T);
	void set_second_derivative(T);
	void power(ADT &, int);
	void sine(ADT &);
	void add(ADT &, ADT &);
	void add(ADT &);
	void add(ADT &,T);
	void substract(ADT &, ADT &);
	void substract(ADT &,T);
	void cosin(ADT &);
	void exponential(ADT &);
	void log(ADT &);
	void division(ADT &, ADT &);
	//	ADT operator + (ADT &);


};

template <class T> ADT<T> :: ADT(){
	fun_val = 0;
	diff1_val = 0;
	diff2_val = 0;
}

template <class T> void ADT<T> :: multiply(T value){
	fun_val = value * fun_val;
	diff1_val = value * diff1_val;
	diff2_val = value*diff2_val;

}

template <class T> void ADT<T> :: multiply(ADT<T> & adt1,T value){
	fun_val = value * adt1.fun_val;
	diff1_val = value * adt1.diff1_val;
	diff2_val = value * adt1.diff2_val;

}

template <class T> void ADT<T> :: multiply(ADT <T> & adt1, ADT<T> & adt2){
	fun_val = adt1.get_fun_value()*adt2.get_fun_value();
	diff1_val = adt1.get_fun_value()*adt2.get_first_derivative() + adt1.get_first_derivative()*adt2.get_fun_value();
	diff2_val = adt1.get_second_derivative()*adt2.get_fun_value() + 2.0*adt1.get_first_derivative()*adt2.get_first_derivative() + adt1.get_fun_value()*adt2.get_second_derivative();
}

template <class T> void ADT<T> :: multiply(ADT <T> & adt1){

	diff1_val = fun_val*adt1.get_first_derivative() + diff1_val*adt1.get_fun_value();
	diff2_val = diff2_val*adt1.get_fun_value() + 2.0*diff1_val*adt1.diff1_val + fun_val*adt1.get_second_derivative();
	fun_val = fun_val*adt1.get_fun_value();
}
template <class T> T ADT<T> :: get_first_derivative(){
	return diff1_val;
}

template <class T> void ADT<T> :: ADvar(T val){
	fun_val = val;
	diff1_val = 1.0;
	diff2_val = 0.0;
}

template <class T> void ADT<T> :: ADvar(){
	diff1_val = 1.0;
}

template <class T> void ADT<T> :: ADConst(){
	diff1_val = 0.0;
	diff2_val = 0.0;
}
template<class T> T ADT<T> :: get_fun_value(){
	return fun_val;
}

template <class T> T ADT<T> :: get_second_derivative(){
	return diff2_val;
}


template<class T> void ADT<T>::ADconst(T val){
	fun_val = val;
	diff1_val = 0.0;
	diff2_val = 0.0;
}

template<class T> void ADT<T>::ADconst(ADT<T> & adt1){
	fun_val = adt1.fun_val;
	diff1_val = 0.0;
	diff2_val = 0.0;
}
template<class T> void ADT<T>:: power(ADT<T> & adt, int val){
	fun_val = pow(adt.get_fun_value(),val*1.0);
	diff1_val = val*pow(adt.get_fun_value(),(val - 1)*1.0)*adt.get_first_derivative();
	diff2_val = val*(val - 1)*pow(adt.get_fun_value(),(val - 2)*1.0)*adt.get_first_derivative()*adt.get_first_derivative() + val*pow(adt.get_fun_value(),(val - 1)*1.0)*adt.get_second_derivative();
}

template<class T> void ADT<T>::sine(ADT<T> & adt1){
	fun_val = sin(adt1.get_fun_value());
	diff1_val = cos(adt1.get_fun_value())*adt1.get_first_derivative();
	diff2_val = -sin(adt1.get_fun_value())*adt1.get_first_derivative()*adt1.get_first_derivative() + cos(adt1.get_fun_value())*adt1.get_second_derivative();
}

template<class T> void ADT<T>::add(ADT<T> & adt1, ADT<T> & adt2){
	fun_val = adt1.get_fun_value() + adt2.get_fun_value();
	diff1_val = adt1.get_first_derivative() + adt2.get_first_derivative();
	diff2_val = adt1.get_second_derivative() + adt2.get_second_derivative();
}

template<class T> void ADT<T>::add(ADT<T> & adt1){
	fun_val = adt1.get_fun_value() + fun_val;
	diff1_val = adt1.get_first_derivative() + diff1_val;
	diff2_val = adt1.get_second_derivative() + diff2_val;
}

template<class T> void ADT<T>::add(ADT<T> & adt,T val){
	fun_val = adt.fun_val + val;
	diff1_val = adt.diff1_val;
	diff2_val = adt.diff2_val;
}

template<class T> void ADT<T>::substract(ADT<T> & adt1, ADT & adt2){
	fun_val = adt1.get_fun_value() - adt2.get_fun_value();
	diff1_val = adt1.get_first_derivative() - adt2.get_first_derivative();
	diff2_val = adt1.get_second_derivative() - adt2.get_second_derivative();

}

template<class T> void ADT<T>::substract(ADT<T> & adt, T val){
	fun_val = adt.fun_val - val;
	diff1_val = adt.diff1_val;
	diff2_val = adt.diff2_val;

}

//template<class T> ADT<T> ADT<T>::operator +(ADT<T> & adt1){
//	ADT <T> nadt;
//	nadt.fun_val = fun_val + adt1.fun_val;
//	nadt.diff1_val = diff1_val + adt1.diff1_val;
//	nadt.diff2_val = diff2_val + adt1.diff2_val;
//	return nadt;
//}

template <class T> void ADT<T>::cosin(ADT<T> & adt1){
	fun_val = cos(adt1.get_fun_value());
	diff1_val = -sin(adt1.get_fun_value())*adt1.get_first_derivative();
	diff2_val = -cos(adt1.get_fun_value())*adt1.get_first_derivative()*adt1.get_first_derivative() - cos(adt1.get_fun_value())*adt1.get_second_derivative();
}

template <class T> void ADT<T>::exponential(ADT<T> & adt1){
	fun_val = exp(adt1.get_fun_value());
	diff1_val = exp(adt1.get_fun_value())*adt1.get_first_derivative();
	diff2_val = exp(adt1.get_fun_value())*adt1.get_first_derivative()*adt1.get_first_derivative() + exp(adt1.get_fun_value())*adt1.get_second_derivative();
}

template <class T> void ADT<T>::log(ADT & adt1){
	fun_val = adt1.fun_val;
	diff1_val = adt1.diff1_val/adt1.fun_val;
	diff2_val = (adt1.fun_val*adt1.diff2_val - adt1.diff1_val*adt1.diff1_val)/(adt1.fun_val*adt1.fun_val);
}

template <class T> void ADT<T>::division(ADT<T> & adt1, ADT<T> & adt2){
	fun_val = adt1.fun_val/adt2.fun_val;
	diff1_val = (adt1.diff1_val*adt2.fun_val - adt1.fun_val*adt2.diff1_val) /(adt2.fun_val*adt2.fun_val);
	diff2_val = ((adt2.fun_val*(adt2.fun_val*adt1.diff2_val - adt1.fun_val*adt2.diff2_val)) - 2*adt2.diff1_val*(adt1.diff1_val*adt2.fun_val - adt1.fun_val*adt2.diff1_val))/pow(adt2.fun_val,3.0);
}


#endif /* ADT_H_ */
