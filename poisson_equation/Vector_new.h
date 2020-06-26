#pragma once
#ifndef VECTOR_NEW_H
#define VECTOR_NEW_H

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include <string>

extern double	EPS;
extern double	PI;

struct Data
{
	double	L1 = 1;//test1,test2,test3
	double	L2 = 1;

	double	T = 0.1;

	bool	condition[4] = { 1, 1, 0, 0 };//левое/правое/нижнее/верхнее условия; 0 -- первый род, 1 -- второй

	std::function<double(double, double, Data)>	initial;
	std::function<double(double, Data)>	left_boarder;
	std::function<double(double, Data)>	right_boarder;
	std::function<double(double, Data)>	front_boader;
	std::function<double(double, Data)>	back_boader;
	std::function<double(double, double, Data)>	f_function;
};

void print(std::vector<double> x);

// переодпределение операций под вектора
std::vector<double> operator * (double a, std::vector<double> b);
std::vector<double> operator + (std::vector<double> a, std::vector<double> b);
std::vector<double> operator - (std::vector<double> a, std::vector<double> b);
std::vector<double> operator / (std::vector<double> a, double b);


double	init(double x, double y, Data my_data);
double	f(double x, Data my_data);
double	g(double x, Data my_data);
double	phi(double t, Data my_data);
double	psi(double t, Data my_data);
double	r(double x, double y, Data my_data);

void	poisson_equation(double h1, double h2, double tau, Data my_data, std::string order, double start);



#endif //VECTOR_NEW_H