#pragma once
#include "Vector_new.h"

void print(std::vector<double> x)
{
	for (int i = 0; i != x.size(); i++)
		std::cout << x[i] << '\t';
	std::cout << std::endl;
}

std::vector<double> operator * (double a, std::vector<double> b)
{
	std::vector<double> c(b);
	for (int i = 0; i != c.size(); i++)
		c[i] = a * b[i];
	return c;
}

std::vector<double> operator + (std::vector<double> a, std::vector<double> b)
{
	std::vector<double> c(b);
	for (int i = 0; i != c.size(); i++)
		c[i] = a[i] + b[i];
	return c;
}

std::vector<double> operator - (std::vector<double> a, std::vector<double> b)
{
	std::vector<double> c(b);
	for (int i = 0; i != c.size(); i++)
		c[i] = a[i] - b[i];
	return c;
}

std::vector<double> operator / (std::vector<double> a, double b)
{
	std::vector<double> c(a);
	for (int i = 0; i != c.size(); i++)
		c[i] = a[i] / b;
	return c;
}