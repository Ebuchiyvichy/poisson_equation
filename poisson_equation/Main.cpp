#include "vector_new.h"

double	EPS = 1.e-3;
double	PI = 3.1415926535;

int main()
{
	Data	DATA;
	double	h1 = 0.01;
	double	h2 = 0.01;
	double	tau = 0.01;
	double	start = 0.0;

	DATA.initial = &init;
	DATA.left_boarder = &phi;
	DATA.right_boarder = &psi;
	DATA.front_boader = &f;
	DATA.back_boader = &g;
	DATA.f_function = &r;

	poisson_equation(h1, h2, tau, DATA, "", 0.0);

	//error_check(0.1, 0.01, DATA, start);


	system("pause");
	return (0);
}