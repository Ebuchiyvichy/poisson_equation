#include "vector_new.h"

double	EPS = 1.e-3;
double	PI = 3.1415926535;

int main()
{
	Data	DATA;
	double	h1 = 0.1;
	double	h2 = 0.1;
	double	tau = 0.1;
	double	start = 0.1;
	std::ofstream		fout_err;

	DATA.initial = &init;
	DATA.left_boarder = &phi;
	DATA.right_boarder = &psi;
	DATA.front_boader = &f;
	DATA.back_boader = &g;
	DATA.f_function = &r;

	fout_err.open("Error on diff step.txt");


//	poisson_equation(h1, h2, tau, DATA, "", 0.0);

	for (int i = 0; i <= 4; i++)
		fout_err << i << '\t' << poisson_equation(h1 / pow(2, i), h2 / pow(2, i), tau / pow(2, i), DATA, "", 0.0) << '\n';


	system("pause");
	return (0);
}