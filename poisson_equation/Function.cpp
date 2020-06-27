#pragma once
#include "Vector_new.h"

double	init(double x, double y, Data my_data)
{
//	return 1;
	return 1 + x;//test1
}
double	f(double x, Data my_data)
{
//	return 1;//test1
	return -1;//test2
//	return x * x;//test3
//	return -sin(x) * sin(x);//var2
	return -sqrt(1 + x);//var6
	return 0.0;//var10
}
double	g(double x, Data my_data)
{
//	return 1;//test1
	return 1;//test2
//	return 1 + x * x;//test3
//	return sin(x) * sin(x);//var2
	return sqrt(2 + x);//var6
	return x * (x - 1);//var10
}
double	phi(double x, Data my_data)
{
//	return 1;//test1
	return 1 + x;//test2
//	return 0.0;//test3
//	return 0.;//var2
	return 2.0 / 3 * pow(x + 1, 1.5);//var6
	return x * x;//var10
}
double	psi(double x, Data my_data)
{
//	return 1;//test1
	return 1 + x;//test2
//	return 2.;//test3
//	return 0.75 * x;//var2
	return 2.0 / 3 * pow(x + 2, 1.5);//var6
	return 4 * x - x * x;//var10
}
double	r(double x, double y, Data my_data)
{
	return 0.0;//test1,test2
//	return -4.0;//test3
//	return -y * cos(2 * x);//var2
	return -1 / sqrt(x + y + 1);//var6
	return -(2 * y - x);//var10
}

double	test(double x, double y)
{
//	return 1.0;//test1
	return 1 + y;//test2
//	return y * sin(x) * sin(x);//var2
	return 2. / 3 * pow(x + y + 1, 1.5);//var6
	return x * x*y - y * y*x;//var10
}

// правая прогонка (половинная прогонка)
void progon_05(int n1, int n2, int k, double **y, double *a, double *b, double *c, double *d, Data my_data)
{
	std::vector<double> alfa(n1 + 1);
	std::vector<double> betta(n1 + 1);

	alfa[1] = c[0] / b[0]; betta[1] = d[0] / b[0];
//	std::cout << c[0] << '\t' << d[0] << '\t';
	for (int i = 1; i < n1; i++) {
		alfa[i + 1] = c[1] / (b[1] - a[1] * alfa[i]);
		betta[i + 1] = (a[1] * betta[i] + d[i]) / (-a[1] * alfa[i] + b[1]);
//		std::cout << alfa[i + 1] << '\t' << betta[i + 1] << '\t';
	}
//	std::cout << n1;
//	std::cout << '\n' << '\n';

	y[n1][k] = (d[n1] + a[2] * betta[n1]) / (b[2] - a[2] * alfa[n1]);
	for (int i = n1 - 1; i >= 0; i--)
		y[i][k] = alfa[i + 1] * y[i + 1][k] + betta[i + 1];
	
//	for (int i = 0; i <= n1; i++)
//	{
//		std::cout << y[i][k] << '\n';
//	}
}

void progon_1(int n1, int n2, int k, double **y, double *a, double *b, double *c, double *d, Data my_data)
{
	std::vector<double> alfa(n2 + 1);
	std::vector<double> betta(n2 + 1);

	alfa[1] = c[0] / b[0]; betta[1] = d[0] / b[0];
//	std::cout << alfa[1] << '\t' << betta[1] << '\t';

	for (int i = 1; i < n2; i++) {
		alfa[i + 1] = c[1] / (b[1] - a[1] * alfa[i]);
		betta[i + 1] = (a[1] * betta[i] + d[i]) / (-a[1] * alfa[i] + b[1]);
	//	std::cout << alfa[i+1] << '\t' << betta[i+1] << '\t';

	}
//	std::cout << '\n' << '\n';

	y[k][n2] = (d[n2] + a[2] * betta[n2]) / (b[2] - a[2] * alfa[n2]);

	for (int i = n2 - 1; i >= 0; i--)
	{
		y[k][i] = alfa[i + 1] * y[k][i + 1] + betta[i + 1];
//		std::cout << y[k][i] << '\n';
	}
//	for (int i = 0; i <= n2; i++)
//	{
//		std::cout << y[k][i] << '\n';
//	}
//	std::cout << '\n';
}

double	max_el(double **y1, double **y2, int n1, int n2)
{
	double max = -100;

	for (int i = 0; i <= n1; i++)
		for (int j = 0; j <= n2; j++)
			if (fabs(y1[i][j] - y2[i][j]) > max)
				max = fabs(y1[i][j] - y2[i][j]);
	return max;
}

double	max_err(double **y1, int n1, int n2, double h1, double h2)
{
	double max = -100;

	for (int i = 0; i <= n1; i++)
		for (int j = 0; j <= n2; j++)
			if (fabs(y1[i][j] - test(i * h1, j * h2)) > max)
				max = fabs(y1[i][j] - test(i * h1, j * h2));
	return max;
}

//для ГУ только первого рода
double	poisson_equation(double h1, double h2, double tau, Data my_data, std::string order, double start)
{
	int					n1 = my_data.L1 / h1;
	int					n2 = my_data.L2 / h2;
	double				max;
	//коэффициенты прогонки
	double				A1[3] = { 0.0, 
									1.0 / (h1 * h1), 
									my_data.condition[1] * 2.0 / (h1 *h1) };
	double				B1[3] = { -(1.0 - my_data.condition[0]) + 2. * my_data.condition[0] * (1. / (h1 * h1) + 1. / tau),
									2. * (1.0 / (h1 * h1) + 1.0 / tau),
									 -(1.0 - my_data.condition[1]) + 2. * my_data.condition[1] * (1. / (h1 * h1) + 1. / tau) };
	double				C1[3] = { my_data.condition[0] * 2.0 / (h1 * h1), 
									A1[1], 
									0.0 };

	//коэффициенты прогонки
	double				A2[3] = { 0.0, 
									1.0 / (h2 * h2), 
									my_data.condition[3] * 2.0 / (h2 * h2) };
	double				B2[3] = { -(1.0 - my_data.condition[2]) + 2. * my_data.condition[2] * (1. / (h2 * h2) + 1. / tau),
									2. * (1.0 / (h2 * h2) + 1.0 / tau),
									 -(1.0 - my_data.condition[3]) + 2. * my_data.condition[3] * (1. / (h2 * h2) + 1. / tau) };
	double				C2[3] = { my_data.condition[2] * 2.0 / (h2 * h2),
									A2[1], 
									0.0 };

	double				*F1 = new double[n1 + 1];
	double				*F2 = new double[n2 + 1];
	//файловые переменные
	std::ofstream		fout;
	std::ofstream		fout_err;

	//массивы по временным шагам
	double				**y1 = new double *[n1 + 1];	//нулевой слой
						for (int i = 0; i <= n1; i++)
							y1[i] = new double[n2 + 1];

	double				**y2 = new double *[n1 + 1];	//половинный слой
						for (int i = 0; i <= n1; i++)
							y2[i] = new double[n2 + 1];

	double				**y3 = new double *[n1 + 1];	//конечный слой
						for (int i = 0; i <= n1; i++)
							y3[i] = new double[n2 + 1];

	std::string			str = "Pois_eq";

	str += order + ".txt";
	fout.open(str);
	fout_err.open("Err.txt");
	// инициализация начальными данными
	for (int i = 0; i <= n1; i++)			//инициализация нулевого слоя
		for (int j = 0; j <= n2; j++)
			y1[i][j] = my_data.initial(i * h1, j * h2 , my_data);
	int it = 0;
	for (double t = tau; t <= my_data.T; t += tau)
//	do
	{

		//прогонка вдоль (половинчатый слой)	/*гранусловия левая и правая граница*/
		if (my_data.condition[2] == 0)	//начальное гранусловие -- первого рода
			for (int i = 0; i <= n1; i++)
				y2[i][0] = my_data.front_boader(i * h1, my_data);
		else							//начальное условие -- второго рода
		{
			F1[0] = -(1.0 - my_data.condition[0]) * my_data.left_boarder(0.0, my_data) + 
				my_data.condition[0] * (2. / tau * y1[0][0] + 2. / (h2 * h2) * (y1[0][1] - y1[0][0]) + 2. / h2 * my_data.front_boader(0.0, my_data) + 
					my_data.f_function(0.0, 0.0, my_data) + 2. / h1 * my_data.left_boarder(0.0, my_data));
			for (int i = 1; i < n1; i++)
				F1[i] = 2. / tau * y1[i][0] + 2. / (h2 * h2) * (y1[i][1] - y1[i][0]) + 2. / h2 * my_data.front_boader(i * h1, my_data) + my_data.f_function(i * h1, 0.0, my_data);
			
			F1[n1] = -(1.0 - my_data.condition[1]) * my_data.right_boarder(0.0, my_data) + 
				my_data.condition[1] * (2. / tau * y1[n1][0] + 2. / (h2 * h2) * (y1[n1][1] - y1[n1][0]) + 2. / h2 * my_data.front_boader(my_data.L1, my_data) + 
					my_data.f_function(my_data.L1, 0.0, my_data) + 2. / h1 * my_data.right_boarder(0.0, my_data));
			

			progon_05(n1, n2, 0, y2, A1, B1, C1, F1, my_data);
		}

		if (my_data.condition[3] == 0)  //начальное гранусловие -- первого рода
			for (int i = 0; i <= n1; i++)
				y2[i][n2] = my_data.back_boader(i * h1, my_data);
		else							//начальное условие -- второго рода
		{
			F1[0] = -(1.0 - my_data.condition[0]) * my_data.left_boarder(my_data.L2, my_data) + 
				my_data.condition[0] * (2. / tau * y1[0][n2] + 2. / (h2 * h2) * (y1[0][n2 - 1] - y1[0][n2]) + 
					2. / h2 * my_data.back_boader(0.0, my_data) + my_data.f_function(0.0, my_data.L2, my_data) + 2. / h1 * my_data.left_boarder(my_data.L2, my_data));

			for (int i = 1; i < n1; i++)
				F1[i] = 2. / tau * y1[i][n2] + 2. / (h2 * h2) * (y1[i][n2 - 1] - y1[i][n2]) + 2. / h2 * my_data.back_boader(i * h1, my_data) + my_data.f_function(i * h1, my_data.L2, my_data);

			F1[n1] = -(1.0 - my_data.condition[1]) * my_data.right_boarder(my_data.L2, my_data) + 
				my_data.condition[1] * (2. / tau * y1[n1][n2] +	2. / (h2 * h2) * (y1[n1][n2 - 1] - y1[n1][n2]) + 
					2. / h2 * my_data.back_boader(my_data.L1, my_data) + my_data.f_function(my_data.L1, my_data.L2, my_data) + 2. / h1 * my_data.right_boarder(my_data.L2, my_data));
			


			progon_05(n1, n2, n2, y2, A1, B1, C1, F1, my_data);
		}

		for (int i = 1; i != n2; i++)	//ходим поперек
		{
			F1[0] = -(1. - my_data.condition[0]) * my_data.left_boarder(i * h2, my_data) + 
				my_data.condition[0] * (2. / tau * y1[0][i] + (y1[0][i + 1] - 2. * y1[0][i] + y1[0][i - 1]) / (h2 * h2)
				+ my_data.f_function(0.0, i * h2, my_data) + 2. / h1 * my_data.left_boarder(i * h2, my_data));

			for (int j = 1; j < n1; j++)	//(по предыдущему слою)
				F1[j] = 2. / tau * y1[j][i] + (y1[j][i + 1] - 2. * y1[j][i] + y1[j][i - 1]) / (h2 * h2) + my_data.f_function(j * h1, i * h2, my_data);

			
			F1[n1] = -(1. - my_data.condition[1]) * my_data.right_boarder(i * h2, my_data) + 
				my_data.condition[1] * (2. / tau * y1[n1][i] + (y1[n1][i + 1] - 2. * y1[n1][i] + y1[n1][i - 1]) / (h2 * h2)
				+ my_data.f_function(my_data.L1, i * h2, my_data) + 2. / h1 * my_data.right_boarder(i * h2, my_data));

			progon_05(n1, n2, i, y2, A1, B1, C1, F1, my_data);
		}

		//прогонка поперек --- искомый слой /*гранусловия -- нижняя и верхняя граница*/
		if (my_data.condition[0] == 0)	//начальное гранусловие -- первого рода
			for (int i = 0; i <= n2; i++)	
				y3[0][i] = my_data.left_boarder(i * h2, my_data);
		else							//начальное условие -- второго рода
		{
			F2[0] = -(1.0 - my_data.condition[2]) * my_data.front_boader(0.0, my_data) + 
				my_data.condition[2] * (2. / tau * y2[0][0] + 2. / (h1 * h1) * (y2[1][0] - y2[0][0]) + 2. / h2 * my_data.front_boader(0.0, my_data) + 
					my_data.f_function(0.0, 0.0, my_data) + 2. / h1 * my_data.left_boarder(0.0, my_data));
	
			for (int i = 1; i < n2; i++)
				F2[i] = 2. / tau * y2[0][i] + 2. / (h1 * h1) * (y2[1][i] - y2[0][i]) + 2. / h1 * my_data.left_boarder(i * h2, my_data) + my_data.f_function(0.0, i*h2, my_data);

			F2[n2] = -(1.0 - my_data.condition[3]) * my_data.back_boader(0., my_data) + 
				my_data.condition[3] * (2. / tau * y2[0][n2] + 2. / (h1 * h1) * (y2[1][n2] - y2[0][n2]) + 2. / h2 * my_data.back_boader(0.0, my_data) + 
					my_data.f_function(0.0, my_data.L2, my_data) + 2. / h1 * my_data.left_boarder(my_data.L2, my_data));



			progon_1(n1, n2, 0, y3, A2, B2, C2, F2, my_data);
		}

		if (my_data.condition[1] == 0)  //начальное гранусловие -- первого рода
			for (int i = 0; i <= n2; i++)	//конечное условие
				y3[n1][i] = my_data.right_boarder(i * h2, my_data);
		else							//начальное условие -- второго рода
		{
			F2[0] = -(1.0 - my_data.condition[2]) * my_data.front_boader(my_data.L1, my_data) +
				my_data.condition[2] * (2. / tau * y2[n1][0] + 2. / (h1 * h1) * (y2[n1 - 1][0] - y2[n1][0]) + 2. / h2 * my_data.front_boader(my_data.L1, my_data) + 
					my_data.f_function(my_data.L1, 0.0, my_data) + 2. / h1 * my_data.right_boarder(0.0, my_data));

			for (int i = 1; i < n2; i++)
				F2[i] = 2. / tau * y2[n1][i] + 2. / (h1 * h1) * (y2[n1 - 1][i] - y2[n1][i]) + 2. / h1 * my_data.right_boarder(i * h2, my_data) + my_data.f_function(my_data.L1, i * h2, my_data);

			F2[n2] = -(1.0 - my_data.condition[3]) * my_data.back_boader(my_data.L1, my_data) + 
				my_data.condition[3] * (2. / tau * y2[n1][n2] + 2. / (h1 *h1) * (y2[n1 - 1][n2] - y2[n1][n2]) + 2. / h2 * my_data.back_boader(my_data.L1, my_data) + 
					my_data.f_function(my_data.L1, my_data.L2, my_data) + 2. / h1 * my_data.right_boarder(my_data.L2, my_data));

			progon_1(n1, n2, n1, y3, A2, B2, C2, F2, my_data);
		}

		for (int i = 1; i != n1; i++)	//ходим вдоль
		{
			F2[0] = -(1.0 - my_data.condition[2]) * my_data.front_boader(i * h1, my_data) + 
				my_data.condition[2] * (2. / tau * y2[i][0] + (y2[i + 1][0] - 2. * y2[i][0] + y2[i - 1][0]) / (h1 * h1)
				+ my_data.f_function(i * h1, 0., my_data) + 2. / h2 * my_data.front_boader(i * h1, my_data));

			for (int j = 1; j < n2; j++)	//(по предыдущему слою)
				F2[j] = 2. / tau * y2[i][j] + (y2[i + 1][j] - 2. * y2[i][j] + y2[i - 1][j]) / (h1 * h1) + my_data.f_function(i * h1, j * h2, my_data);

			F2[n2] = -(1.0 - my_data.condition[3]) * my_data.back_boader(i * h1, my_data) + 
				my_data.condition[3] * (2. / tau * y2[i][n2] + (y2[i + 1][n2] - 2. * y2[i][n2] + y2[i - 1][n2]) / (h1 * h1)
				+ my_data.f_function(i * h1, my_data.L2, my_data) + 2. / h2 * my_data.back_boader(i * h1, my_data));



			progon_1(n1, n2, i, y3, A2, B2, C2, F2, my_data);
		}
		
//		max = max_el(y1, y3, n1, n2);
		max = max_err(y3, n1, n2, h1, h2);

		for (int i = 0; i <= n1; i++)
			for (int j = 0; j <= n2; j++)
				y1[i][j] = y3[i][j];
	//	fout_err << max << '\n';
		it++;
	//	std::cout << t << " is over max is " << max << '\n';
	}//while (max > 0.1);

		for (int i = 0; i <= n1; i++)
			for (int j = 0; j <= n2; j++)
				fout_err << i * h1 << '\t' << j * h2 << '\t' << fabs(y3[i][j] - test(i * h1, j * h2)) << '\n';

		for (int i = 0; i <= n1; i++)
			for (int j = 0; j <= n2; j++)
				fout << i * h1 << '\t' << j * h2 << '\t' << y3[i][j] << '\n';
	fout.close();
	fout_err.close();
	std::cout << "Your file is ready" << std::endl;
	return max;
}
