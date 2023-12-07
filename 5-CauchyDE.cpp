#include <iostream>
#include "Euler.h"
#include "RungeKutta.h"
#include "AdamsBashforth.h"
#include <cmath>
#include <fstream>
#include "Shooting.h"

double f1(double x) { return cosh(x); }

double f2(double x) { return sinh(x); }

double f3(double x) { return cosh(x) + x * sinh(x); }

double u0(double x) { return exp(-sinh(x)) + x; }

int main()
{
	
	double a(0), b(1);
	std::vector<double> h = { 0.01, 0.02, 0.05, 0.1 };
	std::vector<std::vector<double>> ans;
	for (int i(0); i < h.size(); i++)
	{
		Shooting sys(&f1, &f2, &f3, 0, 1, 0, 6, 1, 8.3761, a, b, h[i], 1e-9);
		ans.insert(ans.begin() + i, sys.solve());
	}

	std::ofstream out11;
	std::ofstream out12;
	std::ofstream out13;
	std::ofstream out14;
	out11.open("shooting1.txt");
	out12.open("shooting2.txt");
	out13.open("shooting3.txt");
	out14.open("shooting4.txt");

	for (int i(0); i < ans[0].size(); i++)
	{
		out11 << a + h[0] * i << "\t" << ans[0][i] << "\n";
	}
	for (int i(0); i < ans[1].size(); i++)
	{
		out12 << a + h[1] * i << "\t" << ans[1][i] << "\n";
	}
	for (int i(0); i < ans[2].size(); i++)
	{
		out13 << a + h[2] * i << "\t" << ans[2][i] << "\n";
	}
	for (int i(0); i < ans[3].size(); i++)
	{
		out14 << a + h[3] * i << "\t" << ans[3][i] << "\n";
	}

	return 0;
}