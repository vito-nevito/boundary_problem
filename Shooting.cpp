#include "Shooting.h"
#include <iostream>
double Shooting::getLambdai(int i)
{
	return 0 + pow(-1, i) * 10 * i;
}
double Shooting::rightCond(double x)
{
	RungeKutta system(*f1, *f2, *f3, x, (gamma1 - alpha1 * x) / betha1, a, b, step);
	std::vector<double> u1 = system.solve();
	std::vector<double> v1 = system.getV();
	return alpha2 * u1[u1.size() - 1] + betha2 * v1[v1.size() - 1] - gamma2;
}
std::vector<double> Shooting::solve()
{
	int i(1);
	while (rightCond(getLambdai(0)) * rightCond(getLambdai(i)) > 0)
	{
		i++;
	}
	double al, bl;
	if (i % 2 == 1)
	{
		al = getLambdai(i);
		bl = 0;
	}
	else
	{
		al = 0;
		bl = getLambdai(i);
	}
	double xnew((al + bl) / 2), xprev(al);
	while ( not ( ( abs(xnew - xprev) <= precision ) or (abs( rightCond(xnew) ) <= precision) ) )
	{
		xprev = xnew;
		if (rightCond(xnew) * rightCond(al) < 0)
		{
			bl = xnew;
		}
		else if (rightCond(xnew) * rightCond(bl) < 0)
		{
			al = xnew;
		}
		xnew = (al + bl) / 2;
	}
	RungeKutta goodSys(*f1, *f2, *f3, xnew, (gamma1 - alpha1 * xnew) / betha1, a, b, step);
	return goodSys.solve();
}