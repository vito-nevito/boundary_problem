#ifndef SHOOTING_H
#define SHOOTING_H
#include <vector>
#include "RungeKutta.h"
class Shooting
{
protected:
	func f1, f2, f3;
	double alpha1, betha1, alpha2, betha2;
	double gamma1, gamma2;
	double a, b;
	double step;
	double precision;
	std::vector<double> u{}, v{};
public:
	Shooting(func f, func g, func k, double a1, double b1, double g1, double a2, double b2, double g2, double l, double r, double h, double eps) :
		f1(f), f2(g), f3(k),
		alpha1(a1), betha1(b1), gamma1(g1),
		alpha2(a2), betha2(b2), gamma2(g2),
		a(l), b(r),
		step(h),
		precision(eps)
	{}
	double getXi(int i);
	double getLambdai(int i);
	double rightCond(double lambda);
	std::vector<double> solve();
};

#endif // !SHOOTING_H


