#include<iostream>
#include<math.h>
#include <fstream>
#include<vector>
#include <algorithm>

#define _USE_MATH_DEFINES

double p(double x)
{
    return cosh(x);
}

double q(double x)
{
    return sinh(x);
}

double f(double x)
{
    return cosh(x) + x * sinh(x);
}

double sol(double x)
{
    return exp(-sinh(x)) + x;
}

int main()
{
    // начальные данные
    double h = 0.05;
    double a = 0;
    double b = 1;

    double g_a = 0;
    double n_a = 1;
    double y_a = 0;

    double g_b = 6;
    double n_b = 1;
    double y_b = 8.3761;

    // Равномерная сетка
    int N = int((b - a)/h);
    double x_lin[N + 1] {0};
    for(int i = 0; i < N+1; i++)
        x_lin[i] = a + h*i;

    // Коэффициенты матрицы системы
    double A[N + 1] {0};
    A[0] = g_a - n_a*(2 - q(x_lin[0])*h*h)/(2 - h*p(x_lin[0]))/h;
    for(int i = 1; i < N; i++)
        A[i] = 1/h/h - p(x_lin[i])/2/h;
    A[N] = 0;

    double B[N + 1] {0};
    B[0] = n_a*(1 + (2 + p(x_lin[0])*h)/(2 - p(x_lin[0])*h))/2/h;
    for(int i = 1; i < N; i++)
        B[i] = -2/h/h + q(x_lin[i]);
    B[N] = -n_b*(1 + (2 - p(x_lin[N])*h)/(2 + p(x_lin[N])*h))/2/h;

    double C[N + 1] {0};
    C[0] = 0;
    for(int i = 1; i < N; i++)
        C[i] = 1/h/h + p(x_lin[i])/2/h;
    C[N] = g_b + n_b*(2 - q(x_lin[N])*h*h)/h/(2 + h*p(x_lin[N]));

    double D[N + 1] {0};
    D[0] = y_a + n_a*h*f(x_lin[0])/(2 - h*p(x_lin[0]));
    for(int i = 1; i < N; i++)
        D[i] = f(x_lin[i]);
    D[N] = y_b - n_b*h*f(x_lin[N])/(2 + h*p(x_lin[N]));

    // Прямой ход прогонки
    double a_l[N] {0};
    double b_l[N] {0};
    a_l[0] = -B[0]/A[0];
    b_l[0] = D[0]/A[0];
    for(int i = 1; i < N; i++)
    {
        b_l[i] = (D[i] - A[i]*b_l[i-1])/(A[i]*a_l[i-1] + B[i]);
        a_l[i] = -C[i]/(A[i]*a_l[i-1] + B[i]);
    }

    // Обратный ход прогонки
    std::vector<double> res {};
    res.push_back((D[N] - B[N]*b_l[N-1])/(B[N]*a_l[N-1] + C[N]));
    for(int i = N-1; i >= 0; i--)
        res.push_back(res.back() * a_l[i] + b_l[i]);
    std::reverse(res.begin(), res.end());

    // запись результатов в файл
    std::ofstream file;
    file.open("progonka.csv");
    file << "x,sol,progon_2\n";
    for(int i = 0; i < N+1; i++)
    {
        file << x_lin[i] << ",";
        file << sol(x_lin[i]) << ",";
        file << res[i] << "\n";
    }
    file.close();
    std::cout << "OK";
}

