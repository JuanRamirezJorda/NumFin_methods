#pragma once

#ifndef BASICFUNC
#define BASICFUNC


#include <iostream>
#include <cmath>
#include <stdexcept>
#include <limits>

int recursivelySumBetween1AndN(int n);

double normcdf(double x);

double hornerFunction(double x, double a0, double a1);

double hornerFunction(double x, double a0, double a1, double a2);

double hornerFunction(double x, double a0, double a1, double a2, double a3);

double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4);

double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5);

double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6);

double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7);

double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8);

double normcdfHorner(double x);

double norminv(double x);

double blackScholesCallPrice(double strike, double timeToMaturity, double spot, double volatility, double riskFreeRate);

double integrateSin(double a, double b, int N);

double infiniteIntegral(double x);

double norminvNoCheck(double x);

double norminv(double x, bool checkRange);




#endif // !BASICFUNC
