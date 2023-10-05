#include "basicfunc.h"


int recursivelySumBetween1AndN(int n) {
	if (n == 1) {
		return 1;
	}
	return n + recursivelySumBetween1AndN(n - 1);
}


const double root2Pi = sqrt(2.0 * 3.141592653589793);
double normcdf(double x) {
	if (x < 0) {
		return 1 - normcdf(-x);
	}
	double k = 1 / (1 + 0.2316419 * x);
	double poly = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937
		+ k * (-1.821255978 + 1.330274429 * k))));
	double approx = 1.0 - 1.0 / root2Pi * exp(-0.5 * x * x) * poly;
	return approx;
}


double hornerFunction(double x, double a0, double a1) {
	return a0 + x * a1;
}

double hornerFunction(double x, double a0, double a1, double a2) {
	return a0 + x * hornerFunction(x, a1, a2);
}

double hornerFunction(double x, double a0, double a1, double a2, double a3) {
	return a0 + x * hornerFunction(x, a1, a2, a3);
}

double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4) {
	return a0 + x * hornerFunction(x, a1, a2, a3, a4);
}

double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4,
	double a5) {
	return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5);
}

double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4,
	double a5, double a6) {
	return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5, a6);
}

double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4,
	double a5, double a6, double a7) {
	return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5, a6, a7);
}

double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4,
	double a5, double a6, double a7, double a8) {
	return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5, a6, a7, a8);
}

/**
 *  Arguably this is a little easier to read than the original normcdf
 *  function as it makes the use of horner's method obvious.
 */
double normcdfHorner(double x) {
	if (x <= 0) {
		return 1 - normcdf(-x);
	}
	double k = 1 / (1 + 0.2316419 * x);
	double poly = hornerFunction(k,
		0.0, 0.319381530, -0.356563782,
		1.781477937, -1.821255978, 1.330274429);
	double approx = 1.0 - 1.0 / root2Pi * exp(-0.5 * x * x) * poly;
	return approx;
}

/**
 *   Data sample to test Moro & Horner's algorithms
 */
const double a0 = 2.50662823884;
const double a1 = -18.61500062529;
const double a2 = 41.39119773534;
const double a3 = -25.44106049637;
const double b1 = -8.47351093090;
const double b2 = 23.08336743743;
const double b3 = -21.06224101826;
const double b4 = 3.13082909833;
const double c0 = 0.3374754822726147;
const double c1 = 0.9761690190917186;
const double c2 = 0.1607979714918209;
const double c3 = 0.0276438810333863;
const double c4 = 0.0038405729373609;
const double c5 = 0.0003951896511919;
const double c6 = 0.0000321767881768;
const double c7 = 0.0000002888167364;
const double c8 = 0.0000003960315187;


double blackScholesCallPrice(double strike, double timeToMaturity,
	double spot, double volatility, double riskFreeRate) {
	double numerator = log(spot / strike)
		+ (riskFreeRate + volatility * volatility * 0.5) * timeToMaturity;
	double denominator = volatility * sqrt(timeToMaturity);
	double d1 = numerator / denominator;
	double d2 = d1 - denominator;
	double t1 = normcdf(d1) * spot;
	double t2 = normcdf(d2) * strike * exp(-riskFreeRate * timeToMaturity);
	return t1 - t2;
}

double integrateSin(double a, double b, int N) {
	double h = (b - a) / N;
	double total = 0.0;
	for (int i = 0; i < N; i++) {
		double x = (i + 0.5) * h + a;
		double f = sin(x);
		total += f;
	}
	return total / N;
}


double infiniteIntegral(double x) {
	// we perform the substitution
	// x + 1 - 1/s;
	// to change the infinite integral to an integral between 0 and 1
	double a = 0;
	double b = 1;
	int N = 1000;
	double h = (b - a) / N;
	double total = 0.0;
	for (int i = 0; i < N; i++) {
		double s = (i + 0.5) * h + a;
		double t = x + 1 - 1 / s;
		double f = pow(s, -2) * exp(-0.5 * t * t);
		total += f;
	}
	return total / N;
}



double norminvNoCheck(double x) {
	// We use Moro's algorithm
	double y = x - 0.5;
	if (y<0.42 && y>-0.42) {
		double r = y * y;
		return y * hornerFunction(r, a0, a1, a2, a3) / hornerFunction(r, 1.0, b1, b2, b3, b4);
	}
	else {
		double r;
		if (y < 0.0) {
			r = x;
		}
		else {
			r = 1.0 - x;
		}
		double s = log(-log(r));
		double t = hornerFunction(s, c0, c1, c2, c3, c4, c5, c6, c7, c8);
		if (x > 0.5) {
			return t;
		}
		else {
			return -t;
		}
	}
}

double norminv(double x, bool checkRange = true) {
	// Note the #include <stdexcept> at the top of the file
	if (checkRange && (x < 0 || x>1.0)) {
		throw std::logic_error("Parameter x is out of range for norminv. It should be between 0 and 1");
	} 
	else {
		return norminvNoCheck(x);
	}
}



static void testNormCdf() {
	// test bounds
	ASSERT(normcdf(0.3) > 0);
	ASSERT(normcdf(0.3) < 1);
	// test extreme values
	ASSERT_APPROX_EQUAL(normcdf(-1e10), 0, 0.001);
	ASSERT_APPROX_EQUAL(normcdf(1e10), 1.0, 0.001);
	// test increasing
	ASSERT(normcdf(0.3) < normcdf(0.5));
	// test symmetry
	ASSERT_APPROX_EQUAL(normcdf(0.3),
		1 - normcdf(-0.3), 0.0001);
	ASSERT_APPROX_EQUAL(normcdf(0.0), 0.5, 0.0001);
	// test inverse
	ASSERT_APPROX_EQUAL(normcdf(norminv(0.3)),
		0.3, 0.0001);
	// test well known value
	ASSERT_APPROX_EQUAL(normcdf(1.96), 0.975, 0.001);
}

static void testNormInv() {
	ASSERT_APPROX_EQUAL(norminv(0.975), 1.96, 0.01);
}

void testSample() {
	TEST(testNormInv);
	TEST(testNormCdf);
}