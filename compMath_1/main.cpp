#include <iomanip>
#include <fstream>
#include <conio.h>
#include <iostream>
#include <chrono>

#include "PLUMatrix.h"
#include "QRMatrix.h"
#include "IterationES.h"
#include <cassert>
#include "NewtonMethods.h"
#include "BFS.h"
#include <tuple>

void t3p1(double(*func)(double x), double(*driv)(double x), double eps)
{
	assert(eps > 0);

	const double ans1 = AbsoluteLessRoot1(func, driv, eps);

	std::cout << ans1 << ' ' << func(ans1) << "\n\n";

	const double ans2 = AbsoluteLessRoot2(func, driv, eps);

	std::cout << ans2 << ' ' << func(ans2) << "\n\n";
}

double f(const double x)
{
	return x * x * x - exp(x) + 1;
}

double d(const double x)
{
	return 3 * x * x - exp(x);
}

void t3p2a(MyVector(*func)(const MyVector& x), Matrix(*driv)(const MyVector& x), const MyVector& x,
	double eps, bool show = true)
{
	MyVector x1;
	int time, oper, calc;

	std::tie(x1, time, oper, calc) = NewtonES1(func, driv, x, eps, show);

	std::cout << "was " << oper << " operations and " << calc << " calcualtions\n";
	std::cout << "it needs " << time << " milliseconds\n";
	if(show)
		std::cout << '\n' << func(x1)<<'\n';
	std::cout << '\n';
}

void t3p2b(MyVector(*func)(const MyVector& x), Matrix(*driv)(const MyVector& x), const MyVector& x,
	double eps, bool show = true)
{
	MyVector x1;
	int time, oper, calc;

	std::tie(x1, time, oper, calc) = NewtonES2(func, driv, x, eps, show);

	std::cout << "was " << oper << " operations and " << calc << " calcualtions\n";
	std::cout << "it needs " << time << " milliseconds\n";
	if (show)
		std::cout << '\n' << func(x1) << '\n';
	std::cout << '\n';
}

void t3p2cp(MyVector(*func)(const MyVector& x), Matrix(*driv)(const MyVector& x), const MyVector& x,
	int mx, double eps, bool show = true)
{
	MyVector x1;
	int time, oper, calc;
	int mn = INT_MAX, ans=-1,at=-1,ao=-1,ap=-1;

	for (int k = 1;k <= mx;++k) {
		std::tie(x1, time, oper, calc) = NewtonES3(func, driv, x, k, eps, show);
		if(show)
			std::cout << "for k = " << k << " : " << time << ", " << oper << ", " << calc << "\n";
		if(oper < mn)
		{
			mn = oper;
			ans = k;
			at = time;
			ao = oper;
			ap = calc;
		}
	}

	std::cout << "the best in k = " << ans << " : " << at << ", " << ao << ", " << ap << "\n\n";
}

void t3p2c(MyVector(*func)(const MyVector& x), Matrix(*driv)(const MyVector& x), const MyVector& x,
	int k ,double eps, bool show = true)
{
	MyVector x1;
	int time, oper, calc;

	std::tie(x1, time, oper, calc) = NewtonES3(func, driv, x, k, eps, show);

	std::cout << "was " << oper << " operations and " << calc << " calcualtions\n";
	std::cout << "it needs " << time << " milliseconds\n";
	if (show)
		std::cout << '\n' << func(x1) << '\n';
	std::cout << '\n';
}

void t3p2dp(MyVector(*func)(const MyVector& x), Matrix(*driv)(const MyVector& x), const MyVector& x,
	int mx, double eps, bool show = true)
{
	MyVector x1;
	int time, oper, calc;
	int mn = INT_MAX, ans = -1, at = -1, ao = -1, ap = -1;

	for (int k = 1;k <= mx;++k) {
		std::tie(x1, time, oper, calc) = NewtonES4(func, driv, x, k, eps, show);
		if(show)
			std::cout << "for k = " << k << " : " << time << ", " << oper << ", " << calc << "\n";
		if (oper < mn)
		{
			mn = oper;
			ans = k;
			at = time;
			ao = oper;
			ap = calc;
		}
	}

	std::cout << "the best in k = " << ans << " : " << at << ", " << ao << ", " << ap << "\n\n";
}

void t3p2d(MyVector(*func)(const MyVector& x), Matrix(*driv)(const MyVector& x), const MyVector& x,
	int k, double eps, bool show = true)
{
	MyVector x1;
	int time, oper, calc;

	std::tie(x1, time, oper, calc) = NewtonES4(func, driv, x, k, eps, show);

	std::cout << "was " << oper << " operations and " << calc << " calcualtions\n";
	std::cout << "it needs " << time << " milliseconds\n";
	if (show)
		std::cout << '\n' << func(x1) << '\n';
	std::cout << '\n';
}

void  t3p2e(MyVector(*func)(const MyVector& x), Matrix(*driv)(const MyVector& x), double eps,
	bool show = true)
{
	const MyVector x({ -0.2, 0.5, 1.5, -1.0, -0.5, 1.5, 0.5, -0.5, 1.5, -1.5 });

	t3p2a(func, driv, x, eps, show);

	t3p2b(func, driv, x, eps, show);

	t3p2cp(func, driv, x, 17, eps, show);

	t3p2dp(func, driv, x, 52, eps, show);
}

int main()
{
	std::cout.precision(10);
	std::cout << std::fixed;

	const MyVector x({ 0.5, 0.5, 1.5, -1.0, -0.5, 1.5, 0.5, -0.5, 1.5, -1.5 });

	//t3p1(f, d, 1e-4);

	t3p2a(bfs, bfd, x, 1e-9);

	//t3p2b(bfs, bfd, x, 1e-9);

	//t3p2cp(bfs, bfd, x, 18, 1e-9);

	//t3p2c(bfs, bfd, x, 3, 1e-9);

	//t3p2dp(bfs, bfd, x, 76, 1e-9);

	//t3p2d(bfs, bfd, x, 13, 1e-9);

	//t3p2e(bfs, bfd, 1e-9);
	
	_getch();

	return 0;
}
