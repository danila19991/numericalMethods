#include "Calculations.h"

#include<algorithm>
#include <cassert>

const double PI = acos(-1.);

double mysqrt(const double x, const double eps) noexcept
{
	assert(x >= 0.);
	double ansPrev = std::max(x, 1.);
	double ans = (ansPrev + x / ansPrev) / 2;
	while (abs(ans - ansPrev) >= eps) 
	{
		ansPrev = ans;
		ans = (ansPrev + x / ansPrev) / 2;
	}
	return ans;
}

double mysin(double x, const double eps) noexcept 
{
	const auto k = static_cast<int>(x / (PI / 2));
	x -= k * (PI/2);
	if (k%4 == 1) 
	{
		return mycos(x, eps);
	}
	if (k%4 == 2)
	{
		return -1.*mysin(x, eps);
	}
	if(k%4 == 3)
	{
		return -1.*mycos(x, eps);
	}
	double ans = x, del = x;
	for (int i = 1;abs(del) >= eps; i++) 
	{
		del = del * x * x / (2 * i) / (2 * i + 1);
		if (i & 1)
		{
			ans -= del;
		}
		else
		{
			ans += del;
		}
	}
	return ans;
}

double mycos(double x, const double eps) noexcept
{
	const auto k = static_cast<int>(x / (PI / 2));
	x -= k * (PI / 2);
	if (k % 4 == 1)
	{
		return mysin(x, eps);
	}
	if (k % 4 == 2)
	{
		return -1.*mycos(x, eps);
	}
	if (k % 4 == 3)
	{
		return -1.*mysin(x, eps);
	}
	double ans = 1, del = 1;
	for (int i = 1;abs(del) >= eps; i++) 
	{
		del = del * x * x / (2 * i) / (2 * i - 1);
		if (i & 1)
		{
			ans -= del;
		}
		else
		{
			ans += del;
		}
	}
	return ans;
}

double myatan(const double x, const double eps) noexcept 
{
	double ans,a;
	if(x > 1. - eps || x < -1. + eps)
	{
		a = 1 / x;
		if (x > 1. - eps)
		{
			ans = PI / 2 - a;
		}
		else
		{
			ans = -PI / 2 - a;
		}
		for (int i = 1;abs(a / (2 * i - 1)) >= eps;i++)
		{
			a = a / x / x;
			if (i & 1)
			{
				ans += a / (2 * i + 1);
			}
			else
			{
				ans -= a / (2 * i + 1);
			}
		}
	}
	else
	{
		ans = a = x;
		for(int i=1; abs(a / (2 * i - 1))>=eps;i++)
		{
			a = a * x * x;
			if (i & 1)
			{
				ans -= a / (2 * i + 1);
			}
			else
			{
				ans += a / (2 * i + 1);
			}
		}
	}
	return ans;
}

double myexp(const double x, const double eps) noexcept
{
	if(x < 0)
	{
		return 1 / myexp(-x, eps);
	}
	double ans=1., a = 1.;
	for(int i=1;abs(a)>=eps;i++)
	{
		a = a * x / i;
		ans += a;
	}
	return ans;
}

