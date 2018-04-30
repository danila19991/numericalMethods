#ifndef NEWTON_METHODS_H
#define NEWTON_METHODS_H

double Newton1(double(*func)(double x),double(*driv)(double x), double l, double r, double eps);

double AbsoluteLessRoot1(double(*func)(double x), double(*driv)(double x),double eps);

double Newton2(double(*func)(double x), double(*driv)(double x), double l, double r, double eps);

double AbsoluteLessRoot2(double(*func)(double x), double(*driv)(double x), double eps);

#endif