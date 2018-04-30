#ifndef NEWTON_METHODS_H
#define NEWTON_METHODS_H

/**
 * \brief			Newton method realisation 1.
 * \param[in] func	Function.
 * \param[in] driv	Derivative.
 * \param[in] l		First side of segment.
 * \param[in] r		Second side of segment.
 * \param[in] eps	Error for solution.
 * \return			Root from this segment.
 */
double Newton1(double(*func)(double x),double(*driv)(double x), double l, double r, double eps);

/**
 * \brief			Function for finding absolutely minimal root.
 * \param[in] func	Function.
 * \param[in] driv	Derivative.
 * \param[in] eps	Error for solution.
 * \return			Minimal root by absolute value(not zero).
 */
double AbsoluteLessRoot1(double(*func)(double x), double(*driv)(double x),double eps);

/**
 * \brief			Newton method realisation 2.
 * \param[in] func	Function.
 * \param[in] driv	Derivative.
 * \param[in] l		First side of segment.
 * \param[in] r		Second side of segment.
 * \param[in] eps	Error for solution.
 * \return			Root from this segment.
 */
double Newton2(double(*func)(double x), double(*driv)(double x), double l, double r, double eps);

/**
 * \brief			Function for finding absolutely minimal root.
 * \param[in] func	Function.
 * \param[in] driv	Derivative.
 * \param[in] eps	Error for solution.
 * \return			Minimal root by absolute value(not zero).
 */
double AbsoluteLessRoot2(double(*func)(double x), double(*driv)(double x), double eps);

#endif