#ifndef NEWTON_METHODS_H
#define NEWTON_METHODS_H
#include "MyVector.h"
#include <functional>

/**
* \brief			Newton method realisation 2.
* \param[in] func	Function.
* \param[in] driv	Derivative.
* \param[in] l		First side of segment.
* \param[in] r		Second side of segment.
* \param[in] eps	Error for solution.
* \return			Root from this segment.
*/
double Newton(const std::function<double(double x)>& func, 
	const std::function<double(double x)>& driv, double l, double r, double eps);

/**
 * \brief			Function for finding absolutely minimal root.
 * \param[in] func	Function.
 * \param[in] driv	Derivative.
 * \param[in] eps	Error for solution.
 * \return			Minimal root by absolute value(not zero).
 */
double AbsoluteLessRoot(const std::function<double(double x)>& func, 
	const std::function<double(double x)>& driv, double eps);

/**
 * \brief			First implimentation of Newton method.
 * \param[in] func	Function for calculating functiond in point.
 * \param[in] driv	Function for calculating Jakobi matrix in point.
 * \param[in] fx	First x for calculations. 
 * \param[in] eps	Error for calculations.
 * \param[in] show	Flag if details should be shown.
 * \return			Answer with needed error.
 */
std::tuple<MyVector,int,int,int> NewtonES1(
	const std::function<MyVector(const MyVector& x)>& func,
	const std::function<Matrix(const MyVector& x)>& driv,
	MyVector fx, double eps, bool show = true);

/**
 * \brief			First implimentation of Newton method.
 * \param[in] func	Function for calculating functiond in point.
 * \param[in] driv	Function for calculating Jakobi matrix in point.
 * \param[in] fx	First x for calculations.
 * \param[in] eps	Error for calculations.
 * \param[in] show	Flag if details should be shown.
 * \return			Answer with needed error.
 */
std::tuple<MyVector, int, int, int> NewtonES2(
	const std::function<MyVector(const MyVector& x)>& func,
	const std::function<Matrix(const MyVector& x)>& driv, 
	MyVector fx, double eps, bool show = true);

/**
 * \brief			First implimentation of Newton method.
 * \param[in] func	Function for calculating functiond in point.
 * \param[in] driv	Function for calculating Jakobi matrix in point.
 * \param[in] fx	First x for calculations.
 * \param[in] k		Number for recalculating Jakobi matrix.
 * \param[in] eps	Error for calculations.
 * \param[in] show	Flag if details should be shown.
 * \return			Answer with needed error.
 */
std::tuple<MyVector, int, int, int> NewtonES3(
	const std::function<MyVector(const MyVector& x)>& func,
	const std::function<Matrix(const MyVector& x)>& driv,
	MyVector fx, int k, double eps, bool show = true);

/**
* \brief			First implimentation of Newton method.
* \param[in] func	Function for calculating functiond in point.
* \param[in] driv	Function for calculating Jakobi matrix in point.
* \param[in] fx	First x for calculations.
* \param[in] k		Number for recalculating Jakobi matrix.
* \param[in] eps	Error for calculations.
* \param[in] show	Flag if details should be shown.
* \return			Answer with needed error.
*/
std::tuple<MyVector, int, int, int> NewtonES4(
	const std::function<MyVector(const MyVector& x)>& func,
	const std::function<Matrix(const MyVector& x)>& driv,
	MyVector fx, int k, double eps, bool show = true);


/**
 * \brief			First implimentation of Newton method.
 * \param[in] func	Function for calculating functiond in point.
 * \param[in] driv	Function for calculating Jakobi matrix in point.
 * \param[in] fx	First x for calculations.
 * \param[in] eps	Error for calculations.
 * \param[in] show	Flag if details should be shown.
 * \return			Answer with needed error.
 */
std::tuple<MyVector, int, int, int> NewtonESh(
	const std::function<MyVector(const MyVector& x)>& func,
	const std::function<Matrix(const MyVector& x)>& driv,
	MyVector fx, double eps, bool show = true);

#endif