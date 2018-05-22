#ifndef DIFFERENTIAL_EQUATIONS
#define DIFFERENTIAL_EQUATIONS

#include "Matrix.h"
#include "MyVector.h"
#include <functional>

/**
 * \brief			Solving differntal Equation using modified Eiler metod.
 * \param[in] func	System of differential equations.
 * \param[in] y0	Starting condition in left side of segment.
 * \param[in] len	Length of step for calculating.
 * \param[in] l		Left side of segmant.
 * \param[in] r		Rieght side of segment.
 * \return			Set of points describing plot of solution.
 */
std::vector<std::pair<double,MyVector>> DifferentialEquationsOpp(
	const std::function<MyVector(double,const MyVector&)>& func,MyVector y0, double len,
	double l, double r);

/**
 * \brief			Solving differntal Equation using Runge-Kutt metod.
 * \param[in] func	System of differential equations.
 * \param[in] y0	Starting condition in left side of segment.
 * \param[in] len	Length of step for calculating.
 * \param[in] l		Left side of segmant.
 * \param[in] r		Rieght side of segment.
 * \param[in] a		Matrix with constants
 * \param[in] b		Vector with constants.
 * \param[in] c		Vector with constants.
 * \return			Set of points describing plot of solution.
 */
std::vector<std::pair<double, MyVector>> RungeKutt(
	const std::function<MyVector(double, const MyVector&)>& func, MyVector y0, double len,
	double l, double r, const Matrix& a, const MyVector& b, const MyVector& c);

/**
 * \brief			Solving differntal Equation using modified Eiler metod with auto-step.
 * \param[in] func	System of differential equations.
 * \param[in] y0	Starting condition in left side of segment.
 * \param[in] l		Left side of segmant.
 * \param[in] r		Rieght side of segment.
 * \param[in] rtol	Relative tolerance for one step.
 * \param[in] atol	Absolute tolerance for one step.
 * \param[in] mlen	Maximum len of step.
 * \return			Set of points describing plot of solution, Set of steps which was made, 
 *					Set of steps which may be made, set of calculated error on each step, 
 *					number of calculating func
 */
std::tuple<std::vector<std::pair<double, MyVector>>, std::vector<std::pair<double, double>>, 
	std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>, size_t>  
	DifferentialEquationsOppAuto(const std::function<MyVector(double, const MyVector&)>& func,
	MyVector y0, double l, double r, double rtol, double atol, double mlen);

/**
 * \brief			Solving differntal Equation using Runge-Kutt metod with auto-step.
 * \param[in] func	System of differential equations.
 * \param[in] y0	Starting condition in left side of segment.
 * \param[in] l		Left side of segmant.
 * \param[in] r		Rieght side of segment.
 * \param[in] rtol	Relative tolerance for one step.
 * \param[in] atol	Absolute tolerance for one step.
 * \param[in] mlen	Maximum len of step.
 * \param[in] a		Matrix with constants
 * \param[in] b		Vector with constants.
 * \param[in] c		Vector with constants.
 * \return			Set of points describing plot of solution, Set of steps which was made, 
 *					Set of steps which may be made, set of calculated error on each step, 
 *					number of calculating func
 */
std::tuple<std::vector<std::pair<double, MyVector>>, std::vector<std::pair<double, double>>,
	std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>, size_t> 
	RungeKuttAuto( const std::function<MyVector(double, const MyVector&)>& func, MyVector y0,
	double l, double r,double rtol, double atol, double mlen, const Matrix&  a, const MyVector& b,
	const MyVector& c);

#endif // DIFFERENTIAL_EQUATIONS