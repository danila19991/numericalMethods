#ifndef INTEGRALS_H
#define INTEGRALS_H
#include <utility>

/**
 * \brief			Found integral of function.
 * \param[in] func	Function for calculating integral.
 * \param[in] l		Left side of segment.
 * \param[in] r		Right ide of segment.
 * \param[in] len	Maximal length of step.
 * \return			Value of this intageal.
 */
double RectanglesMethod(double(*func)(double x), double l, double r, double len);

/**
 * \brief				Found integral of function.
 * \param[in] func		Function for calculating integral.
 * \param[in] l			Left side of segment.
 * \param[in] r			Right ide of segment.
 * \param[in] b			Parametr of weight function.
 * \param[in] k			Parametr of weight function.
 * \param[in] maxDer	Absolutly maximal value of n's derivative.
 * \return				Value of this intageal and error.
*/
std::pair<double, double> IQAIntagral(double(*func)(double x), double l, double r, double b,
	double k, double maxDer);

/**
 * \brief			Found integral of function.
 * \param[in] func	Function for calculating integral.
 * \param[in] l		Left side of segment.
 * \param[in] r		Right side of segment.
 * \param[in] b		Parametr of weight function.
 * \param[in] k		Parametr of weight function.
 * \param[in] eps	Error of naswer
 * \return			Value of this intageal and error.
 */
double IntagralNewtonKotsAcurate(double (*func)(double x), double l, double r, double b, double k,
                          double eps = 1e-6);

/**
* \brief			Found integral of function.
* \param[in] func	Function for calculating integral.
* \param[in] l		Left side of segment.
* \param[in] r		Right side of segment.
* \param[in] b		Parametr of weight function.
* \param[in] k		Parametr of weight function.
* \param[in] len	Lehgth of step.
* \return			Value of this intageal and error.
*/
double IntagralNewtonKotsPartial(double (*func)(double x), double l, double r, double b, double k,
                              double len = 0.1);

/**
* \brief			Found integral of function.
* \param[in] func	Function for calculating integral.
* \param[in] l		Left side of segment.
* \param[in] r		Right side of segment.
* \param[in] b		Parametr of weight function.
* \param[in] k		Parametr of weight function.
* \param[in] len	Lehgth of step.
* \return			Value of this intageal and error.
*/
double GausIntagralPartial(double(*func)(double x), double l, double r, double b, double k, 
	double len);

/**
* \brief				Found integral of function.
* \param[in] func		Function for calculating integral.
* \param[in] l			Left side of segment.
* \param[in] r			Right ide of segment.
* \param[in] b			Parametr of weight function.
* \param[in] k			Parametr of weight function.
* \param[in] maxDer		Absolutely maximal value of 2n derivative.
* \return				Value of this intageal and error.
*/
std::pair<double, double> GausIntagral(double(*func)(double x), double l, double r, double b, 
	double k, double maxDer);

/**
* \brief			Found integral of function.
* \param[in] func	Function for calculating integral.
* \param[in] l		Left side of segment.
* \param[in] r		Right side of segment.
* \param[in] b		Parametr of weight function.
* \param[in] k		Parametr of weight function.
* \param[in] eps	Error of taking integral.
* \return			Value of this intageal and error.
*/
double GausIntagralAcurate(double(*func)(double x), double l, double r, double b, double k, 
	double eps);
#endif // INTEGRALS_H