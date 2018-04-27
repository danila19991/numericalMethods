#ifndef ZEIDEL_H
#define ZEIDEL_H

#include "Matrix.h"
#include "MyVector.h"

/**
 * \brief			Jakobi method for solving equation system.
 * \param[in] a		Matrix of equation system.
 * \param[in] b		Vector-column of equation system.
 * \param[out] flag	Flag if this equation was solved in 100'000 steps.
 * \param[in] eps	Error of answer.
 * \return			Answer for equation system.
 */
MyVector Jakobi(Matrix a, MyVector b,bool& flag, double eps = 1e-6);

/**
* \brief			Zeidel method for solving equation system.
* \param[in] a		Matrix of equation system.
* \param[in] b		Vector-column of equation system.
* \param[out] flag	Flag if this equation was solved in 100'000 steps.
* \param[in] eps	Error of answer.
* \return			Answer for equation system.
*/
MyVector Zeidel(Matrix a, MyVector b,bool& flag, double eps = 1e-6);

#endif // ZEIDEL_H