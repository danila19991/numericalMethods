#ifndef BFS_H
#define BFS_H
#include "MyVector.h"

/**
 * \brief		Big function system.
 * \param[in] x	Point for calculating result of function.
 * \return		Result of function.
 */
MyVector bfs(const MyVector& x);

/**
 * \brief		Function for calculating Jakobi matrix in point.
 * \param[in] x	Point for calculating Jakobi matrix in point.
 * \return		Jakobi matrix in point.
 */
Matrix bfd(const MyVector& x);

#endif // BFS_H
