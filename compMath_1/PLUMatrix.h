#ifndef PLU_MATRIX
#define PLU_MATRIX

#include "Matrix.h"
#include "MyVector.h"

/**
 * \brief Class of PLU decomposition.
 */
class PLUMatrix
{
	/**
	 * \brief Data storrage.
	 */
	Matrix _a;

	/**
	 * \brief Columns permutation.
	 */
	std::vector<int> _p;

	/**
	 * \brief Rows permutation.
	 */
	std::vector<int> _q;

public:
	/**
	 * \brief		Constructor for PLU matrix decomposition.
	 * \param[in] a Matrix for decomposition.
	 */
	PLUMatrix(Matrix a);

	/**
	 * \brief	Gatter for primal matrix.
	 * \return	Primal matrix.
	 */
	Matrix getA() const noexcept;

	/**
	 * \brief	Calculate determinant of matrix.
	 * \return	Determinant of matrix.
	 */
	double getDeterminant()const noexcept;

	/**
	 * \brief	Calcualte rank of matrix.
	 * \return	Rank of matrix.
	 */
	size_t getRank() const noexcept;

	/**
	 * \brief	Make obrat matrix.
	 * \return	Obrat matrix.
	 */
	Matrix obrat() const noexcept;

	/**
	 * \brief	Check if equation system has solution.
	 * \param b Vector-column of free members.
	 * \return	Answer for equation system.
 	 */
	bool hasSolution(MyVector b) const noexcept;

	/**
	 * \brief	Solve equation system.
	 * \param b Vector-column of free members.
	 * \return	Answer for equation system.
	 */
	MyVector solve(MyVector b) const noexcept;

	/**
	 * \brief	Return condition number of matrix.
	 * \return	Result.
	 */
	double getConditionNumber() const noexcept;

	/**
	 * \brief			Print all local fields in stream.
	 * \param[in] out	Stream for outputting.
	 */
	void printPLU(std::ostream& out) const noexcept;

	/**
	* \brief				Get random square posistive difinite matrix.
	* \param[in] size		Number of rows in matrix.
	* \param[in] minElem	Minimum value of element.
	* \param[in] maxElem	Maximum value of element.
	* \return				Random matrix.
	*/
	static Matrix getRandomPositiveDifinite(size_t size, double minElem, double maxElem);

	/**
	 * \brief			Operator for outputting matrix to stream.
	 * \param[in] out	Stream for outputting.
	 * \param[in] m		Matrix for outputting.
	 * \return			Stream for outputting.
	 */
	friend std::ostream& operator << (std::ostream& out, const PLUMatrix& m);

	/**
	* \brief			Operator for inputting matrix to stream.
	* \param[in] in		Stream for inputting.
	* \param[in] m		Matrix for inputting.
	* \return			Stream for inputting.
	*/
	friend std::istream& operator >> (std::istream& in, PLUMatrix& m);
};

//todo
//optimise swapping of columns and rows.

#endif // PLU_MATRIX