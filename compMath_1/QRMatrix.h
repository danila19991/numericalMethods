#ifndef QR_MATRIX
#define QR_MATRIX

#include "Matrix.h"
#include "MyVector.h"

/**
 * \brief Class for solving equation system using QR decomposition.
 */
class QRMatrix
{
	/**
	 * \brief Q matrix in decomposition.
	 */
	Matrix _q;

	/**
	 * \brief R matrix in decomposition.
	 */
	Matrix _r;

public:
	/**
	 * \brief Default constructor.
	 */
	QRMatrix() = default;

	/**
	 * \brief		Constructor from matrix.
	 * \param[in] m Matrix for decomposition.
	 */
	QRMatrix(Matrix m);

	/**
	 * \brief	Getter for primal matrix.
	 * \return	primal matrix.
	 */
	Matrix getA() const noexcept;

	/**
	 * \brief		Method for getting answer of equation system.
	 * \param[in] b Vector-column of free members.
	 * \return		Answer for equation system.
	 */
	MyVector solve(MyVector b) const;

	/**
	 * \brief			Printer for Q and R matrix.
	 * \param[in] out	Stream for outputting.
	 */
	void printQR(std::ostream& out) const;

	/**
	 * \brief			Input operator.
	 * \param[in] in	Input stream.
	 * \param[in] m		Matrix for inputting.
	 * \return			Input Stream.
	 */
	friend std::istream& operator >> (std::istream& in, QRMatrix& m);

	/**
	 * \brief			Oputput operator.
	 * \param[in] out	Output stream.
	 * \param[in] m		Matrix for outputting.
	 * \return			Output stream.
	 */
	friend std::ostream& operator << (std::ostream& out, QRMatrix& m);
};

#endif // QR_MATRIX