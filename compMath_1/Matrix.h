#ifndef MATRIX
#define MATRIX

#include <vector>

#include "MyTable.h"

/**
 * \brief Class fo matrix.
 */
class Matrix
{
private:

	/**
	 * \brief Data storrage.
	 */
	std::vector<std::vector<double>> _data;

public:
	/**
	 * \brief Default constructor.
	 */
	explicit Matrix() = default;

	/**
	 * \brief			Constructor of square matrix.
	 * \param[in] size	Size of matrix(rows ans columns).
	 */
	explicit Matrix(size_t size);

	/**
	 * \brief				Constructor of rectangle matrix.
	 * \param[in] numRows	Number of rows in matrix.
	 * \param[in] numCols	Number of columns in matrix.
	 */
	explicit Matrix(size_t numRows, size_t numCols);

	/**
	 * \brief			Constructor of matrix.
	 * \param[in] data	Matrixes data.
	 */
	Matrix(std::vector<std::vector<double>> data);

	/**
	 * \brief			Operator for accesing row.
	 * \param[in] pos	Number of row.
	 * \return			Row of matrix.
	 */
	std::vector<double>& operator [] (size_t pos) noexcept;

	/**
	* \brief			Operator for accesing row.
	* \param[in] pos	Number of row.
	* \return			Row of matrix.
	*/
	const std::vector<double>& operator [] (size_t pos) const noexcept;

	/**
	 * \brief	Calculates max of sum of absolute elements of row
	 * \return  Norma inf.
	 */
	double getNorm() const;

	/**
	 * \brief			Getter of column of matrix.
	 * \param[in] pos	Number of column.
	 * \return			Matrix column.
	 */
	std::vector<double> getCol(size_t pos) const noexcept;

	/**
	 * \brief	Getter for number of rows.
	 * \return	Number of rows.
	 */
	size_t numRows() const noexcept;

	/**
	 * \brief	Getter for number of columns.
	 * \return	Number of columns.
	 */
	size_t numCols() const noexcept;

	/**
	 * \brief		Swap rows.
	 * \param[in] l First row for swapping.
	 * \param[in] r Second row for swapping.
	 */
	void swapRows(size_t l, size_t r) noexcept;

	/**
	 * \brief		Swap columns.
	 * \param[in] l First column for swapping.
	 * \param[in] r Second column for swapping.
	 */
	void swapCols(size_t l, size_t r) noexcept;

	/**
	 * \brief			Output matrix in table.
	 * \param[in] pres	Presicion for outputting.
	 * \param[in] out	Stream for outputting.
	 */
	void print(size_t pres = 9,std::ostream& out = std::cout) const noexcept;

	/**
	 * \brief	Make transposed matrix.
	 * \return	Transponed matrix.
	 */
	Matrix trans() const noexcept;

	/**
	 * \brief		Matrix multiplication with new zero columns and zero rows.
	 * \param[in] l Left matrix for multiplication.
	 * \param[in] r Right matrix for miltiplication.
	 * \return		Result matrix.
	 */
	friend  Matrix operator ^ (const Matrix& l, const Matrix& r);

	/**
	* \brief		Matrix multiplication with checking number of columns and number of rows.
	* \param[in] l Left matrix for multiplication.
	* \param[in] r Right matrix for miltiplication.
	* \return		Result matrix.
	*/
	friend  Matrix operator * (const Matrix& l, const Matrix& r);

	/**
	* \brief			Operator for normalising.
	* \param[in] vec	Vector for normalising.
	* \param[in] k		Coefficient.
	* \return			Result vector.
	*/
	friend Matrix operator / (const Matrix& vec, double k);

	/**
	 * \brief			Check if matrix consist of zeroes.
	 * \param[in] eps	Accurancy.
	 * \return			Result.
	 */
	bool isZero(double eps) const;

	/**
	 * \brief			Return E-matrix of dimension dim.
	 * \param[in] dim	Dimension of matrix.
	 * \return			E-matrix.
	 */
	static Matrix onnes(size_t dim);

	/**
	 * \brief				Get random matrix.
	 * \param[in] numRows	Number of rows in matrix.
	 * \param[in] numCols	Number of columns in matrix.
	 * \param[in] minElem	Minimum value of element.
	 * \param[in] maxElem	Maximum value of element.
	 * \return				Random matrix.
	 */
	static Matrix getRandom(size_t numRows, size_t numCols, double minElem, double maxElem);

	/**
	 * \brief				Get random square matrix with diagonal acquisition.
	 * \param[in] size		Number of rows in matrix.
	 * \param[in] minElem	Minimum value of element.
	 * \param[in] maxElem	Maximum value of element.
	 * \return				Random matrix.
	 */
	static Matrix getRandomWithDiagonalAcpuisition(size_t size, double minElem, double maxElem);

	/**
	 * \brief			Inputting operator.
	 * \param[in] in	Inputting stream.
	 * \param[in] m		Matrix for inputting.
	 * \return			Inputting stream.
	 */
	friend std::istream& operator >> (std::istream& in, Matrix& m);

	/**
	 * \brief			Putputting operator.
	 * \param[in] out	Outputting stream.
	 * \param[in] m		Matrix for outputting.
	 * \return			Outputting stream.
	 */
	friend std::ostream& operator << (std::ostream& out, const Matrix& m);
};

#endif //MATRIX