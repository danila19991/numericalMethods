#ifndef MY_VECTOR
#define MY_VECTOR

#include <vector>
#include <iostream>
#include "Matrix.h"

/**
 * \brief Class for adapting vectors with Matrix.
 */
class MyVector
{
private:

	/**
	 * \brief Basic data storrage.
	 */
	std::vector<double> _data;

public:
	/**
	 * \brief Default constructor.
	 */
	MyVector() = default;

	/**
	 * \brief			Constructor of many dimensional vector.
	 * \param[in] size	Dimension of vector.
	 */
	explicit MyVector(size_t size);

	/**
	 * \brief			Constructor of vector.
	 * \param[in] data	new storrage.
	 */
	MyVector(std::vector<double> data);

	/**
	 * \brief			Operator for accesing to any cordinate.
	 * \param[in] pos	Number of cordinates.
	 * \return			Element of vector.
	 */
	double& operator [] (size_t pos) noexcept;

	/**
	* \brief			Operator for accesing to any cordinate.
	* \param[in] pos	Number of cordinates.
	* \return			Element of vector.
	*/
	const double& operator [] (size_t pos) const noexcept;

	/**
	 * \brief		Operator for multipluing matrix on vector.
	 * \param[in] m Matrix for multipluing.
	 * \param[in] v Vector for multipluing.
	 * \return		Result vector.
	 */
	friend MyVector operator *(const Matrix& m, const MyVector& v);

	/**
	 * \brief		Operator for making matrix from vectors.
	 * \param[in] l Left vector.
	 * \param[in] r Righr vector.
	 * \return		Resulting matrix.
	 */
	friend Matrix operator *(const MyVector& l, const MyVector& r);

	/**
	 * \brief		Operator for scalar multiplication.
	 * \param[in] l Left vector.
	 * \param[in] r Righr vector.
	 * \return		Resulting matrix.
	 */
	friend double operator ^(const MyVector& l, const MyVector& r);

	/**
	 * \brief		Operator for subtraction.
	 * \param[in] l Left operand.
	 * \param[in] r Right operand.
	 * \return		Result vector.
	 */
	friend MyVector operator -(const MyVector& l, const MyVector& r);

	/**
	 * \brief		Operator for addiction.
	 * \param[in] l Left operand.
	 * \param[in] r Right operand.
	 * \return		Result vector.
	 */
	friend MyVector operator +(const MyVector& l, const MyVector& r);

	/**
	* \brief			Operator for normalising.
	* \param[in] vec	Vector for normalising.
	* \param[in] k		Coefficient.
	* \return			Result vector.
	*/
	friend MyVector operator / (const MyVector& vec, double k);

	/**
	 * \brief			Operator for normalising.
	 * \param[in] k		Coefficient.
	 * \param[in] vec	Vector for normalising.
	 * \return			Result vector.
	 */
	friend MyVector operator * (double k,const MyVector& vec);

	/**
	 * \brief	Method for get inf-mesuare.
	 * \return	Inf-mesuare.
	 */
	double getNorm() const noexcept;

	/**
	 * \brief	Get size of vector.
	 * \return	Size of vector.
	 */
	size_t size() const noexcept;

	/**
	 * \brief			Check if this vector consist of zeros.
	 * \param[in] eps	Acyracity.
	 * \return			Answer.
	 */
	bool isZero(double eps) const noexcept;

	/**
	 * \brief	Getter of storrage.
	 * \return	Storrage of vector.
	 */
	std::vector<double> getData() const noexcept;

	/**
	 * \brief				Getter for random vector.
	 * \param[in] dim		Size of vector.
	 * \param[in] minElem	Minimum value of element.
	 * \param[in] maxElem	Maximum value of element.
	 * \return				Random vector.
	 */
	static MyVector getRand(size_t dim, double minElem, double maxElem);

	/**
	 * \brief			Operator for outputting vector.
	 * \param[in] out	Stream for outputting.
	 * \param[in] vec	Vector for outputting.
	 * \return			Stream for outputting.
	 */
	friend std::ostream& operator << (std::ostream& out, const MyVector& vec);
	
	/**
	* \brief			Operator for inputting vector.
	* \param[in] in		Stream for inputting.
	* \param[in] vec	Vector for inputting.
	* \return			Stream for inputting.
	*/
	friend std::istream& operator >> (std::istream& in, MyVector& vec);
};


#endif // MY_VECTOR