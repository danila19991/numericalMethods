#include <iomanip>
#include <fstream>
#include <conio.h>
#include <iostream>

#include "PLUMatrix.h"
#include "QRMatrix.h"
#include "IterationES.h"
#include <cassert>

void first(const std::string& fileName)
{
	Matrix l;

	MyVector vec;

	std::ifstream in(fileName);

	in >> l >> vec;

	PLUMatrix p(l);

	p.printPLU(std::cout);

	std::cout << '\n' << p << '\n';

	std::cout << p.hasSolution(vec) << '\n' << '\n';

	if (p.hasSolution(vec))
	{
		const MyVector ans = p.solve(vec);

		std::cout << ans << '\n' << '\n';

		std::cout << vec - l * ans << '\n';
	}
	else
	{
		std::cout << "There is no solution.\n\n";
	}
}

void second(const std::string& fileName)
{
	Matrix m;

	MyVector vec;

	std::ifstream in(fileName);

	in >> m >> vec;

	QRMatrix qr(m);

	qr.printQR(std::cout);

	std::cout << qr.getA() << "\n\n";

	const MyVector ans = qr.solve(vec);

	std::cout << ans << "\n\n";

	std::cout << vec - (qr.getA()*ans) << "\n\n";
}

void thirdda(size_t size, const double minElem, const double maxElem)
{
	assert(minElem < maxElem);
	
	const Matrix m = Matrix::getRandomWithDiagonalAcpuisition(size, minElem, maxElem);

	const MyVector vec = MyVector::getRand(size,minElem, maxElem);

	std::cout << "A:\n" << m << "\n\n";

	bool isJakobi, isZeidel;

	const MyVector ansj = Jakobi(m, vec, isJakobi, 1e-5);

	const MyVector ansz = Zeidel(m, vec, isZeidel, 1e-5);

	if(isJakobi)
		std::cout << "Jakobi x:\n" << ansj << "\n\n";

	if(isZeidel)
		std::cout << "Zeidel x:\n" << ansz << "\n\n";

	std::cout << "b:\n" << vec << "\n\n";

	if(isJakobi)
		std::cout << "Jakobi b:\n" << m * ansj << "\n\n";

	if(isZeidel)
		std::cout << "Zeidel ans:\n" << m * ansz << "\n\n";
}

void thirdpd(size_t size, const double minElem, const double maxElem)
{
	assert(minElem < maxElem);

	const Matrix m = PLUMatrix::getRandomPositiveDifinite(size, minElem, maxElem);

	const MyVector vec = MyVector::getRand(size, minElem, maxElem);

	std::cout << "A:\n" << m << "\n\n";

	bool isZeidel;

	const MyVector ansz = Zeidel(m, vec, isZeidel, 1e-5);

	if (isZeidel)
		std::cout << "Zeidel x:\n" << ansz << "\n\n";
	if(isZeidel)
		std::cout << "ans - b:\n" << m * ansz -vec << "\n\n";
}
//*/
int main()
{
	std::cout.precision(5);
	std::cout << std::fixed;

	first("input.txt");
	
	//second("input.txt");

	//thirdda(8,-1000,1000);

	//thirdpd(10, -1000, 1000);

	_getch();

	return 0;
}
