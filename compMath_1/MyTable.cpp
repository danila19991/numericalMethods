#include "MyTable.h"

#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>

void printDelimetr(const std::vector<size_t>& lens, std::ostream& out) noexcept
{
	out << '+';
	for(auto &len:lens)
	{
		for (size_t i = 0;i < len;i++)
			out << '-';
		out << '+';
	}
	out << '\n';
}

template<typename T>
void addLine(const std::vector<size_t>& lens,const std::vector<T>& data, std::ostream& out) noexcept
{
	out << '|';
	for(size_t i=0;i<lens.size();i++)
	{
		out << std::setw(lens[i]) << data[i] << '|';
	}
	out << '\n';
	printDelimetr(lens,out);
}

std::vector<size_t> getLens(const std::vector<std::vector<double>>& data, const size_t pres)
{
	std::vector<size_t> lens(data.begin()->size());

	for (size_t i = 0;i<lens.size();++i)
	{
		double mx = data[0][i];
		for (auto it : data)
			mx = std::max(mx, it[i]);
		size_t k = 0;
		while (abs(mx) >= 10.)
		{
			mx /= 10;
			++k;
		}
		lens[i] = std::max(lens[i], pres + k + 3);
	}

	return lens;
}

void outputInTableWithHeader(const std::vector<std::string>& header,
							 const std::vector<std::vector<double>>& data, const size_t pres,
							 std::ostream& out) noexcept
{
	assert(!data.empty());
	for(auto& it: data)
	{
		assert(it.size() == header.size());
	}

	std::vector<unsigned int> lens = getLens(data,pres);

	for(size_t i=0;i<header.size();++i)
	{
		lens[i] = std::max(lens[i], header[i].size());
	}

	out.precision(pres);
	out << std::fixed;

	printDelimetr(lens, out);

	addLine(lens,header, out);

	for(auto& line:data)
	{
		addLine(lens, line, out);
	}
}

void outputTable(const std::vector<std::vector<double>>& data, const size_t pres, std::ostream& out)
				 noexcept
{
	assert(!data.empty());

	std::vector<unsigned int> lens = getLens(data, pres);

	out.precision(pres);
	out << std::fixed;

	printDelimetr(lens, out);

	for (auto& line : data)
	{
		addLine(lens, line, out);
	}
}
