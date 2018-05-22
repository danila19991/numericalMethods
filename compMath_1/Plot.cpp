#include "Plot.h"
#include <cassert>
#include <fstream>

plotData::plotData(std::string fileName, std::string mode, std::string title, std::string labelX,
	std::string labelY) :_fileName(std::move(fileName)), _mode(std::move(mode)),
	_title(std::move(title)), _lableX(std::move(labelX)), _lableY(std::move(labelY))
{
}

SinglePlot::SinglePlot(std::string fileName, std::string mode, std::string title,
	std::string labelX, std::string labelY):
	plotData(std::move(fileName), std::move(mode), std::move(title), std::move(labelX), 
	std::move(labelY))
{
}

SinglePlot::~SinglePlot()
{
	_storrage.clear();
}


void SinglePlot::outPut()
{
	assert(!_storrage.empty());

	std::ofstream out(_fileName);

	assert(out.is_open());

	out << _storrage.size() << "p|" << _mode << '|' << _title << '|' << _lableX << '|'
		<< _lableY << '\n';

	bool was = false;

	for(const auto& plot: _storrage)
	{
		if (was)
			out << "#\n";
		else
			was = true;

		for(const auto& p: plot)
		{
			out << p.first << ' ' << p.second << '\n';
		}
	}
}

void SinglePlot::addLine(const std::vector<std::pair<double, double>>& line)
{
	_storrage.emplace_back(line);
}
void MultyPlot::outPut()
{
	assert(!_storrage.empty());

	std::ofstream out(_fileName);

	assert(out.is_open());

	out << _storrage.size()*4 << "p|" << _mode << '|' << _title << '|' << _lableX << '|'
		<< _lableY << '\n';

	bool was = false;

	for (const auto& plot : _storrage)
	{
		if (was)
			out << "#\n";
		else
			was = true;

		const size_t sz = plot.begin()->second.size();

		for (const auto& p : plot)
		{
			out << p.first;

			assert(p.second.size() == sz);
			
			for(size_t j=0;j<p.second.size();++j)
			{
				out << ' ' << p.second[j];
			}
			
			out << '\n';
		}
	}
}


void MultyPlot::addLine(const std::vector<std::pair<double, MyVector>>& line)
{
	_storrage.emplace_back(line);
}

MultyPlot::MultyPlot(std::string fileName, std::string mode, std::string title,
	std::string labelX, std::string labelY) :
	plotData(std::move(fileName), std::move(mode), std::move(title), std::move(labelX),
		std::move(labelY))
{
}

MultyPlot::~MultyPlot()
{
	_storrage.clear();
}
