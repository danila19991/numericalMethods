#ifndef PLOT_H
#define PLOT_H

#include <string>
#include <vector>
#include "MyVector.h"

/**
 * \brief Storrage of data for plot.
 */
class plotData
{
protected:

	/**
	 * \brief Name of file for writing data.
	 */
	std::string _fileName;

	/**
	 * \brief Description of mode for drawing plots.
	 */
	std::string _mode;

	/**
	 * \brief Name for plots.
	 */
	std::string _title;

	/**
	 * \brief Name of horisontal axis.
	 */
	std::string _lableX;

	/**
	 * \brief Name of vertical axis.
	 */
	std::string _lableY;

public:

	/**
	 * \brief				Constructor.
	 * \param[in] fileName	Name of file for writing data.
	 * \param[in] mode		Description of mode for drawing plots.
	 * \param[in] title		Name for plots.
	 * \param[in] labelX	Name of horisontal axis.
	 * \param[in] labelY	Name of vertical axis.
	 */
	plotData(std::string fileName, std::string mode, std::string title, std::string labelX,
		std::string labelY);

	/**
	 * \brief Assignment constructor.
	 */
	plotData(plotData&) = default;

	/**
	 * \brief Assignment operator.
	 */
	plotData& operator = (const plotData&) = default;

	/**
	 * \brief Copy constructor.
	 */
	plotData(plotData&&) = default;

	/**
	 * \brief Copy operator.
	 */
	plotData& operator =(plotData&&) = default;

	/**
	 * \brief Destructor.
	 */
	virtual ~plotData() = default;

	/**
	 * \brief Metod for outputting plots.
	 */
	virtual void outPut() = 0;
};

class SinglePlot : plotData
{
	/**
	 * \brief Storrage of plots.
	 */
	std::vector<std::vector<std::pair<double, double>>> _storrage;

public:
	
	/**
	 * \brief				Constructor.
	 * \param[in] fileName	Name of file for writing data.
	 * \param[in] mode		Description of mode for drawing plots.
	 * \param[in] title		Name for plots.
	 * \param[in] labelX	Name of horisontal axis.
	 * \param[in] labelY	Name of vertical axis.
	 */
	SinglePlot(std::string fileName, std::string mode, std::string title, std::string labelX,
		std::string labelY);

	/**
	 * \brief Assignment constructor.
	 */
	SinglePlot(SinglePlot&) = default;

	/**
	 * \brief Assignment operator.
	 */
	SinglePlot& operator = (const SinglePlot&) = default;

	/**
	 * \brief Copy constructor.
	 */
	SinglePlot(SinglePlot&&) = default;

	/**
	 * \brief Copy operator.
	 */
	SinglePlot& operator =(SinglePlot&&) = default;

	/**
	 * \brief Destructor.
	 */
	~SinglePlot();

	/**
	 * \brief Metod for outputting plots.
	 */
	void outPut() override;

	/**
	 * \brief			Metof for adding new plot.
	 * \param[in] line	New plot.
	 */
	void addLine(const std::vector<std::pair<double, double>>& line);
};

class MultyPlot : public plotData
{
	/**
	 * \brief Storrage of plots.
	 */
	std::vector<std::vector<std::pair<double, MyVector>>> _storrage;

public:

	/**
	 * \brief				Constructor.
	 * \param[in] fileName	Name of file for writing data.
	 * \param[in] mode		Description of mode for drawing plots.
	 * \param[in] title		Name for plots.
	 * \param[in] labelX	Name of horisontal axis.
	 * \param[in] labelY	Name of vertical axis.
	 */
	MultyPlot(std::string fileName, std::string mode, std::string title, std::string labelX,
		std::string labelY);

	/**
	 * \brief Assignment constructor.
	 */
	MultyPlot(MultyPlot&) = default;

	/**
	 * \brief Assignment operator.
	 */
	MultyPlot& operator = (const MultyPlot&) = default;

	/**
	 * \brief Copy constructor.
	 */
	MultyPlot(MultyPlot&&) = default;

	/**
	 * \brief Copy operator.
	 */
	MultyPlot& operator =(MultyPlot&&) = default;

	/**
	 * \brief Destructor.
	 */
	~MultyPlot();

	/**
	 * \brief Metod for outputting plots.
	 */
	void outPut() override;

	/**
	 * \brief			Metof for adding new plot.
	 * \param[in] line	New plot.
	*/
	void addLine(const std::vector<std::pair<double, MyVector>>& line);
};

#endif // PLOT_H