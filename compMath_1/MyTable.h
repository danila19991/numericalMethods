#ifndef MY_TABLE
#define MY_TABLE

#include <vector>
#include <string>
#include <iostream>

/**
 * \brief			 Function for outputting table with data and header.
 * \param[in] header Header of the table.
 * \param[in] data	 Result of calcilations.
 * \param[in] pres	 Presicion for outputing.
 * \param[in] out	 Stream for outputting.
 */
void outputInTableWithHeader(const std::vector<std::string>& header,
							 const std::vector<std::vector<double>>& data, size_t pres = 9,
							 std::ostream& out = std::cout);

/**
 * \brief			Function for outputting table with data.
 * \param[in] data	Data for outputting.
 * \param[in] pres	Precision for outputting.
 * \param[in] out	Stream for outputting.
 */
void outputTable(const std::vector<std::vector<double>>& data, size_t pres = 9,
				 std::ostream& out = std::cout);

#endif // MY_TABLE
