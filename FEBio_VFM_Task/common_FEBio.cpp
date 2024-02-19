#include "common_FEBio.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

void write_to_log_2(FEModel* fem, const std::string& logs_string, std::ofstream& outFile)
{
	write_log(fem, 0, logs_string.c_str());
	outFile << logs_string;
}