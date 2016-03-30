#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <string>
#include <regex>
#include <numeric>      // std::partial_sum in readNPV_OW function

//************************* 2.READ INPUT DATA ********************************
// This function is simply for reading text file with one value each line
std::vector<double> readSingleValue(std::string filename)
{
	std::ifstream infile(filename);
	std::string buff;
	std::vector<double> valueVec;	//the vector of well control
	while (true){
		if (!getline(infile, buff))
			break;
		valueVec.push_back(std::stod(buff));
	}
	return valueVec;
}