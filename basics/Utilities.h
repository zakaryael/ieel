//
//  Utilities.hpp
//  i-eel
//
//  Created by Jeremie Bec on 25/11/2020.
//

#ifndef Utilities_h
#define Utilities_h

#define _USE_MATH_DEFINES

// loading standard libraries

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <vector>

// loading specific libraries

#include <armadillo>
#include <netcdf.h>

#include "ieel/config.h"


// definition of basic functions

using namespace std;

inline void ErrorMsg(const string& message);

inline string i2s(int n, int length = 0, char pad = '0');

inline std::string r2s(double r);

inline bool fileExists(const std::string& filename);

inline void save(double c, const std::string& filebase);

inline void load(double& c, const std::string& filebase);

inline std::string pathfix(const std::string& path);

inline std::string appendSuffix(const std::string& filename, const std::string& extension);

inline void mkdir(const std::string& dirname);

// source code functions

inline void ErrorMsg(const std::string& message) {
	cerr << message << endl;
	exit(1);
}

inline string i2s(int n, int length, char pad) {
	stringstream ss;
	ss << n;
	string s = ss.str();
	int l = s.length();
	for (int j = l; j < length; ++j) {
		s = pad + s;
	}
	return s;
}

inline std::string r2s(double r) {
	const int Nbuf = 32;
	char buff[Nbuf];
	sprintf(buff, "%g", r);
	return std::string(buff);
}

inline bool fileExists(const std::string& filename) {
	bool res = false;
	struct stat st;
	res = (stat(filename.c_str(), &st) == false) ? true : false;
	return res;
}

inline void save(double c, const std::string& filebase) {
	std::string filename = appendSuffix(filebase, ".asc");
	std::ofstream os(filename.c_str());
	if (!os.good())
		ErrorMsg("save(double, filebase) :  can't open file " + filename);
	os << std::setprecision(17);
	os << c << '\n';  // format can be read by matlab.
}

inline void load(double& c, const std::string& filebase) {
	std::string filename = appendSuffix(filebase, ".asc");
	std::ifstream is(filename.c_str());
	if (!is.good()) {
		std::cerr << "load(double, filebase) :  can't open file " + filename << std::endl;
		exit(1);
	}
	double r = 0;
	is >> r;
	c = r;
}

inline std::string pathfix(const std::string& path) {
	std::string rtn = path;
	if (rtn.length() > 0 && rtn[rtn.length() - 1] != '/')
		rtn += "/";
	return rtn;
}

inline std::string appendSuffix(const std::string& filebase, const std::string& extension) {
	int Lbase = filebase.length();
	int Lext = extension.length();
	std::string filename = filebase;
	if (Lbase < Lext || filebase.substr(Lbase - Lext, Lext) != extension)
		filename += extension;
	return filename;
}

inline void mkdir(const std::string& dirname) {
	::mkdir(dirname.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
}

#endif /* Utilities_h */
