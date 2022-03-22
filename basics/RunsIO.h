//
//  RunsIO.h
//  i-eel
//
//  Created by Jeremie Bec on 25/11/2020.
//

#ifndef RunsIO_h
#define RunsIO_h

#include "basics/Utilities.h"

using namespace std;

void WriteProcessInfo(int argc, char* argv[], string filename = "processinfo", ios::openmode mode = ios::out);

void WriteProcessInfo(int argc, char* argv[], string filename, ios::openmode mode) {
	ofstream fout(filename.c_str(), mode);
	// Save command-line
	fout << "Command:  ";
	for (int n = 0; n < argc; ++n)
	fout << argv[n] << ' ';
	fout << endl;
	
	// Save current path
	fout << "PWD:      ";
	char currentPath[1024];
	if (getcwd(currentPath, 1023) == NULL)
		ErrorMsg("Error in getcwd())");
	fout << currentPath << endl;
	
	// Save host and pid
	fout << "Host:     ";
	char hostname[1024];
	gethostname(hostname, 1023);
	fout << hostname << ", PID: " << getpid() << endl;
	
	// Save current time
	time_t rawtime;
	struct tm* timeinfo;
	char buffer[80];
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(buffer, 80, "%Y-%m-%d %I:%M:%S", timeinfo);
	fout << "Time:     " << buffer << endl;
	
	// Save code revision and compiler version
	fout << "Version:  " << IEEL_VERSION << endl;
	fout << "Compiler: " << COMPILER_VERSION << endl;
	
}


#endif /* RunsIO_h */
