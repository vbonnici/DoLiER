/*
 * dolier-merge-columns.cpp
 *
 *  Created on: Jan 28, 2014
 *      Author: vbonnici
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>

#include "FlatFreqsIterator.h"
#include "FlatKmersIterator.h"
#include "FlatColumnFileIterator.h"
#include "FlatKmersValuesIterator.h"

using namespace dolierlib;



char* scmd;
void pusage(){
	std::cout<<"Merge 1 column files in a single n columns file.\n";
	std::cout<<"Usage: "<<scmd<<" <ofile>  <file1> <file2> ... <filen>\n";
	std::cout<<"\t<ofile> output file.\n";
}

int main(int argc, char* argv[]){
	scmd = argv[0];
	if(argc<3){
		pusage();
		exit(1);
	}
	std::string ofile = argv[1];

	std::ofstream ofs;
	ofs.open(ofile.c_str(), std::ios::out);
	if(!ofs.is_open() || ofs.bad()){
		std::cout<<"Error on opening output file : "<<ofile<<" \n";
		exit(1);
	}


	int nof_cols = argc-2;

	FlatColumnFileIterator<std::string> **its = new FlatColumnFileIterator<std::string>*[nof_cols];
	for(int i=0; i<nof_cols; i++){
		std::string ifile = argv[2+i];
		its[i] = new FlatColumnFileIterator<std::string>(ifile);
	}
	std::string *records = new std::string[nof_cols];


	bool ok = true;
	while(ok){
		ok = true;
		for(int i=0; i<nof_cols && ok; i++){
			ok &= its[i]->next();
			if(ok)
				records[i] = its[i]->get_record();
		}
		if(ok){
			for(int i=0; i<nof_cols-1; i++)
				ofs<<records[i]<<"\t";
			ofs<<records[nof_cols-1]<<"\n";
		}

	}


	delete [] records;
	for(int i=0; i<nof_cols; i++){
		delete its[i];
	}

//	for(int i=2)
//
//

	ofs.flush();
	ofs.close();
	exit(0);
}
