/*
 * dolier-get-column.cpp
 *
 *  Created on: Feb 19, 2014
 *      Author: vbonnici
 */




#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "data_ts.h"
#include "trim.h"
#include "pars_t.h"

#include <set>
#include <vector>

using namespace dolierlib;

char* scmd;
void pusage(){
	std::cout<<"Extract a specific column from a multi-columns file.\n";
	std::cout<<"Usage: "<<scmd<<" <ifile> <nof_cols> <index_col_num> <ofile>\n";
	std::cout<<"\t<ifile> input file.\n";
	std::cout<<"\t<nof_cols> number of columns in the input file.\n";
	std::cout<<"\t<index_col_num> index of the column to be extracted.\n";
	std::cout<<"\t<ofile> output file.\n";
}

class mpars_t : public pars_t{
public:
	std::string ifile;
	int nof_cols;
	int index_col;
	std::string ofile;


	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		nof_cols = 0;
		index_col = 0;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][ifile]["<<ifile<<"]\n";
		std::cout<<"[ARG][nof_cols]["<<nof_cols<<"]\n";
		std::cout<<"[ARG][index_col]["<<index_col<<"]\n";
		std::cout<<"[ARG][ofile]["<<ofile<<"]\n";
	}

	virtual void usage(){
		pusage();
	}

	virtual void check(){
		print();
		std::cout<<"check...";
		ifile = trim(ifile);
		ofile = trim(ofile);

		if(		(ifile.length() == 0)			||
				(ofile.length() == 0)			||
				(nof_cols < 1)					||
				(index_col < 0)					||
				(index_col >= nof_cols)
		){
			usage();
			exit(1);
		}

	}

	virtual void parse(){
		ifile = next_string();
		nof_cols = next_int();
		index_col = next_int();
		ofile = next_string();

		check();
	}

};




int main(int argc, char* argv[]){
	scmd = argv[0];
	mpars_t pars(argc, argv);
	pars.parse();



	std::ifstream ifs;
	ifs.open(pars.ifile.c_str(), std::ios::in);
	if(!ifs.is_open() || !ifs.good()){
		std::cout<<"Warning can not open file "<<pars.ifile<<"\n";
		exit(1);
	}

	std::ofstream ofs;
	ofs.open(pars.ofile.c_str(), std::ios::out);
	if(!ofs.is_open() || ofs.bad()){
		std::cout<<"Error on opening output file : "<<pars.ofile<<" \n";
		exit(1);
	}


	bool ok = true;
	while(ok){
		std::string cell;
		std::string key;

		for(int i=0; i<pars.nof_cols && ok; i++){
			if(ifs >> cell){
				if(i==pars.index_col)
					key = cell;
			}
			else{
				ok = false;
			}
		}

		if(ok){
			ofs<<key<<"\n";
		}
	}

	ofs.flush();
	ofs.close();
	ifs.close();
	exit(0);
}
