/*
 * dolier-times.cpp
 *
 *  Created on: Feb 14, 2014
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
#include <map>

using namespace dolierlib;

char* scmd;
void pusage(){
	std::cout<<"Usage: "<<scmd<<" <ifile> <nof_cols> <index_key_col> <index_value_col> <ofile>\n";
}

class mpars_t : public pars_t{
public:
	std::string ifile;
	int nof_cols;
	int key_col;
	int value_col;
	std::string ofile;


	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		nof_cols = 0;
		key_col = 0;
		value_col = 0;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][ifile]["<<ifile<<"]\n";
		std::cout<<"[ARG][nof_cols]["<<nof_cols<<"]\n";
		std::cout<<"[ARG][key_col]["<<key_col<<"]\n";
		std::cout<<"[ARG][value_col]["<<value_col<<"]\n";
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
				(key_col < 0)					||
				(key_col >= nof_cols)			||
				(value_col < 0)					||
				(value_col >= nof_cols)
		){
			usage();
			exit(1);
		}

	}

	virtual void parse(){
		ifile = next_string();
		nof_cols = next_int();
		key_col = next_int();
		value_col = next_int();
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



	std::map<std::string, double> times;
	std::map<std::string, double>::iterator IT;

	bool ok = true;
	while(ok){
		std::string cell;
		std::string key;
		double value;

		for(int i=0; i<pars.nof_cols && ok; i++){
			if(ifs >> cell){
				if(i==pars.key_col)
					key = cell;
				if(i==pars.value_col){
					std::istringstream iss(cell);
					if(!(iss>>value)){
						std::cout<<"wrong column value "<<key<<"\n";
						exit(1);
					}
				}
			}
			else{
				ok = false;
			}
		}

		if(ok){
			IT = times.find(key);
			if(IT == times.end()){
				times.insert(std::pair<std::string, double>(key, value));
			}
			else{
				times.erase(IT);
				times.insert(std::pair<std::string, double>(key, value + IT->second));
			}
		}
	}

	for(IT=times.begin(); IT!=times.end(); IT++){
		ofs<<IT->first<<"\t"<<IT->second<<"\n";
	}

	ofs.flush();
	ofs.close();
	ifs.close();
	exit(0);
}
