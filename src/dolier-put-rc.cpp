/*
 * dolier-put-rc.cpp
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
#include <fstream>
#include <algorithm>
#include <limits>

#include "data_ts.h"
#include "trim.h"
#include "pars_t.h"
#include "CstyleIndex.h"

#include <set>
#include <vector>

using namespace dolierlib;

char* scmd;
void pusage(){
	std::cout<<"Output kmers and their reverse complements.\n";
	std::cout<<"Usage: "<<scmd<<" <ifile> <ofile> \n";
	std::cout<<"\t<ifile> input file containing a list of kmers.\n";
	std::cout<<"\t<ofile> output file.\n";
	std::cout<<"For each kmers in <ifile>, it simply output the kmer and its reverse complement. It does not check if <ifile> already contains the rc.";
}

class mpars_t : public pars_t{
public:
	std::string ifile;
	std::string ofile;


	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][ifile]["<<ifile<<"]\n";
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
				(ofile.length() == 0)
		){
			usage();
			exit(1);
		}

	}

	virtual void parse(){
		ifile = next_string();
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


	std::string kmer;
	while(ifs>>kmer){
		ofs<<kmer<<"\n";
		ofs<<reverseComplement(kmer)<<"\n";
	}

	ofs.flush();
	ofs.close();
	ifs.close();
	exit(1);
}
