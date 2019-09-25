/*
 * dolier-split-k.cpp
 *
 *  Created on: Jan 30, 2014
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
	std::cout<<"Split a multi-column file according to kmers length. Each output file is a multi-column file containing rows related to a specific kmer length.\n";
	std::cout<<"Usage: "<<scmd<<" <ifile> <nof_cols> <kmers_col_pos>  <opreifx>\n";
	std::cout<<"\t<ifile> input multi-colums file. At least one column must be a kmers list.\n";
	std::cout<<"\t<nof_cols> number of columns in <ifile>.\n";
	std::cout<<"\t<kmers_col_pos> index of the column containing kmers sequences. Indexes start from 0.\n";
	std::cout<<"\t<oprefix> output prefix.\n";
	std::cout<<"For example, considering a 1 column file containing kmers of different lengths, 6,7,8, then 3 different files will be produced named <oprefix>.6 <oprefix>.7 and <oprefix>.8, each one containing only kmers of the same length.\n";
}

class mpars_t : public pars_t{
public:
	std::string ifile;
	int nof_cols;
	int index_col;
	std::string oprefix;


	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		nof_cols = 0;
		index_col = 0;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][ifile]["<<ifile<<"]\n";
		std::cout<<"[ARG][nof_cols]["<<nof_cols<<"]\n";
		std::cout<<"[ARG][index_col]["<<index_col<<"]\n";
		std::cout<<"[ARG][ofile]["<<oprefix<<"]\n";
	}

	virtual void usage(){
		pusage();
	}

	virtual void check(){
		print();
		std::cout<<"check...\n";
		ifile = trim(ifile);
		oprefix = trim(oprefix);

		if(		(ifile.length() == 0)			||
				(oprefix.length() == 0)			||
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
		oprefix = next_string();

		check();
	}

};


class record_t{
public:
	std::string kmer;
	std::string record;
	record_t(std::string _kmer, std::string _record) : kmer(_kmer), record(_record){}
	record_t(const record_t &r){
		this->kmer = r.kmer;
		this->record = r.record;
	}
	friend bool operator== (const record_t& lhs, const record_t& rhs){
		return lhs.kmer == rhs.kmer;
	}

	friend bool operator< (const record_t& lhs, const record_t& rhs) {
		if(lhs.kmer.length() == rhs.kmer.length())
			return lhs.kmer < rhs.kmer;
		return lhs.kmer.length() < rhs.kmer.length();
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
//	ofs.open(pars.oprefix.c_str(), std::ios::out);
//	if(!ofs.is_open() || ofs.bad()){
//		std::cout<<"Error on opening output file : "<<pars.oprefix<<" \n";
//		exit(1);
//	}

	std::vector<record_t> records;

	bool ok = true;
	while(ok){
		std::string cell;
		std::string record = "";
		std::string key;

		for(int i=0; i<pars.nof_cols && ok; i++){
			if(ifs >> cell){
				if(i!=0){
					record += "\t" + cell;
				}
				else{
					record += cell;
				}
				if(i==pars.index_col)
					key = cell;
			}
			else{
				ok = false;
			}
		}

		if(ok){
			records.push_back(record_t(key, record));
		}
	}

	std::sort(records.begin(), records.end());

	size_t plength=0;
//	std::ofstream ofs;
	bool ofs_opened = false;

	for(std::vector<record_t>::iterator IT=records.begin(); IT!=records.end(); IT++){
		if((*IT).kmer.length() != plength){
			plength = (*IT).kmer.length();
			if(ofs_opened){
				ofs.flush();
				ofs.close();
			}
			std::stringstream ss;
			ss<<pars.oprefix<<"."<<plength;
			ofs.open(ss.str().c_str(), std::ios::out);
			if(!ofs.is_open() || ofs.bad()){
				std::cout<<"Error on opening output file : "<<ss.str().c_str()<<" \n";
				exit(1);
			}
			ofs_opened = true;
		}
		ofs<<(*IT).record<<"\n";
	}

	if(ofs_opened == true){
		ofs.flush();
		ofs.close();
	}

//	std::stringstream ss;
//	ss<<pars.prefix<<k<<"."<<m<<pars.suffix;
//	std::cout<<ss.str()<<"\n";


	ifs.close();

	exit(0);
}
