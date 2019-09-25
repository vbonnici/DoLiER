/*
 * dolier-csv-sort.cpp
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

#include "data_ts.h"
#include "trim.h"
#include "pars_t.h"

#include <set>
#include <vector>

using namespace dolierlib;

char* scmd;
void pusage(){
	std::cout<<"Sort an input n columns file according to values of a specific column.\n";
	std::cout<<"Usage: "<<scmd<<" <ifile> <nof_cols> <index_col_num> <index_col_type> <ofile>\n";
	std::cout<<"\t<ifile> input file.\n";
	std::cout<<"\t<nof_cols> number of columns in the input file.\n";
	std::cout<<"\t<index_col_num> index of the column used for sorting. Indexes start form 0.\n";
	std::cout<<"\t<index_col_type> type of values containine din the index column. It can be  string|int|double|usize|size\n";
	std::cout<<"\t<ofile> output file.\n";
	std::cout<<"Example: to sort a 2 column file, containing kmers sequences in the first column and their pvalues in the second column, according to pavlues: ";
	std::cout<<scmd<<" ifile 2 1 double ofile\n";
}

class mpars_t : public pars_t{
public:
	std::string ifile;
	int nof_cols;
	int index_col;
	std::string col_type;
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
		std::cout<<"[ARG][index_col_type]["<<col_type<<"]\n";
		std::cout<<"[ARG][ofile]["<<ofile<<"]\n";
	}

	virtual void usage(){
		pusage();
	}

	virtual void check(){
		print();
		std::cout<<"check...";
		ifile = trim(ifile);
		col_type = trim(col_type);
		ofile = trim(ofile);

		if(		(ifile.length() == 0)			||
				(ofile.length() == 0)			||
				(col_type.length() == 0)		||
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
		col_type = next_string();
		ofile = next_string();

		check();
	}

};


template<typename T>
class record_t{
public:
	T key;
	std::string record;
	record_t(T _key, std::string _record) : key(_key), record(_record){}
	friend bool operator== (const record_t& lhs, const record_t& rhs){
		return lhs.key == rhs.key;
	}

	friend bool operator< (const record_t& lhs, const record_t& rhs) {
		return lhs.key < rhs.key;
	}
};


int main(int argc, char* argv[]){
	scmd = argv[0];
	mpars_t pars(argc, argv);
	pars.parse();

	if(pars.col_type == "string"){}
	else if(pars.col_type == "int"){}
	else if(pars.col_type == "double"){}
	else if(pars.col_type == "usize"){}
	else if(pars.col_type == "size"){}
	else{
		pusage();
		exit(1);
	}


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

	if(pars.col_type == "string"){
		std::vector<record_t<std::string> > records;

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
				records.push_back(record_t<std::string>(key, record));
			}
		}

		std::sort(records.begin(), records.end());
		for(std::vector<record_t<std::string> >::iterator IT=records.begin(); IT!=records.end(); IT++){
			ofs<<(*IT).record<<"\n";
		}
	}
	else if(pars.col_type == "int"){
		std::vector<record_t<int> > records;

		bool ok = true;
		while(ok){
			std::string cell;
			std::string record = "";
			std::string key;
			int kkey;

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
				std::istringstream iss(key);
				if(!(iss>>kkey)){
					std::cout<<"wrong column value "<<key<<"\n";
					exit(1);
				}
				records.push_back(record_t<int>(kkey, record));
			}
		}

		std::sort(records.begin(), records.end());
		for(std::vector<record_t<int> >::iterator IT=records.begin(); IT!=records.end(); IT++){
			ofs<<(*IT).record<<"\n";
		}
	}
	else if(pars.col_type == "double"){
		std::vector<record_t<double> > records;
		bool ok = true;
		while(ok){
			std::string cell;
			std::string record = "";
			std::string key;
			double kkey;

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
				std::istringstream iss(key);
				if(!(iss>>kkey)){
					std::cout<<"wrong column value "<<key<<"\n";
					exit(1);
				}
				records.push_back(record_t<double>(kkey, record));
			}
		}

		std::sort(records.begin(), records.end());
		for(std::vector<record_t<double> >::iterator IT=records.begin(); IT!=records.end(); IT++){
			ofs<<(*IT).record<<"\n";
		}
	}
	else if(pars.col_type == "usize"){
		std::vector<record_t<usize_t> > records;
		bool ok = true;
		while(ok){
			std::string cell;
			std::string record = "";
			std::string key;
			usize_t kkey;

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
				std::istringstream iss(key);
				if(!(iss>>kkey)){
					std::cout<<"wrong column value "<<key<<"\n";
					exit(1);
				}
				records.push_back(record_t<usize_t>(kkey, record));
			}
		}

		std::sort(records.begin(), records.end());
		for(std::vector<record_t<usize_t> >::iterator IT=records.begin(); IT!=records.end(); IT++){
			ofs<<(*IT).record<<"\n";
		}
	}
	else if(pars.col_type == "size"){
		std::vector<record_t<size_t> > records;

		bool ok = true;
		while(ok){
			std::string cell;
			std::string record = "";
			std::string key;
			size_t kkey;

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
				std::istringstream iss(key);
				if(!(iss>>kkey)){
					std::cout<<"wrong column value "<<key<<"\n";
					exit(1);
				}
				records.push_back(record_t<size_t>(kkey, record));
			}
		}

		std::sort(records.begin(), records.end());
		for(std::vector<record_t<size_t> >::iterator IT=records.begin(); IT!=records.end(); IT++){
			ofs<<(*IT).record<<"\n";
		}
	}
	else{
		pusage();
		exit(1);
	}

	ofs.flush();
	ofs.close();
	ifs.close();
	exit(0);
}
