/*
 * dolier-select.cpp
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

#include <set>
#include <vector>

using namespace dolierlib;

char* scmd;
void pusage(){
	std::cout<<"Select kmers according to their p-value and frequency.\n";
	std::cout<<"Usage: "<<scmd<<" <ifile> <ofile> <pvalue_thr> <nof_seqs_thr> [--fully-output]\n";
	std::cout<<"\t<ifile> input 3 columns file in the format: kmer pvalue seq_count.\n";
	std::cout<<"\t<ofile> output file.\n";
	std::cout<<"\t<pvalue_thr> pvalue threshold, only kmers with a pvalue<=<pvalue_thr> will be selected.\n";
	std::cout<<"\t<nof_seqs_thr> frequence threshold, only kmers that appears in at least <nof_seqs_thr> sequences will be selected.\n";
	std::cout<<"\toptional [--fully-output] output file will be a 3 columns file in the format: kmer pvalue seq_count, otherwise it will be a 1 column file where only kmer are listed.\n";
}

class mpars_t : public pars_t{
public:
	std::string ifile;
	std::string ofile;
	double pvalue_thr;
	int nof_seqs_thr;
	bool opt_fully_output;


	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		pvalue_thr = std::numeric_limits<double>::max();
		nof_seqs_thr = 0;
		opt_fully_output = false;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][ifile]["<<ifile<<"]\n";
		std::cout<<"[ARG][pvalue_thr]["<<pvalue_thr<<"]\n";
		std::cout<<"[ARG][nof_seq_thr]["<<nof_seqs_thr<<"]\n";
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
		pvalue_thr = next_double();
		nof_seqs_thr = next_int();


		std::string cmd;
		while(has_more()){
			cmd = next_string();
			if(cmd == "--fully-output"){
				opt_fully_output = true;
			}
			else{
				usage();
				exit(1);
			}
		}

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
	double pvalue;
	int nof_seqs;

	if(pars.opt_fully_output){
		while((ifs>>kmer) && (ifs>>pvalue) && (ifs>>nof_seqs)){
			if(pvalue <= pars.pvalue_thr && nof_seqs >= pars.nof_seqs_thr){
				ofs<<kmer<<"\t"<<pvalue<<"\t"<<nof_seqs<<"\n";
			}
		}
	}
	else{
		while((ifs>>kmer) && (ifs>>pvalue) && (ifs>>nof_seqs)){
			if(pvalue <= pars.pvalue_thr && nof_seqs >= pars.nof_seqs_thr){
				ofs<<kmer<<"\n";
			}
		}
	}


	ofs.flush();
	ofs.close();
	ifs.close();
	exit(0);
}
