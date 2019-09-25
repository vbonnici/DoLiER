/*
 * dolier-fvector.cpp
 *
 *  Created on: Feb 3, 2014
 *      Author: vbonnici
 */




#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include "data_ts.h"
#include "timer.h"
#include "trim.h"
#include "pars_t.h"
#include "dimers.h"

#include "DNA5Alphabet.h"
#include "FASTAReader.h"
#include "CstyleIndex.h"
#include "NSAIterator.h"
#include "Cs5DLIndex.h"
#include "SASearcher.h"
#include "main_common.h"

#include "fvector.h"

using namespace dolierlib;


char* scmd;
void pusage(){
	std::cout<<"Kmers spectra.\n";
	std::cout<<"Usage: "<<scmd<<" <input_file> <output_file> [--normalize] [--keep-only value]\n";
	std::cout<<"\t<input_file> input file containing a set of kmers.\n";
	std::cout<<"\t<output_file> output file. Each row is the spectrum of the kmer at the same row in <input_file>. Only those columns with at least one element different from 0 will be printed. \n";
	std::cout<<"\toptional [--normalize] normalize kmer spectra by columns.\n";
	std::cout<<"\toptional [--keep-only value] keep only those spectra columns having at least an element greater than value.\n";
}
class mpars_t : public pars_t{
public:
	std::string ifile;
	std::string ofile;
	bool opt_normalize;
	double keep_only;


	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		opt_normalize = false;
		keep_only = 0;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][ifile]["<<ifile<<"]\n";
		std::cout<<"[ARG][ofile]["<<ofile<<"]\n";
		std::cout<<"[ARG][normalize]["<<opt_normalize<<"]\n";
		std::cout<<"[ARG][keep_only]["<<keep_only<<"]\n";
	}
	virtual void check(){
		print();
		std::cout<<"check...\n";
		ifile = trim(ifile);
		ofile = trim(ofile);

		if(		(ifile.length() == 0)		||
				(ofile.length() == 0)		||
				(keep_only < 0)
		){
			usage();
			exit(1);
		}
	}

	virtual void usage(){
		pusage();
	}

	virtual void parse(){
		ifile = next_string();
		ofile = next_string();

		std::string cmd;
		while(has_more()){
			cmd = next_string();
			if(cmd == "--normalize"){
				opt_normalize = true;
			}
			else if(cmd == "--keep-only"){
				std::cout<<"keep only!\n";
				keep_only = next_double();
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

	mpars_t pars(argc,argv);
	pars.parse();


	/*
	 * - read kmers form file
	 * - get spectra vectors
	 * - output them on file
	 */


	size_t max_length = 0;

	std::ifstream ifs;
	ifs.open(pars.ifile.c_str(), std::ios::in);
	if(!ifs.is_open() || ifs.bad()){
		std::cout<<"Error on opening input file : "<<pars.ifile<<" \n";
		exit(1);
	}
	std::vector<dna5_t*> f_sequences;
	std::vector<usize_t> f_lengths;
	std::string ikmer;
	while(ifs >> ikmer){
		f_sequences.push_back(to_dna5(ikmer));
		f_lengths.push_back(ikmer.length());
		if(ikmer.length() > max_length)
			max_length = ikmer.length();
	}
	ifs.close();


//	std::cout<<"nof seqs = "<<f_sequences.size()<<"\n";
//	leading_shuffle(f_sequences, f_lengths, f_sequences.size() * 0.1);
//	shuffle(f_sequences, f_lengths, f_sequences.size()*10);
//
//	for(size_t i=0; i<f_sequences.size(); i++){
//		std::cout<<i<<"\t";print_dna5(f_sequences[i], f_lengths[i]);std::cout<<"\n";
//	}

	double **vectors;
	size_t vlength = 0;
	double *weights;
	fvector(f_sequences,f_lengths,&vectors, &vlength, &weights);
//	print_matrix(vectors, f_sequences.size(), vlength);

	std::cout<<"nof sequences "<<f_sequences.size()<<"\n";
	std::cout<<"vector length "<<vlength<<"\n";
	size_t nlength = vlength;

	if(pars.keep_only > 0){
		keep_only(vectors, f_sequences.size(), vlength, 2.0, &nlength);
		std::cout<<"vector length "<<nlength<<"\n";
//		print_matrix(vectors, f_sequences.size(), nlength);
	}

	if(pars.opt_normalize){
		normalize(vectors, f_sequences.size(), nlength);
	}


	std::ofstream ofs;
	//ofs.open(pars.ofile.c_str(), std::fstream::out | std::fstream::app);
	ofs.open(pars.ofile.c_str(), std::fstream::out);
	if(!ofs.is_open() || ofs.bad()){
		std::cout<<"Error on opening output file : "<<pars.ofile<<" \n";
		exit(1);
	}

	for(size_t i=0; i<f_sequences.size(); i++){
		ofs<<i<<",";
		for(size_t j=0; j<nlength-1; j++){
			ofs<<vectors[i][j]<<",";

		}
		ofs<<vectors[i][nlength-1];
		ofs<<"\n";
	}

	ofs.flush();
	ofs.close();

	exit(0);

}
