/*
 * dolier-SWM11.cpp
 *
 *  Created on: Feb 4, 2014
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

#include "SWM11_Point.h"
//#include "SWM11_Comparator.h"
#include "SWM11_SQcomparator.h"
#include "SWM11_Facility.h"
#include "SWM11_frandom.h"

#include "SWM11.h"


using namespace dolierlib;
using namespace dolierlib::SWM11;


char* scmd;
void pusage(){
	std::cout<<"Usage: "<<scmd<<" <input_file> <output_file> <k> [--normalize] [--keep-only value]\n";
}
class mpars_t : public pars_t{
public:
	std::string ifile;
	std::string ofile;
	int k;
	bool opt_normalize;
	double keep_only;


	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		k = 0;
		opt_normalize = false;
		keep_only = 0;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][ifile]["<<ifile<<"]\n";
		std::cout<<"[ARG][ofile]["<<ofile<<"]\n";
		std::cout<<"[ARG][k]["<<k<<"]\n";
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
				(k < 1)						||
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
		k = next_int();

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


	std::cout<<"nof seqs = "<<f_sequences.size()<<"\n";
//	leading_shuffle(f_sequences, f_lengths, f_sequences.size() * 0.1);
//	shuffle(f_sequences, f_lengths, f_sequences.size()*10);
//
//	for(size_t i=0; i<f_sequences.size(); i++){
//		std::cout<<i<<"\t";print_dna5(f_sequences[i], f_lengths[i]);std::cout<<"\n";
//	}

	double **vectors;
	size_t vlength = 0;
	fvector(f_sequences,f_lengths,&vectors, &vlength);
//	print_matrix(vectors, f_sequences.size(), vlength);

	std::cout<<"nof sequences "<<f_sequences.size()<<"\n";
	std::cout<<"vector length "<<vlength<<"\n";
	size_t nlength = vlength;

	if(pars.keep_only > 0){
		keep_only(vectors, f_sequences.size(), vlength, pars.keep_only, &nlength);
		std::cout<<"vector length "<<nlength<<"\n";
	}

	if(pars.opt_normalize)
		normalize(vectors, f_sequences.size(), nlength);


	int maxFacilities =  pars.k;
	int maxSamples = f_sequences.size();
//	maxSamples = 50;

	//SWM11_point_comparator *comp = new SWM11_SQcomparator();
	std::cout<<"_init_\n";
	SWM11_algorithm swm11(vectors, f_sequences.size(), nlength, pars.k, maxFacilities, maxSamples, new SWM11_SQcomparator());
	std::cout<<"_run_\n";
	swm11.run(f_sequences, f_lengths);
	std::cout<<"_print_\n";

	for(int i=0; i<pars.k; i++){
		std::cout<<"Facility["<<i<<"]\n";
		std::cout<<"\t"<<swm11.facilities[i]->p->id<<"\t";
		print_dna5(f_sequences[swm11.facilities[i]->p->id], f_lengths[i]);
		std:cout<<"\n";
		for(std::vector<SWM11_Point*>::iterator IT=swm11.facilities[i]->samples.begin(); IT!=swm11.facilities[i]->samples.end(); IT++){
			std::cout<<"\t"<<(*IT)->id<<"\t";
			print_dna5(f_sequences[(*IT)->id], f_lengths[(*IT)->id]);
			std::cout<<"\n";
		}
	}

}
