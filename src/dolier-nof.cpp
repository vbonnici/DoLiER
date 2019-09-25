/*
 * dolier-nof.cpp
 *
 *  Created on: Dec 19, 2013
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

using namespace dolierlib;



//#define EFREQS_DEBUG


char* scmd;
void pusage(){
	std::cout<<"Number of sequences, kmers or positions.\n";
	std::cout<<"Usage: "<<scmd<<" <seqs_file> [kmers <k>] [positions <k>]\n";
	std::cout<<"\t <seqs_file> intpu FASTA file containing sequences over the alphabet {A,C,G,T,N}.\n";
	std::cout<<"\toptional [kmers <k>] output number of distinct kmers present in input sequences.\n";
	std::cout<<"\toptional [positions <k>] output number of valid positions for length k. Usually, given a sequence of length N, one can have N-k+1 valid positions where kmers reside, but sequences can caontain Ns.\n";
	std::cout<<"\tif none of [kmers] or [positions], then output the number of sequences contained in <seqs_file>.\n";
	std::cout<<"The number of elements is printed in a row in the format count= <nof>\n";
}


class mpars_t : public pars_t{
public:
	std::string iseqs;
	std::string type;
//	std::string ofile;
	int k;

	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		k = 1;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][iseqs]["<<iseqs<<"]\n";
		std::cout<<"[ARG][type]["<<type<<"]\n";
//		std::cout<<"[ARG][ofile]["<<ofile<<"]\n";
		std::cout<<"[ARG][k]["<<k<<"]\n";
	}
	virtual void check(){
		print();
		std::cout<<"check...\n";
		iseqs = trim(iseqs);
		type = trim(type);

		if(		(iseqs.length() == 0)		||
				(type.length() == 0)		||
				(k < 1)
		){
			usage();
			exit(1);
		}

		if(		type!="kmers" 		&&
				type!="positions" 	&&
				type!="sequences"
		){
			usage();
			exit(1);
		}
	}

	virtual void usage(){
		pusage();
	}

	virtual void parse(){
		iseqs = next_string();
		type = next_string();
		if(type == "kmers" || type=="positions")
			k = next_int();
//		ofile = next_string();
		check();
	}

};

int main(int argc, char* argv[]){
	scmd = argv[0];
	mpars_t pars(argc,argv);
	pars.parse();

	usize_t count = 0;

	if(pars.type == "kmers"){
		DNA5MS_t f_ms(true);
			load_dna5ms(pars.iseqs, f_ms);
		CsFullyContainer f_ds(true);
			build_fully_ds(f_ms, f_ds);
//		SASearcher sas(f_ms.seq, f_ds.SA, f_ms.seq_length);

//		dna5_t *it_codes = new dna5_t[pars.k];
		NSAIterator it = NSAIterator::begin(f_ms, f_ds, pars.k);
		count = nof(it);
	}
	else if(pars.type == "positions"){
		DNA5MS_t f_ms(true);
			load_dna5ms(pars.iseqs, f_ms);
		CsFullyContainer f_ds(true);
			build_fully_ds(f_ms, f_ds);
		SASearcher sas(f_ms.seq, f_ds.SA, f_ms.seq_length);

		//dna5_t *it_codes = new dna5_t[pars.k];
		NSAIterator it = NSAIterator::begin(f_ms, f_ds, pars.k);
		while(it.next()){
			//it.get_kmer(it_codes);
			count += it.occ();
		}
	}
	else if(pars.type == "sequences"){
		DNA5MS_t f_ms(true);
		load_dna5ms(pars.iseqs, f_ms);
		count = f_ms.nof_seqs;
	}
	else{
		pusage();
		exit(1);
	}

	std::cout<<"count= "<<count<<"\n";

	exit(0);
}

