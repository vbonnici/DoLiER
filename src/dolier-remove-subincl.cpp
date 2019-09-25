/*
 * dolier-remove-subincl.cpp
 *
 *  Created on: Jan 24, 2014
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


char* scmd;
void pusage(){
	std::cout<<"Remove sub-included kmers.\n";
	std::cout<<"Usage: "<<scmd<<" <input_file> <output_file> [-r]\n";
	std::cout<<"\t<input_file> is a 1 column file containing kmers of different length.\n";
	std::cout<<"\t<output_file> 1 column output file containing only those kmers that are not a sub-string of a longer kmer.\n";
	std::cout<<"\toptional [-r] also remove all those kmers that are a sub-string of the reverse complement of a longer kmer.\n";
	std::cout<<"Duplicated kmers are considered self sub-included and will not be output.\n";

}
class mpars_t : public pars_t{
public:
	std::string ifile;
	std::string ofile;
	bool opt_r;


	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		opt_r = false;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][ifile]["<<ifile<<"]\n";
		std::cout<<"[ARG][ofile]["<<ofile<<"]\n";
		std::cout<<"[ARG][opt_r]["<<opt_r<<"]\n";
	}
	virtual void check(){
		print();
		std::cout<<"check...\n";
		ifile = trim(ifile);
		ofile = trim(ofile);

		if(		(ifile.length() == 0)		||
				(ofile.length() == 0)
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
			if(cmd == "-r"){
				opt_r = true;
			}
			else{
				usage();
				exit(1);
			}
		}
		check();
	}

};




class Ikmer{
public:
	std::string kmer;
	Ikmer(std::string kmer){
		this->kmer = kmer;
	}

	friend bool operator== (const Ikmer& lhs, const Ikmer& rhs){
		return lhs.kmer == rhs.kmer;
	}

	friend bool operator< (const Ikmer& lhs, const Ikmer& rhs) {
		if(lhs.kmer.length() == rhs.kmer.length())
			return lhs.kmer < rhs.kmer;
		return lhs.kmer.length() < rhs.kmer.length();
	}
};

int main(int argc, char* argv[]){
	scmd = argv[0];

	mpars_t pars(argc,argv);
	pars.parse();


	std::ifstream ifs;
	ifs.open(pars.ifile.c_str(), std::ios::in);
	if(!ifs.is_open() || ifs.bad()){
		std::cout<<"Error on opening input file : "<<pars.ifile<<" \n";
		exit(1);
	}

	std::ofstream ofs;
	ofs.open(pars.ofile.c_str(), std::ios::out);
	if(!ofs.is_open()){
		std::cout<<"Error on opening output file : "<<pars.ofile<<" \n";
		exit(1);
	}


	if(!pars.opt_r){
		std::vector<dna5_t*> f_sequences;
		std::vector<usize_t> f_lengths;
		std::string ikmer;
		while(ifs >> ikmer){
			//ikmer = ikmer+DNA5Alphabet::UNDEF_SYMBOL;
			f_sequences.push_back(to_dna5(ikmer));
			f_lengths.push_back(ikmer.length());
		}
		ifs.close();
		DNA5MS_t f_ms(true);
		//concat(f_ms, f_sequences, f_lengths, true);
		concat(f_ms, f_sequences, f_lengths, false);
	//	for(std::vector<dna5_t*>::iterator IT = f_sequences.begin(); IT!=f_sequences.end();IT++)
	//		delete [] (*IT);
		f_ms.areUndefTeminated = true;

		CsFullyContainer f_ds(true);
			//build_fully_ds(f_ms, f_ds);
			build_strict_ds(f_ms, f_ds);
		Cs5DLIndex f_dlindex(f_ms.seq, f_ms.seq_length, f_ds.SA, f_ms.lengths, f_ms.nof_seqs);
		SASearcher sas(f_ms.seq, f_ds.SA, f_ms.seq_length);
		//print_SA_LCP_N(f_ms.seq, f_ms.seq_length, f_ds.SA, f_ds.LCP, f_ds.NS, 6);


		std::set<Ikmer> tset;

		std::vector<usize_t>::iterator lIT = f_lengths.begin();
		for(std::vector<dna5_t*>::iterator IT = f_sequences.begin(); IT!=f_sequences.end(); ){
			if (sas.count(*IT, (*lIT)-1) == 1){
				tset.insert(Ikmer(to_string(*IT, (*lIT))));
				//ofs<< to_string(*IT, (*lIT)-1);
				//ofs<< "\n";
				//print_dna5(*IT, (*lIT)-1);std::cout<<"\n";
			}
			IT++;
			lIT++;
		}


		for(std::set<Ikmer>::iterator IT=tset.begin(); IT!=tset.end(); IT++){
			ofs<<  (*IT).kmer <<"\n";
		}

		for(std::vector<dna5_t*>::iterator IT = f_sequences.begin(); IT!=f_sequences.end();IT++)
			delete [] (*IT);
	}
	else{
		std::vector<dna5_t*> f_sequences;
		std::vector<usize_t> f_lengths;
		std::string ikmer;
		std::string rikmer;
		while(ifs >> ikmer){
			rikmer = reverseComplement(ikmer);
			f_sequences.push_back(to_dna5(ikmer));
			f_lengths.push_back(ikmer.length());
			f_sequences.push_back(to_dna5(rikmer));
			f_lengths.push_back(ikmer.length());
//			std::cout<<ikmer<<"\n";
//			std::cout<<rikmer<<"\n";
		}
		ifs.close();
		DNA5MS_t f_ms(true);
		concat(f_ms, f_sequences, f_lengths, false);
	//	for(std::vector<dna5_t*>::iterator IT = f_sequences.begin(); IT!=f_sequences.end();IT++)
	//		delete [] (*IT);
		f_ms.areUndefTeminated = true;

		CsFullyContainer f_ds(true);
			build_fully_ds(f_ms, f_ds);
			//build_strict_ds(f_ms, f_ds);
		Cs5DLIndex f_dlindex(f_ms.seq, f_ms.seq_length, f_ds.SA, f_ms.lengths, f_ms.nof_seqs);
		SASearcher sas(f_ms.seq, f_ds.SA, f_ms.seq_length);

//		print_dna5(f_ms.seq,f_ms.seq_length); std::cout<<"\n";
//		print_SA_LCP_N(f_ms.seq, f_ms.seq_length, f_ds.SA, f_ds.LCP, f_ds.NS, 6);

		std::set<Ikmer> tset;


		std::vector<usize_t>::iterator lIT = f_lengths.begin();
		for(std::vector<dna5_t*>::iterator IT = f_sequences.begin(); IT!=f_sequences.end(); ){

//			print_dna5(*IT, *lIT); std::cout<<"\n";

			if(RCeq(*IT, *lIT)){
//				std::cout<<"\tRCeq\n";
				if (sas.count(*IT, (*lIT)) == 2){
					tset.insert(Ikmer(to_string(*IT, (*lIT))));
//					ofs<< to_string(*IT, (*lIT));
//					ofs<< "\n";
				}
			}
			else if (sas.count(*IT, (*lIT)) == 1){
				tset.insert(Ikmer(to_string(*IT, (*lIT))));
//				ofs<< to_string(*IT, (*lIT));
//				ofs<< "\n";
				//print_dna5(*IT, (*lIT)-1);std::cout<<"\n";
			}
			IT++; IT++;
			lIT++; lIT++;
		}

		for(std::set<Ikmer>::iterator IT=tset.begin(); IT!=tset.end(); IT++){
			ofs<<  (*IT).kmer <<"\n";
		}

		for(std::vector<dna5_t*>::iterator IT = f_sequences.begin(); IT!=f_sequences.end();IT++)
			delete [] (*IT);
	}

	ofs.close();

	exit(0);
}
