/*
 * dolier-kmer-list.cpp
 *
 *  Created on: Nov 14, 2013
 *      Author: vbonnici
 */




#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include <vector>

#include "data_ts.h"
#include "timer.h"
#include "trim.h"

#include "DNA5Alphabet.h"
#include "FASTAReader.h"
#include "CstyleIndex.h"
#include "NSAIterator.h"
#include "Cs5DLIndex.h"
#include "main_common.h"
#include "SASearcher.h"
#include "DNA4WordGenerator.h"

#include "DNA4Tree.h"

using namespace dolierlib;


//#define FREQS_DEBUG

char* scmd;
void pusage(){
	std::cout<<"Kmers enumeration.\n";
	std::cout<<"Usage: "<<scmd<<" <o_file> <k> [--no-rc|--push-rc] <LIST_TYPE>\n";
	std::cout<<"\t<o_file> output file. Kmers are listed in lexicographic order.\n";
	std::cout<<"\t<k> kmers length\n";
	std::cout<<"\treverse complement handling:\n";
	std::cout<<"\t\t --push-rc output kmers and their rc, even if they are not present in input sequences\n";
	std::cout<<"\t\t --no-rc do not output reverse complement. Given a kmer A and its reverse complement B, output A if A<=B othervise output B even if B is not in input sequences\n";
	std::cout<<"\t\t none of these, output only present kmers\n";
	std::cout<<"\t<LIST TYPE>:\n";
	std::cout<<"\t\tenum <seq1> <seq2> <...> <seqn>: list only the kmers in input sequences\n";
	std::cout<<"\t\tat_dist <m> <seq1> <seq2> <...> <seqn> : list kmers in input sequences and their neighbors at distance <m>\n";
	std::cout<<"\t<seqX> is a FASTA file containing one or more sequences.\n";

}




class pars_t{
public:
	enum LISTYPE{ENUM, ATDIST};
	enum RCTYPE{ASIS,PUSHRC, NORC};

	std::vector<std::string> iseqs;
	std::string ofile;
	int k;
	int m;
	RCTYPE rctype;
	LISTYPE listype;


	pars_t(){
		k = 0;
		m = 0;
		rctype = ASIS;
		listype = ENUM;
	}
	void print(){
		std::cout<<"[ARG][iseqs]["<<iseqs.size()<<"]";
		for(std::vector<std::string>::iterator IT = iseqs.begin(); IT!=iseqs.end(); IT++)
			std::cout<<"["<<(*IT)<<"]";
		std::cout<<"\n";
		std::cout<<"[ARG][ofile]["<<ofile<<"]\n";
		std::cout<<"[ARG][k]["<<k<<"]\n";
		std::cout<<"[ARG][m]["<<m<<"]\n";

		std::cout<<"[ARG][RC_TYPE][";
		switch(rctype){
			case ASIS:
				std::cout<<"ASIS";
				break;
			case PUSHRC:
				std::cout<<"PUSHRC";
				break;
			case NORC:
				std::cout<<"NORC";
				break;
		}
		std::cout<<"]\n";


		std::cout<<"[ARG][LIST_TYPE][";
		switch(listype){
			case ENUM:
				std::cout<<"ENUM";
				break;
			case ATDIST:
				std::cout<<"ATDIST";
			break;
		}
		std::cout<<"]\n";
	}
	void check(){
		ofile = trim(ofile);
		if(		(iseqs.size() <= 0)		||
				(ofile.length() == 0)	||
				(k < 1)					||
				(m < 0)
		){
			pusage();
			exit(1);
		}

		for(std::vector<std::string>::iterator IT = iseqs.begin(); IT!=iseqs.end(); IT++){
			if(trim(*IT).length() <= 0){
				pusage();
				exit(1);
			}
		}

	}

	void assert_ai(int a_i, int argc, char* cmd){
		if(a_i >= argc){
			pusage();
			exit(1);
		}
	}
	void parse_int(int *n, char* str){
		std::istringstream iss(str);
		if(!(iss>>*n)){
			pusage();
			exit(1);
		}
	}


	void parse(int argc, char* argv[]){
		int a_i = 1;
		std::string cmd;

		assert_ai(a_i,argc, argv[0]);
		ofile = argv[a_i];
		a_i++;


		assert_ai(a_i,argc, argv[0]);
		parse_int(&(k), argv[a_i]);
		a_i++;

		assert_ai(a_i,argc, argv[0]);
		cmd = argv[a_i];
		if(cmd == "--no-rc"){
			rctype = NORC;
			a_i++;
		}
		else if(cmd == "--push-rc"){
			rctype = PUSHRC;
			a_i++;
		}


		assert_ai(a_i,argc, argv[0]);
		cmd = argv[a_i];
		if(cmd == "enum"){
			listype = ENUM;
		}
		else if(cmd == "at-dist"){
			listype = ATDIST;
		}
		else{
			pusage();
			exit(1);
		}
		a_i++;


		assert_ai(a_i,argc, argv[0]);
		if(listype == ATDIST){
			parse_int(&(m), argv[a_i]);
			a_i++;
		}

		assert_ai(a_i,argc, argv[0]);
		while(a_i < argc){
			iseqs.push_back(std::string(argv[a_i]));
			a_i++;
		}
	}

};









void run_enum(pars_t &pars);
void run_atdist(pars_t &pars);

int main(int argc, char* argv[]){
	TIMEHANDLE start_t = start_time();
	scmd = argv[0];

	pars_t pars;
	pars.parse(argc, argv);
	pars.check();
	pars.print();


	switch(pars.listype){
		case pars_t::ENUM:
			run_enum(pars);
			break;
		case pars_t::ATDIST:
			run_atdist(pars);
			break;
		default:
			pusage();
			exit(1);
			break;
	}
}






void write_tree(DNA4Tree &tree, pars_t &pars){
	int  k = pars.k;

	std::ofstream ofs;
	ofs.open((pars.ofile).c_str(), std::ios::out);
	if(!ofs.is_open() || ofs.bad()){
		std::cout<<"Error on opening output file : "<<pars.ofile<<" \n";
		exit(1);
	}


	char *it_syms = new char[k+1];
	it_syms[k] = '\0';

	int si = 0;
	int *c = new int[k];
	typedef DNA4Tree::DNA4TreeNode node_t;
	node_t **stack = new node_t*[k];
	int i;

	usize_t gen_kmers = 0;

	stack[0] = tree.root;
	c[0]=-1;
	while(si >= 0){
		c[si]++;

		if(c[si] < 4){
			if(stack[si]->childs[c[si]] != NULL){
				if(si == k-1){
					gen_kmers++;
					for(i=0; i<k; i++){
						it_syms[i] = (char)DNA5Alphabet::symbolFor(c[i]);
					}
					ofs<<it_syms<<"\n";
				}
				else{
					stack[si+1] = stack[si]->childs[c[si]];
					c[si+1] = -1;
					si++;
				}
			}
		}
		else{
			si--;
		}
	}
	ofs.flush();
	ofs.close();

	std::cout<<"nof written kmers "<<gen_kmers<<"\n";

	delete [] it_syms;
	delete [] c;
	delete [] stack;
}





void run_enum(pars_t &pars){
	DNA4Tree tree;

	int k = pars.k;
	std::string seqfile;
	dna5_t *it_codes = new dna5_t[k];
	dna5_t *it_codes_rev = new dna5_t[k];
	TIMEHANDLE start_t;

	for(std::vector<std::string>::iterator IT = pars.iseqs.begin(); IT!=pars.iseqs.end(); IT++){
		seqfile = *IT;
		std::cout<<seqfile<<"\n";

		std::cout<<"ds...	n";
		start_t = start_time();
		DNA5MS_t f_ms(true);
			load_dna5ms(seqfile, f_ms);
		CsFullyContainer f_ds(true);
			build_fully_ds(f_ms, f_ds);
		NSAIterator it = NSAIterator::begin(f_ms, f_ds, k);
		std::cout<<"done "<<end_time(start_t)<<"\n";



		std::cout<<"tree...\n";
		start_t = start_time();
		if(pars.rctype == pars_t::NORC){
			while(it.next()){
				it.get_kmer(it_codes);

				if(leqRC(it_codes, k)){
					tree.insert(it_codes, k);
				}
				else{
					reverseComplement(it_codes, it_codes_rev, k);
					tree.insert(it_codes_rev, k);
				}
			}
		}
		else if(pars.rctype == pars_t::PUSHRC){
			while(it.next()){
				it.get_kmer(it_codes);
				tree.insert(it_codes, k);
				reverseComplement(it_codes, it_codes_rev, k);
				tree.insert(it_codes_rev, k);
			}
		}
		else{
			while(it.next()){
				it.get_kmer(it_codes);
				tree.insert(it_codes, k);
			}
		}
		std::cout<<"done "<<end_time(start_t)<<"\n";
	}

	start_t = start_time();
	std::cout<<"write...\n";
	write_tree(tree,pars);
	std::cout<<"done "<<end_time(start_t)<<"\n";

	delete [] it_codes;
	delete [] it_codes_rev;

}



void run_atdist(pars_t &pars){
	DNA4Tree tree;

	int k = pars.k;
	int m = pars.m;
	std::string seqfile;
	dna5_t *it_codes = new dna5_t[k];
	dna5_t *it_codes_rev = new dna5_t[k];
	TIMEHANDLE start_t;

	for(std::vector<std::string>::iterator IT = pars.iseqs.begin(); IT!=pars.iseqs.end(); IT++){
		seqfile = *IT;

		std::cout<<seqfile<<"\n";

		start_t = start_time();
		std::cout<<"ds...\n";
		DNA5MS_t f_ms(true);
			load_dna5ms(seqfile, f_ms);
		CsFullyContainer f_ds(true);
			build_fully_ds(f_ms, f_ds);
		NSAIterator it = NSAIterator::begin(f_ms, f_ds, k);
		std::cout<<"done "<<end_time(start_t)<<"\n";

		start_t = start_time();
		std::cout<<"tree...\n";
		if(pars.rctype == pars_t::NORC){
			while(it.next()){
				it.get_kmer(it_codes);

				if(leqRC(it_codes, k)){
					tree.insert(it_codes, k);
				}
				else{
					reverseComplement(it_codes, it_codes_rev, k);
					tree.insert(it_codes_rev, k);
				}


				DNA4WordGenerator wgen(it_codes,k,m);
				while(wgen.next()){
					if(leqRC(wgen.word(),k)){
						tree.insert(wgen.word(),k);
					}
					else{
						reverseComplement(wgen.word(), it_codes_rev, k);
						tree.insert(it_codes_rev, k);
					}
				}
			}
		}
		else if(pars.rctype == pars_t::PUSHRC){
			while(it.next()){
				it.get_kmer(it_codes);

				tree.insert(it_codes, k);
				reverseComplement(it_codes, it_codes_rev, k);
				tree.insert(it_codes_rev, k);

				DNA4WordGenerator wgen(it_codes,k,m);
				while(wgen.next()){
					tree.insert(wgen.word(),k);
					reverseComplement(wgen.word(), it_codes_rev, k);
					tree.insert(it_codes_rev, k);
				}
			}
		}
		else{
			while(it.next()){
				it.get_kmer(it_codes);
				tree.insert(it_codes, k);


				DNA4WordGenerator wgen(it_codes,k,m);
				while(wgen.next()){
					tree.insert(wgen.word(),k);
				}
			}
		}
		std::cout<<"done "<<end_time(start_t)<<"\n";
	}


	start_t = start_time();
	std::cout<<"write...\n";
	write_tree(tree,pars);
	std::cout<<"done "<<end_time(start_t)<<"\n";

	delete [] it_codes;
	delete [] it_codes_rev;
}
