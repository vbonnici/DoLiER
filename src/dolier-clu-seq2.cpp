/*
 * dolier-clu-seq1.cpp
 *
 *  Created on: Feb 5, 2014
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
#include "seq2.h"
#include "distances.h"



using namespace dolierlib;
using namespace dolierlib::seq1;


char* scmd;
void pusage(){
	std::cout<<"Usage: "<<scmd<<" <input_file> <output_file> <minDist> [--normalize] [--keep-only value] [--weights] [--log-weights] [--yupac] [-v]\n";
}
class mpars_t : public pars_t{
public:
	std::string ifile;
	std::string ofile;
	double minDist;
	bool opt_normalize;
	double keep_only;
	bool opt_weights;
	bool opt_log_weights;
	bool opt_yupac;
	bool opt_v;


	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		minDist = 0;
		opt_normalize = false;
		keep_only = 0;
		opt_weights = false;
		opt_log_weights = false;
		opt_yupac = false;
		opt_v = false;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][ifile]["<<ifile<<"]\n";
		std::cout<<"[ARG][ofile]["<<ofile<<"]\n";
		std::cout<<"[ARG][minDist]["<<minDist<<"]\n";
		std::cout<<"[ARG][normalize]["<<opt_normalize<<"]\n";
		std::cout<<"[ARG][keep_only]["<<keep_only<<"]\n";
		std::cout<<"[ARG][weights]["<<opt_weights<<"]\n";
		std::cout<<"[ARG][log_weights]["<<opt_log_weights<<"]\n";
		std::cout<<"[ARG][yupac]["<<opt_yupac<<"]\n";
	}
	virtual void check(){
		print();
		std::cout<<"check...\n";
		ifile = trim(ifile);
		ofile = trim(ofile);

		if(		(ifile.length() == 0)		||
				(ofile.length() == 0)		||
				(minDist < 0)				||
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
		minDist = next_double();

		std::string cmd;
		while(has_more()){
			cmd = next_string();
			if(cmd == "--normalize"){
				opt_normalize = true;
			}
			else if(cmd == "--keep-only"){
				keep_only = next_double();
			}
			else if(cmd == "--weights"){
				opt_weights = true;
			}
			else if(cmd == "--log-weights"){
				opt_log_weights = true;
			}
			else if(cmd == "--yupac"){
				opt_yupac = true;
			}
			else if(cmd == "-v"){
				opt_v = true;
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
	leading_shuffle(f_sequences, f_lengths, static_cast<size_t>(ceil(f_sequences.size() * 0.1)));
	shuffle(f_sequences, f_lengths, f_sequences.size()*10);

//	std::cout<<"-----------------------------------------------------------------------\n";
//	for(size_t i=0; i<f_sequences.size(); i++){
//		std::cout<<i<<"\t";print_dna5(f_sequences[i], f_lengths[i]);std::cout<<"\n";
//	}
//	std::cout<<"-----------------------------------------------------------------------\n";

	double **vectors;
	size_t vlength = 0;
	double *weights;
	fvector(f_sequences,f_lengths,&vectors, &vlength, &weights);
//	print_matrix(vectors, f_sequences.size(), vlength);

	//double *weights = new double[vlength];
	if(pars.opt_weights){
		//get_fvectors_weights(weights, vlength);
		if(pars.opt_log_weights){
			for(size_t i=0; i<vlength; i++)
				weights[i] = log(weights[i]);
		}
	}
	else{
		for(size_t i=0 ;i<vlength; i++){
			weights[i] = 1;
		}
	}

	std::cout<<"nof sequences "<<f_sequences.size()<<"\n";
	std::cout<<"vector length "<<vlength<<"\n";
	size_t nlength = vlength;

	if(pars.keep_only > 0){
		keep_only(vectors, f_sequences.size(), vlength, 2.0, &nlength, weights);
		std::cout<<"vector length "<<nlength<<"\n";
	}

	if(pars.opt_normalize)
		normalize(vectors, f_sequences.size(), nlength);


	if(pars.opt_v){
		for(size_t i=0; i<f_sequences.size(); i++){
			for(size_t j=0; j<i; j++){
				print_dna5(f_sequences[i], f_lengths[i]);
				std::cout<<" ";
				print_dna5(f_sequences[j], f_lengths[j]);
				std::cout<<" ";
				std::cout<<dist_tanimoto(vectors[i], vectors[j], weights, nlength);
				std::cout<<"\n";
			}
		}
	}


	seq2_algorithm algo(vectors, f_sequences.size(), nlength, weights, pars.minDist);
	//algo.run();
	algo.run(f_sequences, f_lengths);


	if(pars.opt_v){
		std::vector<size_t> dmetoids;
		algo.get_metoids_by_mindist(dmetoids);
		for(size_t i=0; i<algo.clusters.size(); i++){
			std::cout<<"Cluster["<<i<<"]\n";
			std::cout<<"\tdmetoid\t"; print_dna5(f_sequences[dmetoids[i]],f_lengths[dmetoids[i]]); std::cout<<"\n";
			for(std::set<size_t>::iterator IT = algo.clusters[i].begin(); IT!=algo.clusters[i].end(); IT++){
				print_dna5(f_sequences[(*IT)], f_lengths[(*IT)]); std::cout<<"\n";
			}
		}
	}

//	std::vector<size_t> mmetoids;
//	std::vector<size_t> dmetoids;
//	algo.get_metoids_by_mean(mmetoids);
//	algo.get_metoids_by_mindist(dmetoids);
//	for(size_t i=0; i<algo.clusters.size(); i++){
////		std::cout<<"Cluster["<<i<<"]\t";
////		print_dna5(f_sequences[dmetoids[i]],f_lengths[dmetoids[i]]); std::cout<<"\n";
//		std::cout<<"Cluster["<<i<<"]\n";
//		std::cout<<"\tmmetoid\t"; print_dna5(f_sequences[mmetoids[i]],f_lengths[mmetoids[i]]); std::cout<<"\n";
//		std::cout<<"\tdmetoid\t"; print_dna5(f_sequences[dmetoids[i]],f_lengths[dmetoids[i]]); std::cout<<"\n";
//		//std::vector< std::set<size_t> >
//		for(std::set<size_t>::iterator IT = algo.clusters[i].begin(); IT!=algo.clusters[i].end(); IT++){
//			print_dna5(f_sequences[(*IT)], f_lengths[(*IT)]); std::cout<<"\n";
//		}
//	}

	std::ofstream ofs;
	ofs.open((pars.ofile).c_str(), std::ios::out);
	if(!ofs.is_open() || ofs.bad()){
		std::cout<<"Error on opening output file : "<<pars.ofile<<" \n";
		exit(1);
	}

	if(!pars.opt_yupac){
		std::vector<size_t> dmetoids;
		algo.get_metoids_by_mindist(dmetoids);
		for(size_t i=0; i<dmetoids.size(); i++){
			ofs<<to_string(f_sequences[dmetoids[i]],f_lengths[dmetoids[i]])<<"\n";
		}
	}
	else{

	}

	ofs.flush();
	ofs.close();

}

