/*
 * dolier-clu-seq1-dists.cpp
 *
 *  Created on: Feb 6, 2014
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




double distance0(double *a, double *b, double *weights, size_t dim){
	double dist = 0;

	for(size_t i=0; i<dim; i++){
		dist += ((a[i]-b[i]) * (a[i]-b[i]));
	}
	return dist;
}
double distance1(double *a, double *b, double *weights, size_t dim){
	double dist = 0;

	for(size_t i=0; i<dim; i++){
		dist += ((a[i]-b[i]) * (a[i]-b[i]));
	}
	return sqrt(dist);
}
double distance2(double *a, double *b, double *weights, size_t dim){
	double dist = 0;

	for(size_t i=0; i<dim; i++){
		dist += ((a[i]-b[i]) * (a[i]-b[i]));
	}
	return exp(dist);
}
double distance3(double *a, double *b, double *weights, size_t dim){
	double dist = 0;

	for(size_t i=0; i<dim; i++){
		dist += abs((int)(a[i]-b[i])) / weights[i];
	}
	return dist;
}
double distance4(double *a, double *b, double *weights, size_t dim){
	double dist = 0;

	for(size_t i=0; i<dim; i++){
		dist += ((a[i]-b[i]) * (a[i]-b[i])) / weights[i];
	}
	return sqrt(dist);
}
double distance5(double *a, double *b, double *weights, size_t dim){
	double dist = 0;

	for(size_t i=0; i<dim; i++){
		dist += ((a[i]-b[i]) * (a[i]-b[i])) / weights[i];
	}
	return exp(dist);
}

double distance6(double *a, double *b, double *weights, size_t dim){
	double dist = 0;

	for(size_t i=0; i<dim; i++){
		dist += abs((int)(a[i]-b[i])) * weights[i];
	}
	return dist;
}
double distance7(double *a, double *b, double *weights, size_t dim){
	double dist = 0;

	for(size_t i=0; i<dim; i++){
		dist += ((a[i]-b[i]) * (a[i]-b[i])) * weights[i];
	}
	return sqrt(dist);
}
double distance8(double *a, double *b, double *weights, size_t dim){
	double dist = 0;

	for(size_t i=0; i<dim; i++){
		dist += ((a[i]-b[i]) * (a[i]-b[i])) * weights[i];
	}
	return exp(dist);
}




double tanimoto(double *a, double *b, size_t dim){
	double min = 0;
	double max = 0;
	for(size_t i=0; i<dim; i++){
		if(a[i] < b[i]){
			min += a[i];
			max += b[i];
		}else if(a[i] > b[i]){
			min += b[i];
			max += a[i];
		}
		else{
			min += a[i];
			max += b[i];
		}
	}
	return min/max;
}


//int hamming(std::string &a, std::string &b){
//
//}
int hamming(dna5_t *a, dna5_t *b, size_t a_length, size_t b_length){
	int dist = 0;

	if(a_length == b_length){
		for(size_t i=0; i<a_length; i++){
			if(a[i] != b[i])
				dist++;
		}
	}
	else if(a_length < b_length){
		int mdist = a_length;
		for(size_t i=0; i<b_length - a_length; i++){
			int cdist = 0;
			for(size_t j=0; j<a_length; j++){
				if(a[j] != b[i+j])
					cdist++;
			}
			if(cdist < mdist)
				mdist = cdist;
		}
		dist = mdist;
	}
	else{
		int mdist = b_length;
		for(size_t i=0; i<a_length - b_length; i++){
			int cdist = 0;
			for(size_t j=0; j<b_length; j++){
				if(b[j] != a[i+j])
					cdist++;
			}
			if(cdist < mdist)
				mdist = cdist;
		}
		dist = mdist;
	}


	return dist;
}





size_t tcount(size_t l){
	size_t ret = 0;
	for(size_t i=1; i<=l; i++){
		ret += l -i +1;
	}
	return ret;
}



char* scmd;
void pusage(){
	std::cout<<"Usage: "<<scmd<<" <input_file> [--normalize] [--keep-only value]\n";
}
class mpars_t : public pars_t{
public:
	std::string ifile;
	bool opt_normalize;
	double keep_only;


	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		opt_normalize = false;
		keep_only = 0;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][ifile]["<<ifile<<"]\n";
		std::cout<<"[ARG][normalize]["<<opt_normalize<<"]\n";
		std::cout<<"[ARG][keep_only]["<<keep_only<<"]\n";
	}
	virtual void check(){
		print();
		std::cout<<"check...\n";
		ifile = trim(ifile);

		if(		(ifile.length() == 0)		||
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

	double *weights = new double[vlength];
	get_fvectors_weights(weights, vlength);

	std::cout<<"nof sequences "<<f_sequences.size()<<"\n";
	std::cout<<"vector length "<<vlength<<"\n";
	size_t nlength = vlength;

	if(pars.keep_only > 0){
		keep_only(vectors, f_sequences.size(), vlength, pars.keep_only, &nlength, weights);
		std::cout<<"vector length "<<nlength<<"\n";
	}

	if(pars.opt_normalize)
		normalize(vectors, f_sequences.size(), nlength);




	for(size_t i=0; i<f_sequences.size(); i++){
		for(size_t j=0; j<i; j++){
			print_dna5(f_sequences[i], f_lengths[i]); std::cout<<"\t";
			print_dna5(f_sequences[j], f_lengths[j]); std::cout<<"\t";
			std::cout<<hamming(f_sequences[i], f_sequences[j], f_lengths[i], f_lengths[j])<<"\t";

			std::cout<<tanimoto(vectors[i], vectors[j],nlength)<<"\t";

			std::cout<<distance0(vectors[i], vectors[j],weights, nlength)<<"\t";
			std::cout<<distance1(vectors[i], vectors[j],weights, nlength)<<"\t";
			std::cout<<distance2(vectors[i], vectors[j],weights, nlength)<<"\t";
			std::cout<<distance3(vectors[i], vectors[j],weights, nlength)<<"\t";
			std::cout<<distance4(vectors[i], vectors[j],weights, nlength)<<"\t";
			std::cout<<distance5(vectors[i], vectors[j],weights, nlength)<<"\t";
			std::cout<<distance6(vectors[i], vectors[j],weights, nlength)<<"\t";
			std::cout<<distance7(vectors[i], vectors[j],weights, nlength)<<"\t";
			std::cout<<distance8(vectors[i], vectors[j],weights, nlength)<<"\t";
			std::cout<<"\n";
		}
	}

	std::cout<<tcount(f_lengths[0])<<"\n";
	std::cout<<nlength<<"\n";
//	print_matrix(vectors, f_sequences.size(), nlength);

}
