/*
 * dolier-pwm-scores.cpp
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
#include "CstyleIndex.h"
#include "pwm_supp.h"

#include <set>
#include <vector>

using namespace dolierlib;

char* scmd;
void pusage(){
	std::cout<<"Get best, average and standard deviation of kmers scores over a PWM.\n";
	std::cout<<"Usage: "<<scmd<<" <ifile> <ofile> <pwm_file> <eps> [-r]\n";
	std::cout<<"\t<ifile> input file containing a list of kmers.\n";
	std::cout<<"\t<ofile> best_score avg_score and sd  will be appended at the end of <ofile>.\n";
	std::cout<<"\t<pwm_file> file containing a PWM. Unlikely the standard format, expected PWM is a matrix having 4 columns (for A,C,G,T) and k rows (one for each position). A initial description line is also expected.\n";
	std::cout<<"\t<eps> epsilon used for pseudocounts normalization. If <eps>==0 then no normalization is applied.\n";
	std::cout<<"\toptional [-r] consider reverse complement, too.\n";
}

class mpars_t : public pars_t{
public:
	std::string ifile;
	std::string ofile;
	std::string pwmfile;
	double eps;
	bool opt_r;


	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		eps = 0;
		opt_r = false;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][ifile]["<<ifile<<"]\n";
		std::cout<<"[ARG][ofile]["<<ofile<<"]\n";
		std::cout<<"[ARG][pwmfile]["<<pwmfile<<"]\n";
		std::cout<<"[ARG][eps]["<<eps<<"]\n";
		std::cout<<"[ARG][opt_r]["<<opt_r<<"]\n";
	}

	virtual void usage(){
		pusage();
	}

	virtual void check(){
		print();
		std::cout<<"check...";
		ifile = trim(ifile);
		ofile = trim(ofile);
		pwmfile = trim(pwmfile);

		if(		(ifile.length() == 0)			||
				(ofile.length() == 0)			||
				(pwmfile.length() == 0)
		){
			usage();
			exit(1);
		}

	}

	virtual void parse(){
		ifile = next_string();
		ofile = next_string();
		pwmfile = next_string();
		eps = next_double();
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

	int k_pwm= 0;
//	double **pwm = get_gtpwm(pars.pwmfile, &k_pwm);
//	normalize(pwm, k_pwm);
	double **e_pwm = get_gtpwm(pars.pwmfile, &k_pwm);
	if(pars.eps != 0){
		pseudocount(e_pwm, k_pwm, pars.eps);
		//normalize(e_pwm, k_pwm);
	}
	else{
		normalize(e_pwm, k_pwm);
	}

	std::ofstream ofs;
	ofs.open(pars.ofile.c_str(), std::fstream::out | std::fstream::app);
	if(!ofs.is_open() || ofs.bad()){
		std::cout<<"Error on opening output file : "<<pars.ofile<<" \n";
		exit(1);
	}

	double best_score = 0;
	double avg_score = 0;
	double sd_score = 0;

	//double t_best_score = best_prob_score(pwm,k_pwm, )
//	std::map<int,double> t_best_scores;
	int pk = 0;
	double t_best_score = 0;
	double score;
	int k;
	std::string kmer;
	std::string rckmer;
	usize_t count = 0;
	while(ifs>>kmer){
		k = kmer.length();
		if(t_best_score == 0 || k!=pk){
			t_best_score = best_prob_score(e_pwm,k_pwm, k);
		}
		pk = k;
		score = prob_score(to_dna5(kmer), k,e_pwm, k_pwm) /t_best_score;
		avg_score += score;
		count++;
		if(score > best_score)
			best_score = score;

		if(pars.opt_r){
			rckmer = reverseComplement(kmer);
			score = prob_score(to_dna5(rckmer), k,e_pwm, k_pwm) /t_best_score;
			avg_score += score;
			count++;
			if(score > best_score)
				best_score = score;
		}
	}
	avg_score /= static_cast<double>(count);

	ifs.close();
	ifs.open(pars.ifile.c_str(), std::ios::in);
	if(!ifs.is_open() || !ifs.good()){
		std::cout<<"Warning can not open file "<<pars.ifile<<"\n";
		exit(1);
	}

	pk = 0;
	while(ifs>>kmer){
		k = kmer.length();
		if(t_best_score == 0 || k!=pk){
			t_best_score = best_prob_score(e_pwm,k_pwm, k);
		}
		pk = k;
		score = prob_score(to_dna5(kmer), k,e_pwm, k_pwm) /t_best_score;
		sd_score = (score - avg_score) * (score - avg_score);

		if(pars.opt_r){
			rckmer = reverseComplement(kmer);
			score = prob_score(to_dna5(rckmer), k,e_pwm, k_pwm) /t_best_score;
			sd_score = (score - avg_score) * (score - avg_score);
		}
	}
	sd_score /= static_cast<double>(count);
	sd_score = sqrt(sd_score);

	ofs<<best_score<<"\t"<<avg_score<<"\t"<<sd_score<<"\t";

	ofs.flush();
	ofs.close();
	ifs.close();
	exit(0);
}
