/*
 * dimers.h
 *
 *  Created on: Dec 17, 2013
 *      Author: vbonnici
 */

/*
 * a collection of functions to work with mono-mers (1-mers) and di-mers (2-mers)
 * No longer supported or documented
 */


#ifndef DIMERS_H_
#define DIMERS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <limits>

#include "timer.h"
#include "data_ts.h"
#include "trim.h"


#include "DNA5Alphabet.h"
#include "FASTAReader.h"
#include "CstyleIndex.h"


//#define DIMERS_DEBUG

namespace dolierlib{
//double *monomers;
//double **dimers;
//load_di_mono_mers(pars.dfile, &monomers, &dimers);

void init_modimers(double **monomers, double ***dimers){
	*monomers = new double[4];
	*dimers = new double*[4];

	for(int i=0; i<4; i++){
		(*monomers)[i] = 0;
		(*dimers)[i] = new double[4];
		for(int j=0; j<4; j++){
			(*dimers)[i][j] = 0;
		}
	}
}
void destroy_modimers(double **monomers, double ***dimers){
	delete [] *monomers;
	for(int i=0; i<4; i++){
		delete [] (*dimers)[i];
	}
}


void uniform_modimers(double *monomers, double **dimers){
	for(int i=0; i<4; i++){
		//monomers[i] = 1.0/4.0;
		monomers[i] = 1;
		for(int j=0; j<4; j++){
			dimers[i][j] = 1.0/8.0;
		}
	}
}


void get_modimers(std::string *files, int nfiles, double *monomers, double **dimers){
	double mono_count = 0;
	double di_count = 0;
	dna5_t *window = new dna5_t[2];

	for(int a_i = 0; a_i<nfiles; a_i++){
		std::string ifile = files[a_i];
		std::cout<<ifile<<"\n";
		FASTAReader reader(ifile);
		std::string name, seq;
		while(reader.nextSeq(&name, &seq)){
			if(seq.length()>0){
				window[1] = DNA5Alphabet::codeFor(seq[0]);
				if(window[1] < 4){
					mono_count++;
					monomers[ window[1] ]++;
				}

				for(size_t i=1; i<seq.length(); i++){
					window[0] = window[1];
					window[1] = DNA5Alphabet::codeFor(seq[i]);
					if(window[1] < 4){
						mono_count++;
						monomers[ window[1] ]++;
						if(window[0] < 4){
							di_count++;
							dimers[ window[0] ][ window[1] ]++;
						}
					}

				}

			}
		}
	}

	for(int i=0; i<4; i++){
		monomers[i] /= mono_count;
		for(int j=0; j<4; j++){
			dimers[i][j] /= di_count;
		}
	}
}

void get_modimers(dna5_t *seq, int seq_length, double *monomers, double **dimers){
	double mono_count = 0;
	double di_count = 0;
	dna5_t *window = new dna5_t[2];

	if(seq_length > 0){
		//window[1] = DNA5Alphabet::codeFor(seq[1]);
		window[1] = seq[0];
		if(window[1] < 4){
			mono_count++;
			monomers[ window[1] ]++;
		}

		for(int i=1; i<seq_length; i++){
			window[0] = window[1];
			//window[1] = DNA5Alphabet::codeFor(seq[i]);
			window[1] = seq[i];
			if(window[1] < 4){
				mono_count++;
				monomers[ window[1] ]++;
				if(window[0] < 4){
					di_count++;
					dimers[ window[0] ][ window[1] ]++;
				}
			}

		}
	}

	for(int i=0; i<4; i++){
		monomers[i] /= mono_count;
		for(int j=0; j<4; j++){
			dimers[i][j] /= di_count;
		}
	}
#ifdef DIMERS_DEBUG
	std::cout<<"modimers("<<seq_length<<",";
	print_dna5(seq, seq_length);
	std::cout<<")\n";
	for(int i=0; i<4; i++){
		std::cout<<DNA5Alphabet::symbolFor(i)<<"\t"<<monomers[i]<<"\n";
	}
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			std::cout<<DNA5Alphabet::symbolFor(i)<<DNA5Alphabet::symbolFor(j)<<"\t"<<dimers[i][j]<<"\n";
		}
	}
#endif
}



void load_modimers(std::string &file, double *monomers, double **dimers){
	std::cout<<"loading real frequencies...\n";

	std::ifstream ifs;
	ifs.open((file).c_str(), std::ios::in);
	if(!ifs.is_open() || ifs.bad()){
		std::cout<<"Error on opening input file : "<<file<<" \n";
		exit(1);
	}

//	*monomers = new double[4];
//	*dimers = new double*[4];
//
//	for(int i=0; i<4; i++){
//		(*monomers)[i] = 0;
//		(*dimers)[i] = new double[4];
//		for(int j=0; j<4; j++){
//			(*dimers)[i][j] = 0;
//		}
//	}

	std::string s;
	double freq;


//	std::cout<<"loading monomers...\n";
	for(int i=0; i<4; i++){
		if(( (ifs>>s) && (ifs>>freq) )){
			dna5_t c = DNA5Alphabet::codeFor(s[0]);
			if(c > 3 ){
				std::cout<<"Error on reading monomers: "<<file<<" \n";
				exit(1);
			}
			monomers[c] = freq;
		}
		else{
			std::cout<<"Error on reading monomers: "<<file<<" \n";
			exit(1);
		}
	}
	std::cout<<"done\n";


//	std::cout<<"loading dimers...\n";
	for(int i=0; i<4*4; i++){
		if(( (ifs>>s) && (ifs>>freq) )){
			dna5_t c1 = DNA5Alphabet::codeFor(s[0]);
			dna5_t c2 = DNA5Alphabet::codeFor(s[1]);

			if(c1>3  || c2>3){
				std::cout<<"Error on reading dimers: "<<file<<" \n";
				exit(1);
			}
			dimers[c1][c2] = freq;
		}
		else{
			std::cout<<"Error on reading dimers: "<<file<<" \n";
			exit(1);
		}
	}
	std::cout<<"done\n";

	ifs.close();
}


double exp_modimer_freq(dna5_t *kmer, int k, double *monomers,double **dimers){
	double x = 1;
//	std::cout<<"x= ";
	for(int j=1; j<k; j++){
//		std::cout<<"*"<<dimers[kmer[j-1]] [kmer[j]]<<" ";
		x *= dimers[kmer[j-1]] [kmer[j]];
	}
	for(int j=1; j<k-1; j++){
//		std::cout<<"/"<<monomers[kmer[j]]<<" ";
		x /= monomers[kmer[j]];
	}
//	std::cout<<"\n";
	return x;
}



void print_modimers(double *monomers, double **dimers){
//	std::cout<<monomers<<"\n";
	for(int i=0; i<4; i++){
		std::cout<<DNA5Alphabet::symbolFor(i)<<"\t"<<monomers[i]<<"\n";
	}
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			std::cout<<DNA5Alphabet::symbolFor(i)<<DNA5Alphabet::symbolFor(j)<<"\t"<<dimers[i][j]<<"\n";
		}
	}
}


}


#endif /* DIMERS_H_ */
