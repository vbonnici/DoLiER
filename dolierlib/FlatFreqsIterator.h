/*
 * FlatFreqsIterator.h
 *
 *  Created on: Nov 12, 2013
 *      Author: vbonnici
 */

#ifndef FLATFREQSITERATOR_H_
#define FLATFREQSITERATOR_H_


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


#include "data_ts.h"

namespace dolierlib{


/*
 * a class to iterate over a 2 columns file in the format  "kmer frequency".
 * Usage: while(IT.next())(kmer = IT.get_kmer(); freq = IT.get_freq(); )
 * goto supposes that kmers are lexicographic ordered
 */


class FlatFreqsIterator{
	std::string file;
	std::ifstream ifs;

	std::string kmer;
	usize_t freq;

public:
	FlatFreqsIterator(std::string &file){
		this->file = file;
		this->ifs.open(file.c_str(), std::ios::in);
		if(!ifs.is_open() || !ifs.good()){
			std::cout<<"Warning can not open file "<<file<<"\n";
		}
		kmer = "";
		freq = 0;
	}
	~FlatFreqsIterator(){
		if(ifs.is_open())
			ifs.close();
	}


	std::string
	get_kmer(){
		return kmer;
	}
	usize_t
	get_freq(){
		return freq;
	}


	bool
	next(){
			if(ifs >> kmer){
			//if(kmer.length() > 0){
				if(ifs >> freq)
					return true;
			}
		return false;
	}


	/*
	 * this supposes that kmers are lexicographic ordered
	 */
	bool
	goto_kmer(std::string s){
		bool n = false;
		if(kmer.length() == 0) n = next();
		else n =true;

		if(n && kmer.compare(s) < 0){
			while((n=next()) && kmer.compare(s) < 0);
		}

		if(n){
			if(kmer.compare(s) == 0){
				return true;
			}
			return false;
		}
		else
			return false;
	}


};


}



#endif /* FLATFREQSITERATOR_H_ */
