
#ifndef FLATKMERSVALUESITERATOR_H_
#define FLATKMERSVALUESITERATOR_H_


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


#include "data_ts.h"

#include "FlatKmersIterator.h"
#include "FlatColumnFileIterator.h"

namespace dolierlib{

/*
 * a generic class to iterate over a 2 columns file in the format  "kmer value".
 * Usage: while(IT.next())(kmer = IT.get_kmer(); T value = IT.get_value(); )
 * goto supposes that kmers are lexicographic ordered
 */

template<typename T>
class FlatKmersValuesIterator{
	FlatKmersIterator &kit;
	FlatColumnFileIterator<T> &cit;

	std::string kmer;
	T value;

public:
	FlatKmersValuesIterator<T>(FlatKmersIterator &_kit, FlatColumnFileIterator<T> &_cit) : kit(_kit), cit(_cit){
		kmer = "";
	}
	~FlatKmersValuesIterator(){
	}


	std::string
	get_kmer(){
		return kmer;
	}
	usize_t
	get_value(){
		return value;
	}


	bool
	next(){

		bool ret = (kit.next() && cit.next());
		if(ret){
			kmer = kit.get_kmer();
			value= cit.get_record();
		}
		return ret;
	}


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



#endif /* FLATKMERSVALUESITERATOR_H_ */
