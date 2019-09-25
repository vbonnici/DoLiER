/*
 * FlatFreqsIterator.h
 *
 *  Created on: Nov 12, 2013
 *      Author: vbonnici
 */



/*
 * a generic class to iterate over a 1 column file.
 * Usage:  while(IT.next())(T e = IT.get_record)
 */


#ifndef FLATCOLUMNFILEITERATOR_H_
#define FLATCOLUMNFILEITERATOR_H_


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


#include "data_ts.h"

namespace dolierlib{


template<typename T>
class FlatColumnFileIterator{
	std::string file;
	std::ifstream ifs;

	T record;

public:
	FlatColumnFileIterator<T>(std::string &file){
		this->file = file;
		this->ifs.open(file.c_str(), std::ios::in);
		if(!ifs.is_open() || !ifs.good()){
			std::cout<<"Warning can not open file "<<file<<"\n";
		}
	}
	~FlatColumnFileIterator(){
		if(ifs.is_open())
			ifs.close();
	}



	T get_record(){
		return record;
	}


	bool
	next(){
		if(ifs >> record){
				return true;
		}
		return false;
	}
};


}



#endif /* FLATCOLUMNFILEITERATOR_H_ */
