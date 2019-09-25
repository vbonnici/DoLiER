/*
 * main_common.h
 *
 *  Created on: Oct 16, 2013
 *      Author: vbonnici
 */

/*
 * a non-OOP (c style) collection of most common operations used in executables sources (namely main functions)
 */


#ifndef MAIN_COMMON_H_
#define MAIN_COMMON_H_


#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <vector>
#include <limits>

#include "data_ts.h"
#include "timer.h"

#include "DNA5Alphabet.h"
#include "FASTAReader.h"

#include "CstyleIndex.h"
#include "NSAIterator.h"
#include "Cs5DLIndex.h"
//#include "BEMIndex.h"




namespace dolierlib{


/*
 * load all the sequences in the input FASTA file and return a their concatenation ms
 */
inline
void load_dna5ms(const std::string &file, DNA5MS_t &ms){
	std::cout<<"reading sequences...\n";
	TIMEHANDLE start_t = start_time();
	std::vector<dna5_t*> f_sequences;
	std::vector<usize_t> f_lengths;
	FASTAReader f_reader(file);
	f_reader.readAll(f_sequences, f_lengths,true);
	f_reader.close();
	concat(ms, f_sequences, f_lengths, true);
	for(std::vector<dna5_t*>::iterator IT = f_sequences.begin(); IT!=f_sequences.end();IT++)
		delete [] (*IT);
	ms.areUndefTeminated = true;
	std::cout<<end_time(start_t)<<"\n";
}


/*
 * just build the SA of the concatenation sequence and put it in ds
 */
inline
void build_strict_ds(const DNA5MS_t &ms, CsFullyContainer &ds){
	TIMEHANDLE start_t = start_time();
	std::cout<<"building ds...\n";
	ds.SA = build_SA<dna5_t>(ms.seq, ms.seq_length);
	std::cout<<end_time(start_t)<<"\n";
}


/*
 * Fill a CsFullyContainer object starting from a concatenation sequence
 * Buil SA,LCP and NS   and retrieve first_ncode and first_ncode.
 * See CsFullyContainer
 */
inline
void build_fully_ds(const DNA5MS_t &ms, CsFullyContainer &ds){
	TIMEHANDLE start_t = start_time();
	std::cout<<"building ds...\n";
	ds.SA = build_SA<dna5_t>(ms.seq, ms.seq_length);
	ds.LCP = build_LCP<dna5_t>(ms.seq, ms.seq_length, ds.SA);
	ds.NS = build_NS(ms.seq, ms.seq_length, ds.SA);
	ds.first_ncode = where_first_ncode(ds.SA, ds.NS, ms.seq_length);
	ds.first_ecode = where_first_ecode(ms.seq, ds.SA, ms.seq_length);
	std::cout<<end_time(start_t)<<"\n";
}




//template<typename T>
//usize_t nof(T t);
//
//template<typename T>
//usize_t nof<NSAIterator>(NSAIterator it){
//	usize_t count = 0;
//	while(it.next())
//		count++;
//	return count;
//}
//
//template<typename T>
//usize_t nof<BEMIndex::iterator>(BEMIndex::iterator it){
//	usize_t count = 0;
//	while(it.next())
//		count++;
//	return count;
//}



/*
 * Return number of kmers of length k, after the current pointed kmer.
 * If it = begin, than return the total number of distinct kmers.
 */
inline
usize_t nof(NSAIterator it /*it's a copy, so it does not consume the original one*/){
	usize_t count = 0;
	while(it.next())
		count++;
	return count;
}


/*
 * Return the sum of occurences of those kmers, after the current pointed kmer.
 * If it = begin, than return the total number of occurrences.
 */
inline
usize_t nof_occ(NSAIterator it /*it's a copy, so it does not consume the original one*/){
	usize_t count = 0;
	while(it.next())
		count += it.occ();
	return count;
}

//inline
//usize_t nof(BEMIndex::iterator it){
//	usize_t count = 0;
//	while(it.next())
//		count++;
//	return count;
//}


void open_ofs(std::ofstream &ofs, const std::string file){
	ofs.open(file.c_str(), std::ios::out);
	if(!ofs.good()){
		std::cout<<"ERROR on opening: "<<file<<"\n";
		exit(-1);
	}
}
void close_ofs(std::ofstream &ofs){
	ofs.flush();
	ofs.close();
}


}

#endif /* MAIN_COMMON_H_ */
