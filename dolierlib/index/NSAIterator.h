/*
 * NSAIterator.h
 *
 *  Created on: Oct 15, 2013
 *      Author: vbonnici
 */

#ifndef NSAITERATOR_H_
#define NSAITERATOR_H_


#include "data_ts.h"
#include "CstyleIndex.h"


//#define DEBUG_NSAITERATOR_H_

namespace dolierlib{

/*
 * a suffix array based class to enumerate kmers (in {A,C,G,T}^k) and their occurrences.
 * Usage:  IT = begin(...); while(IT.next())(kmer = IT.get_kmer(); occ = IT.occ();)
 */

class NSAIterator{
public:
	dna5_t *seq; //indexed sequence
	usize_t length; //sequence length == arrays lengths

	usize_t *SA;
	usize_t *LCP;
	usize_t *NS;

	usize_t first_ncode; //min_i(seq[SA[i]] = N)

	//these are the LCP range [start,end[ of the current kmer
	usize_t i_start;
	usize_t i_end;

	usize_t k;

	NSAIterator(
				dna5_t *_seq,
				usize_t _length,
				usize_t *_SA,
				usize_t *_LCP,
				usize_t *_NS,
				usize_t _first_ncode,
				usize_t i_start,
				usize_t i_end,
				usize_t k)
		: seq(_seq),
		  length(_length),
		  SA(_SA),
		  LCP(_LCP),
		  NS(_NS),
		  first_ncode(_first_ncode){
		this->i_start = i_start;
		this->i_end = i_end;
		this->k = k;
	}

	NSAIterator(const NSAIterator& it) {
		this->seq = it.seq;
		this->length = it.length;
		this->SA = it.SA;
		this->LCP = it.LCP;
		this->NS = it.NS;
		this->first_ncode = it.first_ncode;
		this->i_start = it.i_start;
		this->i_end = it.i_end;
		this->k = it.k;
	}

	~NSAIterator(){}

	bool next(){
#ifdef DEBUG_NSAITERATOR_H_
		std::cout<<"NSAIterator::next()\n";
#endif
		i_start = i_end;

		if(i_start < first_ncode){
			while(i_start < first_ncode && NS[i_start] < k){
				i_start++;
			}

			if(i_start == first_ncode){
				i_end = i_start;
#ifdef DEBUG_NSAITERATOR_H_
		std::cout<<"NSAIterator::next()<false\n";
#endif
				return false;
			}

			i_end = i_start + 1;
			while(i_end < first_ncode
					&& NS[i_end] >= k
					&& LCP[i_end] >= k )
				i_end++;

#ifdef DEBUG_NSAITERATOR_H_
		std::cout<<"NSAIterator::next()<true|"<<i_start<<", "<<i_end<<"\n";
#endif
			return true;

		}
#ifdef DEBUG_NSAITERATOR_H_
		std::cout<<"NSAIterator::next()<false_2\n";
#endif
		return false;
	}

	void get_kmer(dna5_t *kmer){
#ifdef DEBUG_NSAITERATOR_H_
		std::cout<<"NSAIterator::get_kmer()\n";
#endif
		memcpy(kmer, &(seq[SA[i_start]]), k*sizeof(dna5_t));
#ifdef DEBUG_NSAITERATOR_H_
		std::cout<<"NSAIterator::get_kmer()<\n";
#endif
	}

	/*
	 * return number of occurrences of the current kmer (with no mismatches)
	 */
	usize_t occ(){
		return i_end - i_start;
	}



	static const NSAIterator begin(
			dna5_t *_seq,
			usize_t _length,
			usize_t *_SA,
			usize_t *_LCP,
			usize_t *_NS,
			usize_t _first_ncode,
			usize_t k){
		return NSAIterator(_seq, _length, _SA, _LCP, _NS, _first_ncode, 0,0, k);
	}


	static const NSAIterator begin(
			DNA5MS_t& ms,
			CsFullyContainer& fc,
			usize_t k){
		return NSAIterator(ms.seq, ms.seq_length, fc.SA, fc.LCP, fc.NS, fc.first_ncode, 0,0, k);
	}

};


}

#endif /* NSAITERATOR_H_ */
