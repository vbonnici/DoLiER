/*
 * sais.h
 *
 *  Created on: Sep 27, 2013
 *      Author: vbonnici
 */


/*
 * sais.h for sais-lite
 * Copyright (c) 2008-2010 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */



/*
 * adapted by vbonnici to deal with dna5_t sequences
 */

#ifndef SAIS_H_
#define SAIS_H_

//#define SAIS_DEBUG

#include <assert.h>
#include <stdlib.h>

//#include "Sequence.h"
#include "DNA5Alphabet.h"


namespace dolierlib{



class sais_base_array{
public:
	virtual ~sais_base_array(){};
	virtual _si64 get(_si64 i) = 0;
    virtual void set(_si64 i, _si64 val) = 0;
    virtual _si64 update(_si64 i, _si64 val) = 0;
	virtual _si64 decrease(_si64 i, _si64 val) = 0;
};
class sais_int_array: public sais_base_array{
private:
	_si64 *m_a;
	_si64 m_pos;
public:
	sais_int_array(_si64 *m_array, _si64 pos) : m_a(m_array), m_pos(pos) {};
	virtual _si64 get(_si64 i) {
#ifdef SAIS_DEBUG
		if(i<0)std::cout<<"int_array:get_negative "<<i<<"\n";
#endif
		return m_a[m_pos + i];
	};
	virtual void set(_si64 i, _si64 val) {
#ifdef SAIS_DEBUG
		if(i<0)std::cout<<"int_array:set_negative_pos "<<i<<"\n";
		if(val<0)std::cout<<"int_array:set_negative_value "<<val<<"\n";
#endif
		m_a[m_pos + i] =  val;
	};
	 virtual _si64 update(_si64 i, _si64 val){
#ifdef SAIS_DEBUG
		 if(i<0)std::cout<<"int_array:update_negative_pos "<<i<<"\n";
		 if(val<0)std::cout<<"int_array:update_negative_value "<<val<<" ("<<m_a[m_pos +i]<<")\n";
		 if(m_a[m_pos +i]+val<0)std::cout<<"int_array:update_negative_result "<<(m_a[m_pos +i]+val)<<" ("<<m_a[m_pos +i]<<")\n";
#endif
		 return m_a[m_pos +i] += val;
	 };

	 virtual _si64 decrease(_si64 i, _si64 val){
#ifdef SAIS_DEBUG
		 if(i<0)std::cout<<"int_array:decrease_negative_pos "<<i<<"\n";
#endif
		 //return m_a[m_pos +i] -= val;
		 m_a[m_pos +i] -= val;
		 return m_a[m_pos +i];
	 };
};

//class sais_sequence_array: public sais_base_array{
//private:
//	Sequence* m_seq;
//	usize_t m_pos;
//public:
//	sais_sequence_array(Sequence* seq, usize_t pos) : m_seq(seq), m_pos(pos) {};
//	virtual int get(int i) {
//#ifdef SAIS_DEBUG
//		if(i<0)std::cout<<"seq_array:get_negative";
//#endif
//		return m_seq->getCode(m_pos + i);
//	};
//	virtual void set(int i, int val) {
//#ifdef SAIS_DEBUG
//		if(i<0)std::cout<<"seq_array:set_negative_pos";
//		if(val<0)std::cout<<"seq_array:set_negative_value";
//#endif
//		m_seq->setCode(m_pos + i, val);
//	};
//	virtual int update(int i, int val){return 0;};
//};


template<typename T>
class sais_generic_array: public sais_base_array{
private:
	T *m_seq;
	_si64 m_pos;
	_si64 _length;
public:
	sais_generic_array(T *seq, _si64 pos, _si64 length) : m_seq(seq), m_pos(pos), _length(length) {};
	virtual _si64 get(_si64 i) {
#ifdef SAIS_DEBUG
		if(i<0)std::cout<<"seq_array:get_negative";
		if(i+m_pos >= _length)std::cout<<"seq_array:get_overflow "<<i<<" "<<m_pos<<" "<<_length;
#endif
//std::cout<<"seq_array:get "<<i<<" "<<m_pos<<" "<<_length;
		return m_seq[m_pos + i];
	};
	virtual void set(_si64 i, _si64 val) {
#ifdef SAIS_DEBUG
		if(i<0)std::cout<<"seq_array:set_negative_pos";
		if(val<0)std::cout<<"seq_array:set_negative_value";
#endif
		m_seq[m_pos + i] =  val;
	};
	virtual _si64 update(_si64 i, _si64 val){return 0;};
	virtual _si64 decrease(_si64 i, _si64 val){return 0;};
};



void
sais_get_counts(sais_base_array* T, sais_base_array* C, _si64 n, _si64 k) {
	for(_si64 i = 0; i < k; ++i) { C->set(i, 0); }
	for(_si64 i = 0; i < n; ++i) { C->update(T->get(i), 1); }
}

void
sais_get_buckets(sais_base_array* C, sais_base_array* B, _si64 k, bool end) {
	_si64 i, sum = 0;
	if (end != false) { for(i = 0; i < k; ++i) { sum += C->get(i); B->set(i, sum); } }
	else { for(i = 0; i < k; ++i) { sum += C->get(i); B->set(i, sum - C->get(i)); } }
}

void
sais_induce_SA(sais_base_array* T, _si64* SA, sais_base_array* C, sais_base_array* B, _si64 n, _si64 k){
	_si64 b, i, j;
	_si64 c0, c1;
	/* compute SAl */
	if(C == B) { sais_get_counts(T, C, n, k); }
	sais_get_buckets(C, B, k, false); /* find starts of buckets */
	j = n - 1;
	b = B->get(c1 = T->get(j));
	SA[b++] = ((0 < j) && (T->get(j - 1) < c1)) ? ~j : j;
	for(i = 0; i < n; ++i) {
	  j = SA[i]; SA[i] = ~j;
	  if(0 < j) {
		if((c0 = T->get(--j)) != c1) { B->set(c1, b); b = B->get(c1 = c0); }
		SA[b++] = ((0 < j) && (T->get(j - 1) < c1)) ? ~j : j;
	  }
	}
	/* compute SAs */
	if(C == B) { sais_get_counts(T, C, n, k); }
	sais_get_buckets(C, B, k, true); /* find ends of buckets */
	for(i = n - 1, b = B->get(c1 = 0); 0 <= i; --i) {
	  if(0 < (j = SA[i])) {
		if((c0 = T->get(--j)) != c1) { B->set(c1, b); b = B->get(c1 = c0); }
		SA[--b] = ((j == 0) || (T->get(j - 1) > c1)) ? ~j : j;
	  } else {
		SA[i] = ~j;
	  }
	}
}


/* find the suffix array SA of T[0..n-1] in {0..k-1}^n
     use a working space (excluding T and SA) of at most 2n+O(1) for a constant alphabet */
_si64
sais_sa_is(sais_base_array* T, _si64 *SA, _si64 fs, _si64 n, _si64 k){
#ifdef SAIS_DEBUG
	std::cout<<">sais_sa_is\n";
#endif
	sais_base_array *C=NULL, *B=NULL, *RA=NULL;
	_si64 i, j, c, m, p, q, plen, qlen, name, pidx = 0;
	_si64 c0, c1;
	bool diff;

	/* stage 1: reduce the problem by at least 1/2
	   sort all the S-substrings */
	if(k <= fs) {
	  C = new sais_int_array(SA, n);
	  B = (k <= (fs - k)) ? new sais_int_array(SA, n + k) : C;
	} else {
	  B = C = new sais_int_array(new _si64[k], 0);
	}
#ifdef SAIS_DEBUG
	std::cout<<">sais_sa_is:1.0\n";
#endif
	sais_get_counts(T, C, n, k); sais_get_buckets(C, B, k, true); /* find ends of buckets */
	for(i = 0; i < n; ++i) { SA[i] = 0; }
	for(i = n - 2, c = 0, c1 = T->get(n - 1); 0 <= i; --i, c1 = c0) {
	  if((c0 = T->get(i)) < (c1 + c)) { c = 1; }
	  //else if(c != 0) { SA[B->update(c1, -1)] = i + 1; c = 0; }
	  else if(c != 0) { SA[B->decrease(c1, 1)] = i + 1; c = 0; }
	}
#ifdef SAIS_DEBUG
	std::cout<<">sais_sa_is:1.0.1\n";
#endif
	sais_induce_SA(T, SA, C, B, n, k);
#ifdef SAIS_DEBUG
	std::cout<<">sais_sa_is:1.0.2\n";
#endif
	if(C != NULL){
		if(B != C)
			delete B;
		delete C;
	}
	//if(B != NULL) delete B;
	C = NULL; B = NULL;

#ifdef SAIS_DEBUG
	std::cout<<">sais_sa_is:1.1\n";
#endif

	/* compact all the sorted substrings into the first m items of SA
	   2*m must be not larger than n (proveable) */
	for(i = 0, m = 0; i < n; ++i) {
	  p = SA[i];
	  if((0 < p) && (T->get(p - 1) > (c0 = T->get(p)))) {
		for(j = p + 1; (j < n) && (c0 == (c1 = T->get(j))); ++j) { }
		if((j < n) && (c0 < c1)) { SA[m++] = p; }
	  }
	}
	j = m + (n >> 1);
	for(i = m; i < j; ++i) { SA[i] = 0; } /* init the name array buffer */
	/* store the length of all substrings */
	for(i = n - 2, j = n, c = 0, c1 = T->get(n - 1); 0 <= i; --i, c1 = c0) {
	  if((c0 = T->get(i)) < (c1 + c)) { c = 1; }
	  else if(c != 0) { SA[m + ((i + 1) >> 1)] = j - i - 1; j = i + 1; c = 0; }
	}
	/* find the lexicographic names of all substrings */
	for(i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
	  p = SA[i]; plen = SA[m + (p >> 1)]; diff = true;
	  if(plen == qlen) {
		for(j = 0; (j < plen) && (T->get(p + j) == T->get(q + j)); ++j) { }
		if(j == plen) { diff = false; }
	  }
	  if(diff != false) { ++name; q = p; qlen = plen; }
	  SA[m + (p >> 1)] = name;
	}

#ifdef SAIS_DEBUG
	std::cout<<">sais_sa_is:2\n";
#endif
	/* stage 2: solve the reduced problem
	   recurse if names are not yet unique */
	if(name < m) {
	  RA = new sais_int_array(SA, n + fs - m);
	  for(i = m + (n >> 1) - 1, j = n + fs - 1; m <= i; --i) {
		if(SA[i] != 0) { SA[j--] = SA[i] - 1; }
	  }
	  sais_sa_is(RA, SA, fs + n - m * 2, m, name);

	  if(RA != NULL) delete RA;
	  RA = NULL;

	  for(i = n - 2, j = m * 2 - 1, c = 0, c1 = T->get(n - 1); 0 <= i; --i, c1 = c0) {
		if((c0 = T->get(i)) < (c1 + c)) { c = 1; }
		else if(c != 0) { SA[j--] = i + 1; c = 0; } /* get p1 */
	  }
	  for(i = 0; i < m; ++i) { SA[i] = SA[SA[i] + m]; } /* get index */
	}

#ifdef SAIS_DEBUG
	std::cout<<">sais_sa_is:3\n";
#endif
	/* stage 3: induce the result for the original problem */
	if(k <= fs) {
	  C = new sais_int_array(SA, n);
	  B = (k <= (fs - k)) ? new sais_int_array(SA, n + k) : C;
	} else {
	  B = C = new sais_int_array(new _si64[k], 0);
	}
	/* put all left-most S characters into their buckets */
	sais_get_counts(T, C, n, k); sais_get_buckets(C, B, k, true); /* find ends of buckets */
	for(i = m; i < n; ++i) { SA[i] = 0; } /* init SA[m..n-1] */
	for(i = m - 1; 0 <= i; --i) {
	  j = SA[i]; SA[i] = 0;
	  //SA[B->update(T->get(j), -1)] = j;
	  SA[B->decrease(T->get(j), 1)] = j;
	}
	sais_induce_SA(T, SA, C, B, n, k);

	if(C != NULL){
		if(B != C)
			delete B;
		delete C;
	}
	C = NULL; B = NULL;

#ifdef SAIS_DEBUG
	std::cout<<"<sais_sa_is\n";
#endif

	return pidx;
}


//usize_t*
//sais_build_SA(Sequence& seq){
//#ifdef SAIS_DEBUG
//	std::cout<<">sais_build_SA\n";
//#endif
//	usize_t *SA = new usize_t[seq.length()];
//	usize_t n = seq.length();
//	usize_t k = seq.alphabet_size() +1;
//
//	if((seq.length() < n) ||
//	   (k <= 0)) {
//#ifdef SAIS_DEBUG
//	std::cout<<"<sais_build_SA\n";
//#endif
//		return SA;
//	}
//	if(n <= 1) { if(n == 1) { SA[0] = 0; }
//#ifdef SAIS_DEBUG
//	std::cout<<"<sais_build_SA\n";
//#endif
//	return SA;
//	}
//
//	//sa_is(sais_base_array* T, usize_t *SA, int fs, int n, int k){
//	sais_sa_is(new sais_sequence_array(&seq, 0), SA, 0, n, k);
//#ifdef SAIS_DEBUG
//	std::cout<<"<sais_build_SA\n";
//#endif
//	return SA;
//}











usize_t*
sais_build_SA(dna5_t *seq, usize_t length){
#ifdef SAIS_DEBUG
	std::cout<<">sais_build_SA\n";
#endif
	_si64 *SA = new _si64[length];
	_si64 n = length;
	_si64 k = DNA5Alphabet::size() + 1;

/*
	if((length < n) ||
	   (k <= 0)) {
#ifdef SAIS_DEBUG
	std::cout<<"<sais_build_SA\n";
#endif
		return SA;
	}
	if(n <= 1) { if(n == 1) { SA[0] = 0; }
#ifdef SAIS_DEBUG
	std::cout<<"<sais_build_SA\n";
#endif
	return SA;
	}
*/
	//sa_is(sais_base_array* T, usize_t *SA, int fs, int n, int k){

	sais_generic_array<dna5_t> gseq(seq, 0, length);

	sais_sa_is(&gseq, SA, 0, n, k);


#ifdef SAIS_DEBUG
	std::cout<<"<sais_build_SA\n";
#endif

	usize_t *SSA = new usize_t[length];
	for(usize_t i=0; i<length; i++){
		SSA[i]=SA[i];
		if(SA[i]<=0){
			std::cout<<"SA <= 0; "<<i<<" "<<SA[i]<<"\n";
		}
	}
	delete [] SA;

	return SSA;
}


}



#endif /* SAIS_H_ */
