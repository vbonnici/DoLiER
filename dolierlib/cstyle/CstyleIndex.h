/*
 * CstyleIndex.h
 *
 *  Created on: Oct 8, 2013
 *      Author: vbonnici
 */

#ifndef CSTYLEINDEX_H_
#define CSTYLEINDEX_H_

#include <stdio.h>
#include <string.h>
#include <vector>
#include <bitset>

#include "data_ts.h"
#include "sais.h"

#include "DNA5Alphabet.h"

namespace dolierlib{


/*
 * print the suffix (of length n) in the suffix array at position pos
 */
void
print_dna5_SA(dna5_t *seq/*the indexed sequence*/, usize_t *SA, usize_t pos, int n){
	for(int i=0; i<n; i++){
		std::cout<<DNA5Alphabet::symbolFor(seq[SA[pos]+i]);
	}
}




//inline
//void
//fill_code4_from_SA(dna5_t *seq, usize_t *SA, usize_t pos, _ui16 *c4){
////	print_dna5_SA(seq,SA,pos,8);std::cout<<"\n";
//	int mask = 0x3;
//	int shift = 14;
//	int cc5;
//	*c4 = 0;
//	for(int i=0; i< 8; i++ ){
//		cc5 = seq[SA[pos]+i];
//		if(cc5 > DNA5Alphabet::T) cc5 = DNA5Alphabet::A;
////		std::cout<<shift<<"|"<<cc5<<"|"<<((cc5&mask)<<shift)<<"\n";
////		std::cout<<std::bitset<16>((cc5&mask)<<shift)<<"\n";
//		*c4 |= (cc5&mask)<<shift;
//		shift -= 2;
//	}
////	std::cout<<(*c4)<<"\n";
////	std::cout<<std::bitset<16>(*c4)<<"\n";
//}


/*
 * convert a nucleotide sequence in a dna5_t array
 */
dna5_t*
to_dna5(std::string &s){
	dna5_t *s5 = new dna5_t[s.length()];
	for(usize_t i=0; i<s.length(); i++){
		s5[i] = DNA5Alphabet::codeFor(s[i]);
	}
	return s5;
}

/*
 * convert a dna5_t array of length k in a string
 */
std::string
to_string(dna5_t *t, int k){
	std::stringstream ss;
	std::string s;
	for(int i=0;i<k;i++){
		s = s + DNA5Alphabet::symbolFor(t[i]);
	}
	return s;
}

void
print_dna5(dna5_t *s, int n){
	for(int i=0;i<n;i++){
		std::cout<<DNA5Alphabet::symbolFor(s[i]);
	}
}

void
print_dna5(const dna5_t *s, int n){
	for(int i=0;i<n;i++){
		std::cout<<DNA5Alphabet::symbolFor(s[i]);
	}
}


void
print_dna5(dna5_t *s, int p, int n){
	for(int i=p;i<p+n;i++){
		std::cout<<DNA5Alphabet::symbolFor(s[i]);
	}
}


/*
 * push the reverse complement of sequence s (of length n) in d
 */
inline
void
reverseComplement(dna5_t *s, dna5_t *d, int n){
	for(int i=0; i<n; i++){
		d[n-i-1] = DNA5Alphabet::complement(s[i]);
	}
}


/*
 * tell if s (of length k) is less/equal to its reverse complement
 */
inline
bool
leqRC(dna5_t *s, int k){
	for(int i=0; i<k; i++){
		if(DNA5Alphabet::complement(s[i]) != s[k-i-1])
			return (DNA5Alphabet::complement(s[i]) >= s[k-i-1]);
	}
	return true;
}


/*
 * tell if s (of length k) is equal to its reverse complement
 */
inline
bool
RCeq(dna5_t *s, int k){
////	if(k%2==0){
//	//std::cout<<"#"<<k<<"|"<<(k/2)<<"\n";
//		//for(int i=0; i<k/2; i++){
//	for(int i=0; i<k; i++){
//			if(DNA5Alphabet::complement(s[i]) != s[k-i-1])
//				return false;
//		}
////	}
////	else{
////		for(int i=0; i<k; i++){
////			if(DNA5Alphabet::complement(s[i]) != s[k-i-1])
////				return false;
////		}
////	}
//	return true;

	if(k%2==0){
		for(int i=0; i<k/2; i++){
			if(DNA5Alphabet::complement(s[i]) != s[k-i-1])
				return false;
		}
		return true;
	}
	return false;
}



/*
 * return the reverse complement of the nucleotide sequence s
 */
inline
std::string
reverseComplement(std::string &s){
	std::string r = "";
	//for(int i=0; i<s.length(); i++){
	for(size_t i=0; i<s.length(); i++){
//		std::cout<<"r:"<<(DNA5Alphabet::symbolFor( DNA5Alphabet::complement( DNA5Alphabet::codeFor(s[i]) ) ))<<"\n";
		r = DNA5Alphabet::symbolFor( DNA5Alphabet::complement( DNA5Alphabet::codeFor(s[i]) ) ) + r;
	}
	return r;
}


/*
 * tell if two nucleotide sequences of length n are equal
 */
inline
bool
equal_dna5(dna5_t *s, dna5_t *d, int n){
	for(int i=0; i<n; i++)
		if(s[i]!=d[i])	return false;
	return true;
}


/*
 * return
 *   0 if a=b
 *  <0 if a<b
 *  >0 if a>b
 *  in the lexicographic order
 */
inline
int
compare(dna5_t *a, dna5_t *b, int k){
	for(int i=0; i<k; i++){
		if(a[i] != b[i]){
			//it works because in DNA5Alphabet  A=0,C=1,G=2,T=3,N=4,$=5
			return a[i] - b[i];
		}
	}
	return 0;
}




/*
 * An object to store the concatenation of a set of sequences plus other data.
 */
class  DNA5MS_t{
public:
	dna5_t *seq; //the concatenation sequence
	usize_t seq_length; //the length of the concatenation sequence
	usize_t nof_seqs; //the number of original sequences
	usize_t *lengths; //the length of each original sequence, in the same order in which they are concatenated
	bool areUndefTeminated; //true if the concatenation is made by putting $ between sequences

	bool freeOnExit; //if true then delete array on deconstruction

	DNA5MS_t(){
		seq = NULL;
		seq_length = 0;
		nof_seqs = 0;
		lengths = NULL;
		freeOnExit = false;
		areUndefTeminated = false;
	}
	DNA5MS_t(bool _freeOnExit){
		seq = NULL;
		seq_length = 0;
		nof_seqs = 0;
		lengths = NULL;
		areUndefTeminated = false;
		freeOnExit = _freeOnExit;
	}

//	DNA5MS_t(const DNA5MS_t &ms){
//		seq = ms.seq;
//		seq_length = ms.seq_length;
//		nof_seqs = ms.nof_seqs;
//		lengths = ms.lengths;
//	}
	~DNA5MS_t(){
		if(freeOnExit){
			if(seq!=NULL)	delete [] seq;
			if(lengths!= NULL)	delete [] lengths;
		}
	}
};



/*
 * concatenate a set of sequences, putting $ between them and at the end of the result sequence
 */
void
concat(DNA5MS_t& gmseq, //output concatenation
		std::vector<dna5_t*>& a_sequences, //input sequences
		std::vector<usize_t>& a_lengths, //input sequences lengths
		bool areUndefTerminated //true if input sequences are already $ terminated
	){


	gmseq.nof_seqs = a_sequences.size();

	gmseq.lengths = new usize_t[gmseq.nof_seqs];

	if(areUndefTerminated){
		gmseq.seq_length = 0;
		for(usize_t i=0; i<gmseq.nof_seqs; i++){
			gmseq.lengths[i] = a_lengths.at(i);
			gmseq.seq_length += gmseq.lengths[i];
		}

		gmseq.seq = new dna5_t[gmseq.seq_length];

		dna5_t *ss;
		int i=0, ii=0, l=0;
		for(std::vector<dna5_t*>::iterator IT = a_sequences.begin(); IT!= a_sequences.end();IT++){
			ss = (*IT);
			l = gmseq.lengths[i];

			memcpy(&(gmseq.seq[ii]), ss, l*sizeof(dna5_t));
			ii += l;

			i++;
		}
	}
	else{
		gmseq.seq_length = 0;
		for(usize_t i=0; i<gmseq.nof_seqs; i++){
			gmseq.lengths[i] = a_lengths.at(i) +1;
			gmseq.seq_length += gmseq.lengths[i];
		}

		gmseq.seq = new dna5_t[gmseq.seq_length];

		dna5_t *ss;
		int i=0, ii=0, l=0;
		for(std::vector<dna5_t*>::iterator IT = a_sequences.begin(); IT!= a_sequences.end();IT++){
			ss = (*IT);
			l = gmseq.lengths[i]-1;

			memcpy(&(gmseq.seq[ii]), ss, l*sizeof(dna5_t));
			ii += l;
			gmseq.seq[ii] = DNA5Alphabet::UNDEF;
			ii++;

			i++;
		}
	}


}


/*
 * An object to store suffix array data
 */
class CsFullyContainer{
public:
	usize_t *SA;
	usize_t *LCP;
	usize_t *NS;
	usize_t first_ncode; //first position in SA where we found a N character (as start of that suffix)
	usize_t first_ecode; //first position in SA where we found a $ character (as start of that suffix)
	bool freeOnExit;//if true then delete data on decostruction

	CsFullyContainer(
			usize_t *_SA,
			usize_t *_LCP,
			usize_t *_NS,
			usize_t _first_ncode,
			usize_t _first_ecode,
			bool _freeOnExit
			)
	: SA(_SA), LCP(_LCP), NS(_NS), first_ncode(_first_ncode), first_ecode(_first_ecode), freeOnExit(_freeOnExit){
	}
	CsFullyContainer(bool _freeOnExit){
		SA =NULL;
		LCP = NULL;
		NS = NULL;
		first_ncode = 0;
		first_ecode = 0;
		freeOnExit = _freeOnExit;
	}

	~CsFullyContainer(){
		if(freeOnExit){
			if(SA!=NULL) delete [] SA;
			if(LCP!=NULL) delete [] LCP;
			if(NS!=NULL) delete [] NS;
		}
	}

};




/*
 * template function to build a suffix array of a sequence (of undefined type). See sais.h, too.
 */
template <typename T>
usize_t*
build_SA(T *seq, usize_t length){
	//usize_t *SA;// = new usize_t[length];
	//SA = sais_build_SA(seq, length);
	return sais_build_SA(seq, length);
}


/*
 * print a suffix array and its suffixes (of length n)
 */
void
print_SA(dna5_t *seq/*indexed sequence*/, usize_t length, usize_t* SA, usize_t n){
	for(usize_t i=0; i<length; i++){
		std::cout<<i<<"\t";
		for(usize_t j=0; j<n; j++){
			std::cout<<DNA5Alphabet::symbolFor(seq[SA[i] + j]);
		}
		std::cout<<"\t"<<SA[i]<<"\n";
	}
}



/*
 * linear algorithm to biuld the LCP array.
 * LCP[i] = lcp( suffix(sa[i]), suffix(sa[i-1]) )
 */
template <typename T>
usize_t*
build_LCP(T *seq, usize_t length, usize_t *SA){
	usize_t *LCP = new usize_t[length];

	_si64* rank = new _si64 [length];

	usize_t ilength = (usize_t)length;

	for (_si64 i = 0; i < ilength; i++)
	  rank[SA[i]] = i;

	_si64 h = 0, k, j;

	for (_si64 i = 0; i < ilength; i++){
	  k = rank[i];
	  if (k == 0){
		  //lcp[k] = -1;
		  LCP[k] = 0;
	  }
	  else{
		  j = SA[k - 1];
		  while (i + h < ilength && j + h < ilength
			  && seq[i + h] == seq[j + h]){
			  h++;
		  }
		  LCP[k] = h;
	  }
	  if (h > 0) h--;
	}
	delete [] rank;

	return LCP;
}


/*
 * linear algorithm to build the NS array.
 * NS[i] = is the distance between the position SA[i] in seq and the closest next N.
 * For example, if seq[sa[i]]=N  => NS[i] = 0
 * if seq[sa[i]]|=N  and  seq[sa[i]+1]=N  => NS[i] = 1
 */
usize_t*
build_NS(dna5_t *seq, usize_t length, usize_t *SA){
	usize_t *NS = new usize_t[length];
	usize_t *_ns = new usize_t[length];
	usize_t last_n = length;
	dna5_t ncode = DNA5Alphabet::N;


	for(usize_t i=length-1; i>=0; i--){
		if(seq[i] >= ncode){
			last_n = i;
		}
		_ns[i] = last_n - i;

		if(i==0){
			break;
		}
	}

	for(usize_t i=0; i<length; i++){
		NS[i] = _ns[SA[i]];
	}

	delete [] _ns;
	return NS;
}


/*
 * print SA,LCP and NS data, plus suffixes of length n
 */
void
print_SA_LCP_N(dna5_t *seq/*indexed sequence*/, usize_t length, usize_t* SA, usize_t *LCP, usize_t *N, usize_t n){
	for(usize_t i=0; i<length; i++){
		std::cout<<i<<"\t";
		for(usize_t j=0; j<n; j++){
			std::cout<<DNA5Alphabet::symbolFor(seq[SA[i] + j]);
		}
		std::cout<<"\t"<<SA[i]<<"\t"<<LCP[i]<<"\t"<<N[i]<<"\n";
	}
}



/*
 * tell min_i(SA[i] = N)
 */
usize_t
where_first_ncode(usize_t *SA, usize_t *NS, usize_t length){
	usize_t i;
	for(i=length -1; i>0 && NS[i]==0; i--);
	return i +1;
}

/*
 * tell min_i(SA[i] = N)
 */
usize_t
where_first_ncode(dna5_t *seq, usize_t *SA, usize_t length){
	usize_t i;
	for(i=length -1; i>0 && seq[SA[i]] > DNA5Alphabet::T; i--);
	return i +1;
}

/*
 * tell min_i(SA[i] > N), namely  min_i(SA[i] = $)
 */
usize_t
where_first_ecode(dna5_t *seq, usize_t *SA, usize_t length){
	usize_t i;
	for(i=length -1; i>0 && seq[SA[i]] > DNA5Alphabet::N; i--);
	return i +1;
}


void
write_nelsa(
	std::string &ofile,
	usize_t iseq_length,
	usize_t *SA,
	usize_t *LCP,
	usize_t *NS
){
	std::ofstream off;
    off.open(ofile, std::ios::out | std::ios::binary);
    off.write(reinterpret_cast<const char *>(&iseq_length), sizeof(usize_t) );
    for(usize_t i=0; i<iseq_length; i++){
        off.write( reinterpret_cast<const char *>(&SA[i]), sizeof(usize_t) );
        off.write( reinterpret_cast<const char *>(&LCP[i]), sizeof(usize_t) );
        off.write( reinterpret_cast<const char *>(&NS[i]), sizeof(usize_t) );
    }
    off.flush();
    off.close();
}


void
read_nelsa(
	std::string &ifile,
	usize_t *length,
	usize_t *&SA,
	usize_t *&LCP,
	usize_t *&NS
){
	std::ifstream iff;
    iff.open(ifile, std::ios::in | std::ios::binary);

	iff.read(reinterpret_cast<char *>(length), sizeof(usize_t) );

	SA = new usize_t[*length];
	LCP = new usize_t[*length];
	NS = new usize_t[*length];

	std::cout<<"allocation done "<< length <<"\n";

	for(usize_t i=0; i<*length; i++){
		iff.read( reinterpret_cast<char *>(SA+i) ,sizeof(usize_t));
		iff.read( reinterpret_cast<char *>(LCP+i) ,sizeof(usize_t));
		iff.read( reinterpret_cast<char *>(NS+i) ,sizeof(usize_t));
	}


    iff.close();
}


}



#endif /* CSTYLEINDEX_H_ */
