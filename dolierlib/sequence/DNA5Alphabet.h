/*
 * DNA5Alphabet.h
 *
 *  Created on: Sep 26, 2013
 *      Author: vbonnici
 */


/*
 * Implementation of the 3bit alphabet {A,C,G,T,N,$}. It's exentially a bijective map from {A,C,G,T,N,$} to [0,6]
 */


#ifndef DNA5ALPHABET_H_
#define DNA5ALPHABET_H_


#include "data_ts.h"

namespace dolierlib{

class DNA5Alphabet {
public:

	static const dna5_t A = 0x0;
	static const dna5_t C = 0x1;
	static const dna5_t G = 0x2;
	static const dna5_t T = 0x3;
	static const dna5_t N = 0x4;
	static const dna5_t UNDEF = 0x5;

	static const symbol_t UNDEF_SYMBOL = '$';
	static const symbol_t SYMBOLS[6]; //map [0,5] in {A,C,G,T,N,$}, see at the end of the class.

	/*
	 * return a 3bit code for a symbol
	 */
	static dna5_t codeFor(symbol_t s){
		switch(s){
			case 'a':
			case 'A':
				return A;
			case 'c':
			case 'C':
				return C;
			case 'g':
			case 'G':
				return G;
			case 't':
			case 'T':
				return T;
			case 'n':
			case 'N':
				return N;
		}
		return UNDEF;
	}

	/* retun the complement of a dn5_t code
	 */
	static dna5_t complement(dna5_t s){
		switch(s){
			case A:
				return T;
			case C:
				return G;
			case G:
				return C;
			case T:
				return A;
			case N:
				return N;
		}
		return UNDEF;
	}

	/*
	 * return the symbol of a 3bit code
	 */
	static symbol_t symbolFor(dna5_t c){
		if(c<6)
			return SYMBOLS[c];
		return UNDEF_SYMBOL;
	}
	/*
	 * return the symbol of a 3bit code stored in a integer field
	 */
	static symbol_t symbolFor(int c){
		if(c<6)
			return SYMBOLS[c];
		return UNDEF_SYMBOL;
	}

	/*
	 * size of this alphabet. $ is not considered.
	 */
	static usize_t size(){
		return 5;
	}


	/*
	 * covert a 3bit sequence (of length n) in a char (symbol) array
	 */
	static void toSymbolArray(dna5_t *codes, char *syms, int n){
		for(int i=0; i<n; i++){
			syms[i] = symbolFor(codes[i]);
		}
	}

};

const symbol_t DNA5Alphabet::SYMBOLS[6] = {'A', 'C', 'G', 'T', 'N', UNDEF_SYMBOL};

}


#endif /* DNA5ALPHABET_H_ */
