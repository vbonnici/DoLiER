/*
 * BEMIndex.h
 *
 *  Created on: Oct 10, 2013
 *      Author: vbonnici
 */

/*
 * it's the first version of the index to deal with the document listing problem.
 * It's quite obsolete but faster.
 * It can store frequencies up to fixed values of k and m.
 */

#ifndef BEMINDEX_H_
#define BEMINDEX_H_

//#define DEBUG_BEMINDEX_H_s

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include "data_ts.h"

#include "timer.h"
#include "CstyleIndex.h"
#include "FASTAReader.h"
#include "DNA5Alphabet.h"

#include "Cs5DLIndex.h"

#include "main_common.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>


namespace dolierlib{

class BEMIndex{
public:
	static const usize_t version;

	usize_t min_k;
	usize_t max_k;
	usize_t max_miss;

	bool verse;
	bool reverse;

	usize_t nof_sequences;
	usize_t *seq_lengths;
	usize_t g_length;
	dna5_t* g_seq;

	//8-mer index


	usize_t SA_length;
	usize_t *SA;
	_ui8 *SA_mink;
	_ui8 *SA_maxk;
	usize_t **dcounts;

	bool freeOnExit;

	BEMIndex(){
		min_k = 0;
		max_k = 0;
		max_miss = 0;
		verse= false;
		reverse = false;
		nof_sequences = 0;
		seq_lengths = NULL;
		g_length = 0;
		g_seq = NULL;
		SA_length = 0;
		SA = NULL;
		SA_mink = NULL;
		SA_maxk = NULL;
		dcounts = NULL;
		freeOnExit = true;

		b_SA = NULL;
		b_dlindex = NULL;
	}
	~BEMIndex(){
#ifdef DEBUG_BEMINDEX_H_
			std::cout<<"BEMIndex::deco>\n";
#endif
		if(freeOnExit){
			if(seq_lengths != NULL)
				delete [] seq_lengths;
			if(g_seq != NULL)
				delete [] g_seq;
			if(SA != NULL)
				delete [] SA;
			if(SA_mink != NULL)
				delete [] SA_mink;
			if(SA_maxk != NULL)
				delete [] SA_maxk;
			if(dcounts != NULL){
				for(usize_t i=0; i<SA_length; i++)
					delete [] dcounts[i];
				delete [] dcounts;
			}
		}
#ifdef DEBUG_BEMINDEX_H_
			std::cout<<"BEMIndex::deco<\n";
#endif
	}


	bool check(usize_t from_k, usize_t to_k, usize_t upto_m, bool is_verse, bool is_reverse){
		bool ret = true;
		if(from_k < min_k){ std::cout<<"BEM index do not support k = "<<from_k<<", minimum value = "<<min_k<<"\n";  ret = false;}
		if(to_k > max_k){ std::cout<<"BEM index do not support k = "<<to_k<<", maximum value = "<<max_k<<"\n";  ret = false;}
		if(upto_m > max_miss){ std::cout<<"BEM index do not support m = "<<upto_m<<", maximum value = "<<max_miss<<"\n";  ret = false;}
		if(is_verse && !verse){ std::cout<<"BEM index does not store verse statistics\n";  ret = false;}
		if(!is_verse && reverse){ std::cout<<"BEM index stores reverse complement statistics\n";  ret = false;}
		if(is_reverse && !reverse){ std::cout<<"BEM index does not store reverse complement statistics\n";  ret = false;}
		if(!is_reverse && reverse){ std::cout<<"BEM index stores reverse complement statistics\n";  ret = false;}
		return ret;
	}

	bool check(usize_t from_k, usize_t to_k, usize_t upto_m){
		bool ret = true;
		if(from_k < min_k){ std::cout<<"BEM index do not support k = "<<from_k<<", minimum value = "<<min_k<<"\n";  ret = false;}
		if(to_k > max_k){ std::cout<<"BEM index do not support k = "<<to_k<<", maximum value = "<<max_k<<"\n";  ret = false;}
		if(upto_m > max_miss){ std::cout<<"BEM index do not support m = "<<upto_m<<", maximum value = "<<max_miss<<"\n";  ret = false;}
		return ret;
	}


	bool load(std::string file){
		std::ifstream ifs;
		ifs.open(file.c_str(), std::ios::in);
		if(!ifs.good()){
			std::cout<<"BEMIndex error on opening file "<<file<<"!!!\n";
			return false;
		}
		bool ret = load(ifs, *this);
		ifs.close();
		std::cout<<"done\n";
		return ret;
	}

	static
	bool load(std::istream& is, BEMIndex& index){
		usize_t vers = 0;
		is.read((char*)&vers,sizeof(usize_t));
		if(vers != version){
			std::cout<<"BEMIndex wrong version!!!\n";
			return false;
		}

		/**
		 * FORMAT:
		 * 		version (usize_t)
		 * 		min_k (usize_t)
		 * 		max_k (usize_t)
		 * 		max_miss (usize_t)
		 * 		verse (bool)
		 * 		reverse (bool)
		 * 		SEQUENCE
		 * 			nof_sequences (usize_t)
		 * 			seq_lenghs	(usize_t*)
		 * 			length (usize_t)
		 * 			data	(dna5_t*)
		 *
		 * 		8-MER INDEX (implict length = 4^8 = 65536)
		 * 			fseek(strictSA_start) (usize_t)
		 * 			//0 if absent
		 *
		 *		STRICT SA
		 * 			StrictSA length	(usize_t)
		 * 			StrictSA + dcounts
		 * 				SA (usize_t)	mink (_ui8)	maxk (_ui8)		dcounts (usize_t*)
		 */

		is.read((char*)&(index.min_k),sizeof(usize_t));
//		std::cout<<index.min_k<<"\n";
		is.read((char*)&(index.max_k),sizeof(usize_t));
//		std::cout<<index.max_k<<"\n";
		is.read((char*)&(index.max_miss),sizeof(usize_t));
//		std::cout<<index.max_miss<<"\n";
		is.read((char*)&(index.verse),sizeof(bool));
//		std::cout<<index.verse<<"\n";
		is.read((char*)&(index.reverse),sizeof(bool));
//		std::cout<<index.reverse<<"\n";


		is.read((char*)&(index.nof_sequences),sizeof(usize_t));
//		std::cout<<index.nof_sequences<<"\n";

		index.seq_lengths = new usize_t[(index.nof_sequences)];
		for(usize_t i=0; i<(index.nof_sequences); i++){
			is.read((char*)&(index.seq_lengths[i]),sizeof(usize_t));
//			std::cout<<index.seq_lengths[i]<<"\t";
		}
		std::cout<<"\n";

		is.read((char*)&(index.g_length),sizeof(usize_t));
//		std::cout<<index.g_length<<"\n";
		index.g_seq = new dna5_t[index.g_length];
		for(usize_t i=0; i<index.g_length; i++)
			is.read((char*)&(index.g_seq[i]),sizeof(dna5_t));

//		usize_t fake = 0;
//		for(int i=0; i<65536+1; i++)
//			is.read((char*)&(fake), sizeof(usize_t));


		is.read((char*)&(index.SA_length),sizeof(usize_t));
//		std::cout<<index.SA_length<<"\n";
		index.SA = new usize_t[index.SA_length];
		index.SA_mink = new _ui8[index.SA_length];
		index.SA_maxk = new _ui8[index.SA_length];
		index.dcounts = new usize_t*[index.SA_length];


		usize_t tmp = 0;

		for(usize_t i=0; i<index.SA_length; i++){
			is.read((char*)&(index.SA[i]),sizeof(usize_t));
//			std::cout<<index.SA[i]<<"\t";
			is.read((char*)&(index.SA_mink[i]),sizeof(_ui8));
//			std::cout<<(int)index.SA_mink[i]<<"\t";
			is.read((char*)&(index.SA_maxk[i]),sizeof(_ui8));
//			std::cout<<(int)index.SA_maxk[i]<<"\t";
			tmp = (index.SA_maxk[i]-index.SA_mink[i]+1) * (index.max_miss+1);
//			std::cout<<tmp<<"\t";
			index.dcounts[i] = new usize_t[tmp];
			for(usize_t j=0; j<tmp;j++){
				is.read((char*)&(index.dcounts[i][j]),sizeof(usize_t));
//				std::cout<<index.dcounts[i][j]<<"\t";
			}
//			std::cout<<"\n";
		}

//		std::cout<<"done\n";

		return true;
	}



#define BEMIndex_codeAt(pos, shift) g_seq[SA[pos] + shift]

	usize_t find_ii_end;
	dna5_t find_cs;
	usize_t find_m;

	/**
	 * return i_end+1 if not found
	 */
	inline
	usize_t find_first(dna5_t c, usize_t i_start, usize_t i_end, usize_t shift){
//std::cout<<"fin_first("<<i_start<<","<<i_end<<")\n";
//		usize_t find_ii_end = i_end;
//		dna5_t find_cs = CsDLIndex_codeAt(i_start, shift);
		find_ii_end = i_end;
		find_cs = BEMIndex_codeAt(i_start, shift);

		if(find_cs > c){
			return i_end+1;
		}
		find_cs = BEMIndex_codeAt(i_end, shift);
		if(find_cs < c){
			return i_end+1;
		}

//		usize_t find_m;
		while(true){
			if(i_start == i_end){
				break;
			}
			find_m = ((i_end - i_start)/2) + i_start;

			find_cs = BEMIndex_codeAt(find_m, shift);
			if(c <= find_cs){
				i_end = find_m;
			}
			else{
				i_start = find_m+1;
			}

		}
		if(c != BEMIndex_codeAt(i_start, shift)){
			return find_ii_end+1;
		}

		return i_start;

	}
	/**
	 * return i_end+1 if not found
	 */
	inline
	usize_t find_last(dna5_t c, usize_t i_start, usize_t i_end, usize_t shift){
//		usize_t find_ii_end = i_end;
//		dna5_t find_cs = CsDLIndex_codeAt(i_start, shift);
		find_ii_end = i_end;
		find_cs = BEMIndex_codeAt(i_start, shift);

		if(find_cs > c){
			return i_end+1;
		}
		find_cs = BEMIndex_codeAt(i_end, shift);
		if(find_cs < c){
			return i_end+1;
		}

//		usize_t find_m;
		while(true){
			if(i_start == i_end){
				break;
			}
			find_m = ((i_end - i_start)/2) + i_start;

			if((i_end - i_start)%2 != 0)
				find_m++;

			find_cs = BEMIndex_codeAt(find_m, shift);
			if(c < find_cs){
				i_end = find_m-1;
			}
			else{
				i_start = find_m;
			}

		}
		if(c != BEMIndex_codeAt(i_start, shift))
			return find_ii_end+1;
		return i_start;
	}




private:
	bool
	count_docs(dna5_t *pattern, usize_t p_length, usize_t mismatches, usize_t *count){
		if(mismatches > max_miss)	return 0;
		if(p_length<min_k || p_length>max_k) return 0;

		usize_t c_start = 0;
		usize_t c_end = SA_length-1;

		usize_t t_start;
		usize_t t_end;

		usize_t last_si = p_length-1;

		for(usize_t si=0; si<p_length; si++){
//			std::cout<<DNA5Alphabet::symbolFor(pattern[si])<<"\n";
			t_start = find_first(pattern[si], c_start, c_end, si);
			if(t_start != c_end+1){
//				std::cout<<si<<" t_start\t";print_dna5_SA(g_seq, SA, t_start, /*SA_maxk[i]*/ max_k);std::cout<<"\t"<<SA[t_start]<<"\n";

				if(si == last_si){
					if(p_length<SA_mink[t_start] || p_length>SA_maxk[t_start])
						std::cout<<"warning wrong SA_mink,SA_maxk";

					*count = dcounts[t_start][  ((p_length-SA_mink[t_start])*(max_miss+1)) + mismatches  ];

//					if(t_start>0)
//					print_dna5_SA(g_seq, SA, t_start-1, /*SA_maxk[i]*/ max_k);std::cout<<"\n";
//					print_dna5_SA(g_seq, SA, t_start, /*SA_maxk[i]*/ max_k);std::cout<<"\t";
//					std::cout<<SA[t_start]<<"\t";
//					std::cout<<(int)SA_mink[t_start]<<"\t";
//					std::cout<<(int)SA_maxk[t_start]<<"\t";
//					usize_t tmp = (SA_maxk[t_start]-SA_mink[t_start]+1) * (max_miss+1);
//					std::cout<<tmp<<"\t";
//					for(usize_t j=0; j<tmp;j++){
//						std::cout<<dcounts[t_start][j]<<"\t";
//					}
//					std::cout<<"\n";
//					if(t_start<SA_length-1)
//					print_dna5_SA(g_seq, SA, t_start+1, /*SA_maxk[i]*/ max_k);std::cout<<"\n";
				}
				else{
					t_end = find_last(pattern[si], t_start, c_end, si);
//					std::cout<<si<<" t_end \t";print_dna5_SA(g_seq, SA, t_end, /*SA_maxk[i]*/ max_k);std::cout<<"\t"<<SA[t_end]<<"\n";

					c_start = t_start;
					c_end = t_end;
				}
			}
			else{
				return false;
			}
		}

		return true;
	}


private:
	usize_t *b_SA;
	Cs5DLIndex *b_dlindex;
public:
	bool
	count_docs_nf(dna5_t *pattern, usize_t p_length, usize_t mismatches, usize_t *count){
		bool ret = count_docs(pattern,p_length,mismatches,count);
		if(!ret){
			if(b_dlindex == NULL){
				b_SA = build_SA<dna5_t>(g_seq, g_length);
				b_dlindex = new Cs5DLIndex(g_seq, g_length,b_SA, seq_lengths, nof_sequences);
			}
			bool* b_idocs = (bool*)calloc( nof_sequences, sizeof(bool));
			dolierlib::ssize_t *D = new dolierlib::ssize_t[p_length];


			if(verse){
				b_dlindex->dl_count_wD(pattern, D, p_length, mismatches, b_idocs);

			}
			if(!verse || reverse){
				dna5_t *it_codes_rev = new dna5_t[p_length];
				reverseComplement(pattern, it_codes_rev, p_length);
				if(!equal_dna5(pattern, it_codes_rev, p_length)){
					b_dlindex->dl_count_wD(it_codes_rev, D, p_length, mismatches, b_idocs);
				}
				delete [] it_codes_rev;
			}


			*count = 0;
			for(usize_t idoc_i=0; idoc_i<nof_sequences; idoc_i++)
				if(b_idocs[idoc_i]) (*count)++;
			free(b_idocs);
			delete [] D;
		}
	}




	void print(){
		std::cout<<min_k<<"\n";
		std::cout<<max_k<<"\n";
		std::cout<<max_miss<<"\n";
		std::cout<<verse<<"\n";
		std::cout<<reverse<<"\n";


		std::cout<<nof_sequences<<"\n";

		for(usize_t i=0; i<(nof_sequences); i++){
			std::cout<<seq_lengths[i]<<"\t";
		}
		std::cout<<"\n";

		std::cout<<g_length<<"\n";

		std::cout<<SA_length<<"\n";

		usize_t tmp = 0;

		std::cout<<"i\tdna\tSA[i]\tmink[i]\tmaxk[i]\tdsize\tdcounts\n";

		for(usize_t i=0; i<SA_length; i++){
			std::cout<<i<<"\t";
			print_dna5_SA(g_seq, SA, i, /*SA_maxk[i]*/ max_k);std::cout<<"\t";
			std::cout<<SA[i]<<"\t";
			std::cout<<(int)SA_mink[i]<<"\t";
			std::cout<<(int)SA_maxk[i]<<"\t";
			tmp = (SA_maxk[i]-SA_mink[i]+1) * (max_miss+1);
			std::cout<<tmp<<"\t";
			for(usize_t j=0; j<tmp;j++){
				std::cout<<dcounts[i][j]<<"\t";
			}
			std::cout<<"\n";
		}

		std::cout<<"done\n";
	}






	class iterator{
	public:
		BEMIndex &_index;
		usize_t _i;
		usize_t _k;
	public:
		iterator(BEMIndex &index, usize_t i, usize_t k) : _index(index), _i(i), _k(k){
		}

		iterator(const BEMIndex::iterator &it) : _index(it._index), _i(it._i), _k(it._k){};

		bool next(){
			_i++;
			while(_i<=_index.SA_length && (_index.SA_mink[_i-1]>_k || _index.SA_maxk[_i-1]<_k)){
				_i++;
			}
			return _i <= _index.SA_length;
		}

		usize_t freq(usize_t miss){
			return _index.dcounts[_i-1][  ((_k - _index.SA_mink[_i-1])*(_index.max_miss+1)) + miss  ];
		}

		void get_kmer(dna5_t *kmer){
			memcpy(kmer, &(_index.g_seq[_index.SA[_i-1]]), _k*sizeof(dna5_t));
		}
	};

	BEMIndex::iterator begin(usize_t k){
		return BEMIndex::iterator(*this, 0, k);
	}
















	static
	void write(
			const std::string	&o_file,
			usize_t min_k,
			usize_t max_k,
			usize_t max_miss,
			bool verse,
			bool reverse,
			const DNA5MS_t &f_cs_gseq,
			const CsFullyContainer &f_container,
			const int *mink_table,
			const int *maxk_table,
			const usize_t * const *dcount,
			usize_t strictSA_size
			){
		std::cout<<"writing...\n";

		/**
		 * saving index
		 */
		std::ofstream os;
		os.open((o_file).c_str(), std::ios::out);
		if(!os.is_open() || os.bad()){
			std::cout<<"Error on opening output file : "<<o_file<<" \n";
			exit(1);
		}



		/**
		 * FORMAT:
		 * 		version (usize_t)
		 * 		min_k (usize_t)
		 * 		max_k (usize_t)
		 * 		max_miss (usize_t)
		 * 		verse (bool)
		 * 		reverse (bool)
		 * 		SEQUENCE
		 * 			nof_sequences (usize_t)
		 * 			seq_lenghs	(usize_t*)
		 * 			length (usize_t)
		 * 			data	(dna5_t*)
		 *
		 *		(removed)
		 * 		8-MER INDEX (implict length = 4^8 = 65536)
		 * 			fseek(strictSA_start) (usize_t)
		 * 			//0 if absent
		 *
		 *		STRICT SA
		 * 			StrictSA length	(usize_t)
		 * 			StrictSA + dcounts
		 * 				SA (usize_t)	mink (_ui8)	maxk (_ui8)		dcounts (usize_t*)
		 */

		usize_t version = BEMIndex::version;

		//version (usize_t)
		os.write((char*)&version, sizeof(usize_t));
	//	std::cout<<version<<"\n";

		//min_k (usize_t)
		os.write((char*)&min_k, sizeof(usize_t));
	//	std::cout<<min_k<<"\n";
		//max_k (usize_t)
		os.write((char*)&max_k, sizeof(usize_t));
	//	std::cout<<max_k<<"\n";
		//max_miss (usize_t)
		os.write((char*)&max_miss, sizeof(usize_t));
	//	std::cout<<max_miss<<"\n";

		//verse (bool)
		os.write((char*)&verse, sizeof(bool));
	//	std::cout<<verse<<"\n";
		//reverse (bool)
		os.write((char*)&reverse, sizeof(bool));
	//	std::cout<<reverse<<"\n";

		//SEQUENCE
		//nof_sequences (usize_t)
		os.write((char*)&(f_cs_gseq.nof_seqs), sizeof(usize_t));
	//	std::cout<<f_cs_gseq.nof_seqs<<"\n";
		//seq_lenghs	(usize_t*)
		for(usize_t i=0; i<f_cs_gseq.nof_seqs; i++){
			os.write((char*)&(f_cs_gseq.lengths[i]), sizeof(usize_t));
	//		std::cout<<f_cs_gseq.lengths[i]<<"\t";
		}
	//	std::cout<<"\n";
		//length
		os.write((char*)&(f_cs_gseq.seq_length), sizeof(usize_t));
	//	std::cout<<seq_length<<"\n";
		//data	(dna5_t*)
		for(usize_t i=0; i<f_cs_gseq.seq_length; i++)
			os.write((char*)&(f_cs_gseq.seq[i]), sizeof(dna5_t));


		usize_t fake = 0;
	//	//8-MER INDEX (implict length = 4^8 = 65536)
	//	usize_t k8index_p = os.tellp();
	//
	//	for(int i=0; i<65536+1; i++)
	//		os.write((char*)&(fake), sizeof(usize_t));
	//	//it will be filled later
	//	_ui16 *k8index = new _ui16[65536+1];memset(k8index,0,(65536+1)*sizeof(_ui16));


	//	int code4_length = sizeof(code4_t)*4;
	//	code4_t c4 = 0;


		_ui8 tmp;
		//STRICT SA
		//StrictSA length	(usize_t)
		os.write((char*)&strictSA_size, sizeof(usize_t));
	//	std::cout<<strictSA_size<<"\n";
		for(usize_t i=0; i<f_cs_gseq.seq_length; i++){
			//we still have to mask for the real strict SA
			if(mink_table[i] != -1){
				//k8index start interval
	//			fill_code4_from_SA(f_cs_gseq.seq, f_container.SA, i, &c4);
	//			std::cout<<c4<<"\n";
	//			if(k8index[c4]==0){
	//				k8index[c4] = os.tellp();
	//			}

				//SA (usize_t)
				os.write((char*)&(f_container.SA[i]), sizeof(usize_t));
	//			std::cout<<f_container.SA[i]<<"\t";
				//mink (usize_t?int)
				tmp = (_ui8)mink_table[i];
				os.write((char*)&(tmp), sizeof(_ui8));
	//			std::cout<<(int)tmp<<"\t";
				//maxk (usize_t)
				tmp = (_ui8)maxk_table[i];
	//			std::cout<<(int)tmp<<"\t";
				os.write((char*)&(tmp), sizeof(_ui8));
				//dcounts (usize_t*)
				fake = (mink_table[i]==-1 ? 0: maxk_table[i]-mink_table[i]+1);
	//			std::cout<<(fake*(max_miss+1))<<"\t";
				for(usize_t j=0; j<fake * (max_miss+1); j++){
					os.write((char*)&(dcount[i][j]), sizeof(usize_t));
	//				std::cout<<dcount[i][j]<<"\t";
				}
	//			std::cout<<"\n";
			}
		}
		//k8index[65536-1] = os.tellp();
	//	k8index[65536] = os.tellp();
	//
	//
	//	//re-writing k8index
	//	os.seekp(k8index_p);
	//	for(int i=0; i<65536; i++)
	//			os.write((char*)&(k8index[i]), sizeof(usize_t));

		os.flush();
		os.close();
	}



	static
	void create(
			std::string	seqs_file,
			std::string	o_file,
			usize_t min_k,
			usize_t max_k,
			usize_t max_miss,
			bool verse,
			bool reverse
			){


		TIMEHANDLE total_t = start_time();
		TIMEHANDLE start_t;

		start_t = start_time();

		std::cout<<"loading sequence...\n";
		std::vector<dna5_t*> f_cs_sequences;
		std::vector<usize_t> f_cs_lengths;
		FASTAReader f_cs_reader(seqs_file);
		f_cs_reader.readAll(f_cs_sequences, f_cs_lengths, true);
		f_cs_reader.close();

		std::cout<<"concatenating...\n";
		DNA5MS_t f_cs_gseq;
		concat(f_cs_gseq, f_cs_sequences, f_cs_lengths, true);
		for(std::vector<dna5_t*>::iterator IT = f_cs_sequences.begin(); IT!=f_cs_sequences.end();IT++)
			delete [] (*IT);

		std::cout<<"building ds...\n";

		CsFullyContainer f_container(true);
		f_container.SA = build_SA<dna5_t>(f_cs_gseq.seq, f_cs_gseq.seq_length);
		f_container.LCP = build_LCP<dna5_t>(f_cs_gseq.seq, f_cs_gseq.seq_length, f_container.SA);
		f_container.NS = build_NS(f_cs_gseq.seq, f_cs_gseq.seq_length, f_container.SA);
		f_container.first_ncode = where_first_ncode(f_container.SA, f_container.NS, f_cs_gseq.seq_length);

		Cs5DLIndex f_cs_dl5index(f_cs_gseq.seq, f_cs_gseq.seq_length, f_container.SA, f_cs_gseq.lengths, f_cs_gseq.nof_seqs);

		std::cout<<"cs build done: "<<end_time(start_t)<<"\n";
		std::cout<<"nof sequences:"<<f_cs_gseq.nof_seqs<<"\n";
		std::cout<<"total length: "<<f_cs_gseq.seq_length<<"\n";


	//	print_SA(f_cs_gseq.seq, f_cs_gseq.seq_length, f_container.SA,max_k+1);




		usize_t seq_length = f_cs_gseq.seq_length;




		/*
		 * mink_table  maxk_table
		 */
		int *mink_table = new int[seq_length];memset(mink_table,-1,seq_length*sizeof(int));
		int *maxk_table = new int[seq_length];memset(maxk_table,-1,seq_length*sizeof(int));


		usize_t size_k = 0;
		for(usize_t k=min_k; k<=max_k; k++){
			dna5_t *cs_it_codes = new dna5_t[k];
			dolierlib::ssize_t *D = new dolierlib::ssize_t[k];

			size_k = 0;

			NSAIterator cs_iterator = NSAIterator::begin(f_cs_gseq, f_container, k);
			while(cs_iterator.next()){
				cs_iterator.get_kmer(cs_it_codes);
	//			print_dna5(cs_it_codes,k);std::cout<<"\t"<<cs_iterator.i_start<<"\n";
				size_k++;


				if(mink_table[cs_iterator.i_start] == -1){
					mink_table[cs_iterator.i_start] = k;
					maxk_table[cs_iterator.i_start] = k;
				}
				else{
					if(k<mink_table[cs_iterator.i_start])	mink_table[cs_iterator.i_start] = k;
					if(k>maxk_table[cs_iterator.i_start])	maxk_table[cs_iterator.i_start] = k;
				}
			}

			std::cout<<k<<"-mers: "<<size_k<<"\n";

			delete [] cs_it_codes;
			delete [] D;
		}



		/**
		 * Strict SA
		 */
		//size
		//mapSA
		usize_t strictSA_size = 0;
		usize_t *mapSA = new usize_t[seq_length];

		for(usize_t i=0, s=0; i<seq_length; i++){
			if(mink_table[i] != -1){
				strictSA_size++;
				mapSA[i] = s;
				s++;
			}
		}

		std::cout<<"strict SA size = "<<strictSA_size<<"\n";


		/**
		 * DCount table
		 */
	//	usize_t ***dcount = new usize_t**[seq_length];
	//	for(usize_t i=0; i<seq_length; i++){
	//
	//		int x = mink_table[i]==-1 ? 0: maxk_table[i]-mink_table[i]+1;
	//
	//		dcount[i] = new usize_t*[ x ];
	//		//dcount[i] = new usize_t*[ max_k-min_k+1 ];
	//
	//		std::cout<<"dcount["<<i<<"] = "<<x<<"; "<<mink_table[i]<<","<<maxk_table[i]<<"\n";
	//
	//		for(usize_t j=0; j< x; j++)
	//		//for(usize_t j=0; j<max_k-min_k+1; j++)
	//			dcount[i][j] = new usize_t[max_miss];
	//	}

		usize_t **dcount = new usize_t*[seq_length];
		for(usize_t i=0; i<seq_length; i++){

			int x = mink_table[i]==-1 ? 0: maxk_table[i]-mink_table[i]+1;

			dcount[i] = new usize_t[ x *  (max_miss+1)];
			memset(dcount[i],0, (x *  (max_miss+1))*sizeof(usize_t));
			//dcount[i] = new usize_t[ (max_k-min_k+1) *  (max_miss+1)];

			//dcount[i] = new usize_t*[ max_k-min_k+1 ];

	//		std::cout<<"dcount["<<i<<"] = "<<x<<"; "<<(x *  max_miss)<<"; "<<mink_table[i]<<","<<maxk_table[i]<<"\n";
	//
	//		for(usize_t j=0; j< x; j++)
	//		//for(usize_t j=0; j<max_k-min_k+1; j++)
	//			dcount[i][j] = new usize_t[max_miss];
		}



		usize_t f_dcount,idoc_i;
		for(usize_t k=min_k; k<=max_k; k++){
			dna5_t *cs_it_codes = new dna5_t[k];
			dna5_t *cs_it_codes_rev = new dna5_t[k];
			dolierlib::ssize_t *D = new dolierlib::ssize_t[k];



	//		NSAIterator<dna5_t> cs_iterator = cs_begin<dna5_t>(f_cs_gseq, f_container, k);
	//		while(cs_iterator.next()){
	//			cs_iterator.get_kmer(cs_it_codes);
	//
	//			bool* f_idocs = (bool*)calloc(f_cs_gseq.nof_seqs, sizeof(bool));
	//
	//			for(idoc_i == cs_iterator.i_start; idoc_i <= cs_iterator.i_end-1; idoc_i++)
	//				f_idocs[ f_cs_dl5index.ss_ids[ idoc_i ] ] = true;
	//
	//			f_dcount = 0;
	//			for(idoc_i=0; idoc_i<f_cs_gseq.nof_seqs; idoc_i++){
	//				if(f_idocs[idoc_i]) f_dcount++;
	//			}
	//
	//			dcount[ cs_iterator.i_start ][ ((k-mink_table[cs_iterator.i_start])*(max_miss+1)) ] = f_dcount;
	//
	//
	//			free(f_idocs);
	//		}








	//		NSAIterator cs_iterator = NSAIterator::begin(f_cs_gseq, f_container, k);
			for(usize_t m=0; m<=max_miss; m++){

				std::cout<<k<<"\t"<<m<<"\t";

				start_t = start_time();

				//NSAIterator<dna5_t> cs_iterator = cs_begin<dna5_t>(f_cs_gseq, f_container, k);
				NSAIterator cs_iterator = NSAIterator::begin(f_cs_gseq, f_container, k);
				while(cs_iterator.next()){
					cs_iterator.get_kmer(cs_it_codes);

					bool* f_idocs = (bool*)calloc(f_cs_gseq.nof_seqs, sizeof(bool));



					if(verse)
						f_cs_dl5index.dl_count_wD(cs_it_codes, D, k, m, f_idocs);


					if(reverse){
						reverseComplement(cs_it_codes, cs_it_codes_rev, k);
						if(!verse || !equal_dna5(cs_it_codes, cs_it_codes_rev, k))
							f_cs_dl5index.dl_count_wD(cs_it_codes_rev, D, k, m, f_idocs);
					}


					f_dcount = 0;
					for(idoc_i=0; idoc_i<f_cs_gseq.nof_seqs; idoc_i++){
						if(f_idocs[idoc_i]) f_dcount++;
					}

	//				std::cout<<f_container.SA[cs_iterator.i_start]<<"\t";print_dna5(cs_it_codes,k);std::cout<<"\t"<<m<<"\t"<<f_dcount<<"\n";


	//				print_dna5(cs_it_codes,k);std::cout<<"\t"<<cs_iterator.i_start<<"\t"<<f_dcount<<"\n";
	//				std::cout<<"i_start = "<<cs_iterator.i_start<<"\n";
	//				std::cout<<"k = "<<k<<"\n";
	//				std::cout<<"m = "<<m<<"\n";
	//				std::cout<<"mink = "<<mink_table[cs_iterator.i_start]<<"\n";
	////				std::cout<<"map = "<< mapSA[ cs_iterator.i_start ] <<"\n";
	//				std::cout<<"K[] = "<<(   ((k-mink_table[cs_iterator.i_start])*max_miss) + m  )<<"\n";

					dcount[ cs_iterator.i_start ][  ((k-mink_table[cs_iterator.i_start])*(max_miss+1)) + m] = f_dcount;
					//dcount[ cs_iterator.i_start ][  (k*(max_miss)) + m] = f_dcount;

	//				GET2D(dcount[ mapSA[ cs_iterator.i_start ] ], max_miss, k, m) = f_dcount;


					free(f_idocs);
				}

				std::cout<<end_time(start_t)<<"\n";
			}

			delete [] cs_it_codes;
			delete [] cs_it_codes_rev;
			delete [] D;
		}


		write(o_file, min_k, max_k, max_miss, verse, reverse, f_cs_gseq, f_container,mink_table, maxk_table, dcount, strictSA_size);
		std::cout<<"done "<<end_time(total_t)<<"\n";




		delete [] mink_table;
		delete [] maxk_table;
		delete [] mapSA;
		for(usize_t i=0; i<seq_length; i++)
			delete [] dcount[i];
		delete [] dcount;
	}



private:
class SearcherThread{
public:
	const usize_t i_start;
	const usize_t i_end;
	const usize_t max_miss;
	const bool verse;
	const bool reverse;
	const DNA5MS_t &f_cs_gseq;
	const CsFullyContainer &f_container;
	Cs5DLIndex *f_cs_dl5index;
	const int *mink_table;
	const int *maxk_table;
	usize_t **dcount;
	double **timers;

	SearcherThread(
			const usize_t _i_start,
			const usize_t _i_end,
			const usize_t _max_miss,
			const bool _verse,
			const bool _reverse,
			const DNA5MS_t &_f_cs_gseq,
			const CsFullyContainer &_f_container,
			const int *_mink_table,
			const int *_maxk_table,
			usize_t **_dcount,
			double **_timers
			)
		: i_start(_i_start), i_end(_i_end),
		  max_miss(_max_miss), verse(_verse),reverse(_reverse),
		  f_cs_gseq(_f_cs_gseq), f_container(_f_container),
		  mink_table(_mink_table), maxk_table(_maxk_table),
		  dcount(_dcount), timers(_timers)
	{
		f_cs_dl5index = new Cs5DLIndex (f_cs_gseq.seq, f_cs_gseq.seq_length, f_container.SA, f_cs_gseq.lengths, f_cs_gseq.nof_seqs);
	}

	static void* run(void* argsptr){
		SearcherThread* st = (SearcherThread*)argsptr;

		TIMEHANDLE total_t = start_time();
		TIMEHANDLE start_t;
		int count = 0;

		usize_t f_dcount=0,idoc_i=0;

		for(usize_t i = st->i_start; i<=st->i_end; i++){
			if(st->mink_table[i] != -1){
				count++;

				for(usize_t k = st->mink_table[i]; k <= st->maxk_table[i]; k++){

					dna5_t *cs_it_codes = new dna5_t[k];
					memcpy(cs_it_codes, &(st->f_cs_gseq.seq[  st->f_container.SA[i]]), k*sizeof(dna5_t));

					dna5_t *cs_it_codes_rev = new dna5_t[k];
					dolierlib::ssize_t *D = new dolierlib::ssize_t[k];



					for(int m=0; m<=st->max_miss; m++){

						start_t = start_time();
						bool* f_idocs = (bool*)calloc(st->f_cs_gseq.nof_seqs, sizeof(bool));



						if(st->verse)
							st->f_cs_dl5index->dl_count_wD(cs_it_codes, D, k, m, f_idocs);


						if(st->reverse){
							reverseComplement(cs_it_codes, cs_it_codes_rev, k);
							if(!st->verse || !equal_dna5(cs_it_codes, cs_it_codes_rev, k))
								st->f_cs_dl5index->dl_count_wD(cs_it_codes_rev, D, k, m, f_idocs);
						}


						f_dcount = 0;
						for(idoc_i=0; idoc_i<st->f_cs_gseq.nof_seqs; idoc_i++){
							if(f_idocs[idoc_i]) f_dcount++;
						}

						st->dcount[ i ][  ((k-st->mink_table[i])*(st->max_miss+1)) + m] = f_dcount;

						free(f_idocs);

						st->timers[k][m] += end_time(start_t);
					}

					delete [] cs_it_codes;
					delete [] cs_it_codes_rev;
					delete [] D;
				}

			}
		}

		std::cout<<"-t "<<end_time(total_t)<<"\t"<<count<<"\n";

		pthread_exit(NULL);
	}
};



public:

	static
	void create_pp(
			std::string	seqs_file,
			std::string	o_file,
			usize_t min_k,
			usize_t max_k,
			usize_t max_miss,
			bool verse,
			bool reverse,
			int nthreads
			){

		TIMEHANDLE total_t = start_time();
		TIMEHANDLE start_t;

		start_t = start_time();

		std::cout<<"loading sequence...\n";
		std::vector<dna5_t*> f_cs_sequences;
		std::vector<usize_t> f_cs_lengths;
		FASTAReader f_cs_reader(seqs_file);
		f_cs_reader.readAll(f_cs_sequences, f_cs_lengths, true);
		f_cs_reader.close();

		std::cout<<"concatenating...\n";
		DNA5MS_t f_ms;
		concat(f_ms, f_cs_sequences, f_cs_lengths, true);
		for(std::vector<dna5_t*>::iterator IT = f_cs_sequences.begin(); IT!=f_cs_sequences.end();IT++)
			delete [] (*IT);

		std::cout<<"building ds...\n";

		CsFullyContainer f_ds(true);
		f_ds.SA = build_SA<dna5_t>(f_ms.seq, f_ms.seq_length);
		f_ds.LCP = build_LCP<dna5_t>(f_ms.seq, f_ms.seq_length, f_ds.SA);
		f_ds.NS = build_NS(f_ms.seq, f_ms.seq_length, f_ds.SA);
		f_ds.first_ncode = where_first_ncode(f_ds.SA, f_ds.NS, f_ms.seq_length);
		f_ds.first_ecode = where_first_ecode(f_ms.seq, f_ds.SA, f_ms.seq_length);

		Cs5DLIndex f_cs_dl5index(f_ms.seq, f_ms.seq_length, f_ds.SA, f_ms.lengths, f_ms.nof_seqs);

		std::cout<<"cs build done: "<<end_time(start_t)<<"\n";
		std::cout<<"nof sequences:"<<f_ms.nof_seqs<<"\n";
		std::cout<<"total length: "<<f_ms.seq_length<<"\n";


	//	print_SA(f_cs_gseq.seq, f_cs_gseq.seq_length, f_container.SA,max_k+1);

		usize_t seq_length = f_ms.seq_length;




		/*
		 * mink_table  maxk_table
		 */
		int *mink_table = new int[seq_length];memset(mink_table,-1,seq_length*sizeof(int));
		int *maxk_table = new int[seq_length];memset(maxk_table,-1,seq_length*sizeof(int));


		usize_t size_k = 0;
		for(usize_t k=min_k; k<=max_k; k++){
			dna5_t *cs_it_codes = new dna5_t[k];
			dolierlib::ssize_t *D = new dolierlib::ssize_t[k];

			size_k = 0;

			NSAIterator cs_iterator = NSAIterator::begin(f_ms, f_ds, k);
			while(cs_iterator.next()){
				cs_iterator.get_kmer(cs_it_codes);
				size_k++;


				if(mink_table[cs_iterator.i_start] == -1){
					mink_table[cs_iterator.i_start] = k;
					maxk_table[cs_iterator.i_start] = k;
				}
				else{
					if(k<mink_table[cs_iterator.i_start])	mink_table[cs_iterator.i_start] = k;
					if(k>maxk_table[cs_iterator.i_start])	maxk_table[cs_iterator.i_start] = k;
				}
			}

			std::cout<<k<<"-mers: "<<size_k<<"\n";

			delete [] cs_it_codes;
			delete [] D;
		}



		/**
		 * Strict SA
		 */
		//size
		//mapSA
		usize_t strictSA_size = 0;
		usize_t *mapSA = new usize_t[seq_length];

		for(usize_t i=0, s=0; i<seq_length; i++){
			if(mink_table[i] != -1){
				strictSA_size++;
				mapSA[i] = s;
				s++;
			}
		}

		std::cout<<"strict SA size = "<<strictSA_size<<"\n";


		/**
		 * DCount table
		 */

		usize_t **dcount = new usize_t*[seq_length];
		for(usize_t i=0; i<seq_length; i++){

			int x = mink_table[i]==-1 ? 0: maxk_table[i]-mink_table[i]+1;

			dcount[i] = new usize_t[ x *  (max_miss+1)];
			memset(dcount[i],0, (x *  (max_miss+1))*sizeof(usize_t));
		}



		std::cout<<"parallel searching...\n";
		start_t =start_time();



		int ratio = floor((double)strictSA_size /(double)nthreads);
		usize_t **t_ranges = new usize_t*[nthreads];
		for(int i=0; i<nthreads; i++)
			t_ranges[i] = new usize_t[2];

		if(nthreads<2){
			t_ranges[0][0] = 0;
			t_ranges[0][1] = seq_length-1;
		}
		else{
			t_ranges[0][0] = 0;
			t_ranges[nthreads-1][1] = seq_length-1;

			int t = 0;
			int count = 0;
			for(usize_t i=0; i< seq_length && t < nthreads-1; i++){
				if(mink_table[i] != -1){
					count++;
				}
				if(count == ratio){
					count = 0;
					t_ranges[t][1] = i;
					t_ranges[t+1][0] = i+1;
					t++;
				}
			}
		}

		std::cout<<"workload ratio "<<ratio<<", s_length "<<seq_length<<"\n";
		std::cout<<"ranges = ";
		for(int i=0; i<nthreads; i++){
			std::cout<<"["<<t_ranges[i][0]<<","<<t_ranges[i][1]<<"] ";
		}
		std::cout<<"\n";



		//double s_timers[nthreads][max_k+1][max_miss+1];
		double ***s_timers = new double**[nthreads];
		for(int i=0; i<nthreads; i++){
			s_timers[i] = new double*[max_k+1];
			for(int j=0; j<=max_k; j++){
				s_timers[i][j] = new double[max_miss+1];
				for(int k=0; k<=max_miss; k++){
					s_timers[i][j][k] = 0;
				}
			}
		}


		pthread_t *searcherPThreads = new pthread_t[nthreads];
		int rc;
		for(int i=0;i<nthreads;i++){
			SearcherThread* st = new SearcherThread(
					t_ranges[i][0],
					t_ranges[i][1],
					max_miss,
					verse, reverse,
					f_ms, f_ds,
					mink_table, maxk_table,
					dcount, s_timers[i]
					);
			rc = pthread_create(&searcherPThreads[i], NULL, SearcherThread::run, (void*)st);
			if(rc){
				printf("ERROR; return code from pthread_create() is %d\n", rc);
				exit(-1);
			}
		}

		for(int i=0;i<nthreads;i++){
			rc = pthread_join(searcherPThreads[i], NULL);
			if(rc){
				printf("ERROR; return code from pthread_join() is %d\n", rc);
				exit(-1);
			}
		}




		std::cout<<end_time(start_t)<<"\n";

		write(o_file, min_k, max_k, max_miss, verse, reverse, f_ms, f_ds,mink_table, maxk_table, dcount, strictSA_size);


		std::cout<<"done "<<end_time(total_t)<<"\n";


		for(int i=min_k; i<=max_k; i++){
			for(int j=0; j<=max_miss; j++){
				double km_timer = 0;
				for(int k=0; k<nthreads; k++)
					km_timer += s_timers[k][i][j];
				std::cout<<i<<"\t"<<j<<"\t"<<km_timer<<"\n";
			}
		}


		for(int i=0; i<nthreads; i++){
			for(int j=0; j<=max_k; j++){
				delete [] s_timers[i][j];
			}
			delete[] s_timers[i];
		}
		delete [] s_timers;

		//delete t_ranges
		for(int i=0; i<nthreads; i++)
			delete [] t_ranges[i];
		delete [] t_ranges;

		delete [] searcherPThreads;



		delete [] mink_table;
		delete [] maxk_table;
		delete [] mapSA;
		for(usize_t i=0; i<seq_length; i++)
			delete [] dcount[i];
		delete [] dcount;
	}





private:
class EnlargerThread{
public:
	const usize_t i_start;
	const usize_t i_end;
	const usize_t p_max_miss;
	const usize_t max_miss;
	const bool verse;
	const bool reverse;
	const DNA5MS_t &f_cs_gseq;
	const CsFullyContainer &f_container;
	Cs5DLIndex *f_cs_dl5index;
	const int p_maxk;
	const int *mink_table;
	const int *maxk_table;
	usize_t **dcount;
	double **timers;

	EnlargerThread(
			const usize_t _i_start,
			const usize_t _i_end,
			const usize_t _min_miss,
			const usize_t _max_miss,
			const bool _verse,
			const bool _reverse,
			const DNA5MS_t &_f_cs_gseq,
			const CsFullyContainer &_f_container,
			const int _mink,
			const int *_mink_table,
			const int *_maxk_table,
			usize_t **_dcount,
			double **_timers
			)
		: i_start(_i_start), i_end(_i_end),
		  p_max_miss(_min_miss), max_miss(_max_miss),
		  verse(_verse),reverse(_reverse),
		  f_cs_gseq(_f_cs_gseq), f_container(_f_container),
		  p_maxk(_mink), mink_table(_mink_table), maxk_table(_maxk_table),
		  dcount(_dcount), timers(_timers)
	{
#ifdef DEBUG_BEMINDEX_H_
			std::cout<<"BEMIndex::EnlargerThread::init>\n";
#endif
		f_cs_dl5index = new Cs5DLIndex (f_cs_gseq.seq, f_cs_gseq.seq_length, f_container.SA, f_cs_gseq.lengths, f_cs_gseq.nof_seqs);
#ifdef DEBUG_BEMINDEX_H_
			std::cout<<"BEMIndex::EnlargerThread::init<\n";
#endif
	}

	static void* run(void* argsptr){
#ifdef DEBUG_BEMINDEX_H_
			std::cout<<"BEMIndex::EnlargerThread::run>\n";
#endif
		EnlargerThread* st = (EnlargerThread*)argsptr;

		TIMEHANDLE total_t = start_time();
		TIMEHANDLE start_t;
		int count = 0;

		usize_t f_dcount=0,idoc_i=0;

		int mm = 0;

		for(usize_t i = st->i_start; i<=st->i_end; i++){
			if(st->maxk_table[i] > -1){

				count++;

				for(usize_t k = st->mink_table[i]; k <= st->maxk_table[i]; k++){

					dna5_t *cs_it_codes = new dna5_t[k];
					memcpy(cs_it_codes, &(st->f_cs_gseq.seq[  st->f_container.SA[i]]), k*sizeof(dna5_t));

					dna5_t *cs_it_codes_rev = new dna5_t[k];
					dolierlib::ssize_t *D = new dolierlib::ssize_t[k];


					k <= st->p_maxk ? mm = st->p_max_miss +1 : mm = 0;

					for(int m=mm; m<=st->max_miss; m++){

						start_t = start_time();
						bool* f_idocs = (bool*)calloc(st->f_cs_gseq.nof_seqs, sizeof(bool));



						if(st->verse)
							st->f_cs_dl5index->dl_count_wD(cs_it_codes, D, k, m, f_idocs);


						if(st->reverse){
							reverseComplement(cs_it_codes, cs_it_codes_rev, k);
							if(!st->verse || !equal_dna5(cs_it_codes, cs_it_codes_rev, k))
								st->f_cs_dl5index->dl_count_wD(cs_it_codes_rev, D, k, m, f_idocs);
						}


						f_dcount = 0;
						for(idoc_i=0; idoc_i<st->f_cs_gseq.nof_seqs; idoc_i++){
							if(f_idocs[idoc_i]) f_dcount++;
						}
#ifdef DEBUG_BEMINDEX_H_
						std::cout<<"dcount["<<i<<"]["<<(((k-st->mink_table[i])*(st->max_miss+1)) + m)<<"] = "<<f_dcount<<"\n";
#endif
						st->dcount[ i ][  ((k-st->mink_table[i])*(st->max_miss+1)) + m] = f_dcount;


						free(f_idocs);

						st->timers[k][m] += end_time(start_t);
					}

					delete [] cs_it_codes;
					delete [] cs_it_codes_rev;
					delete [] D;
				}

			}
		}

		std::cout<<"-t "<<end_time(total_t)<<"\t"<<count<<"\n";

		pthread_exit(NULL);
	}
};

public:


	void
	enlarge_pp(
			std::string	n_o_file,
			usize_t n_max_k,
			usize_t n_max_miss,
			int nthreads
			){
		TIMEHANDLE total_t = start_time();
		TIMEHANDLE start_t;

		start_t = start_time();

		DNA5MS_t f_ms(false);
			f_ms.seq = g_seq;
			f_ms.seq_length = g_length;
			f_ms.nof_seqs = nof_sequences;
			f_ms.lengths = seq_lengths;
		CsFullyContainer f_ds(true);
			build_fully_ds(f_ms, f_ds);
		Cs5DLIndex f_dlindex(f_ms.seq, f_ms.seq_length, f_ds.SA, f_ms.lengths, f_ms.nof_seqs);
		std::cout<<"cs build done: "<<end_time(start_t)<<"\n";
		usize_t seq_length = f_ms.seq_length;

#ifdef DEBUG_BEMINDEX_H_
		print_SA(f_ms.seq, f_ms.seq_length, f_ds.SA, 4);
#endif


//		usize_t min_k;
//			usize_t max_k;
//			usize_t max_miss;
//
//			bool verse;
//			bool reverse;
//
//			usize_t nof_sequences;
//			usize_t *seq_lengths;
//			usize_t g_length;
//			dna5_t* g_seq;

//			usize_t SA_length;
//			usize_t *SA;
//			_ui8 *SA_mink;
//			_ui8 *SA_maxk;
//			usize_t **dcounts;
//
//			bool freeOnExit;

		int *n_mink_table = new int[seq_length];memset(n_mink_table,-1,seq_length*sizeof(int));
		int *n_maxk_table = new int[seq_length];memset(n_maxk_table,-1,seq_length*sizeof(int));


		for(usize_t i = 0, j=0; i<SA_length; i++){
			while(f_ds.SA[j] != SA[i] && j<seq_length) j++;
			n_mink_table[j] = SA_mink[i];
			n_maxk_table[j] = SA_maxk[i];
		}



		usize_t size_k = 0;
		for(usize_t k=max_k+1; k<=n_max_k; k++){
			dna5_t *cs_it_codes = new dna5_t[k];
			dolierlib::ssize_t *D = new dolierlib::ssize_t[k];

			size_k = 0;

			NSAIterator it = NSAIterator::begin(f_ms, f_ds, k);
			while(it.next()){


//				std::cout<<it.i_start<<"\n";

				it.get_kmer(cs_it_codes);
				size_k++;

				if(n_mink_table[it.i_start] == -1){
					n_mink_table[it.i_start] = k;
					n_maxk_table[it.i_start] = k;
				}
				else{
					if(k>n_maxk_table[it.i_start])	n_maxk_table[it.i_start] = k;
				}
			}

			std::cout<<k<<"-mers: "<<size_k<<"\n";

			delete [] cs_it_codes;
			delete [] D;
		}

#ifdef DEBUG_BEMINDEX_H_
		std::cout<<"i\tmink\tmaxk\n";
		for(usize_t i = 0, j=0; i<f_ms.seq_length; i++){
			std::cout<<i<<"\t"<<n_mink_table[i]<<"\t"<<n_maxk_table[i]<<"\n";
		}
#endif


		usize_t n_strictSA_size = 0;
//		usize_t *mapSA = new usize_t[seq_length];

		for(usize_t i=0, s=0; i<seq_length; i++){
			if(n_mink_table[i] != -1){
				n_strictSA_size++;
//				mapSA[i] = s;
				s++;
			}
		}


		std::cout<<"old strict SA size = "<<SA_length<<"\n";
		std::cout<<"new strict SA size = "<<n_strictSA_size<<"\n";


		usize_t **n_dcount = new usize_t*[seq_length];
		for(usize_t i=0; i<seq_length; i++){
			int x = n_mink_table[i]==-1 ? 0: (n_maxk_table[i]-n_mink_table[i]+1) * (n_max_miss+1);
			n_dcount[i] = new usize_t[x];
			memset(n_dcount[i],0, (x)*sizeof(usize_t));
		}

		for(usize_t i = 0, j=0; i<SA_length; i++){
			while(f_ds.SA[j] != SA[i] && j<seq_length) j++;

			for(int k=n_mink_table[j]; k<=max_k;k++){
				for(int m=0; m<=max_miss; m++){
					n_dcount[j][ ((k-n_mink_table[j]) * (n_max_miss+1)) + m] = dcounts[i][ ((k-n_mink_table[j]) * (max_miss+1)) + m];
#ifdef DEBUG_BEMINDEX_H_
					std::cout<<"n_dcount["<<j<<"]["<<k<<"]["<<m<<"]..["<<(((k-n_mink_table[j]) * (n_max_miss+1)) + m)<<"] = dcount["<<i<<"][("<<k<<"-"<<n_mink_table[j]<<") * "<<((max_miss+1))<<"+m] = "<<dcounts[i][ ((k-n_mink_table[j]) * (max_miss+1)) + m]<<"\n";
#endif
				}
			}
		}


		std::cout<<"parallel searching...\n";
		start_t =start_time();



		int ratio = floor((double)n_strictSA_size /(double)nthreads);
		usize_t **t_ranges = new usize_t*[nthreads];
		for(int i=0; i<nthreads; i++)
			t_ranges[i] = new usize_t[2];

		if(nthreads<2){
			t_ranges[0][0] = 0;
			t_ranges[0][1] = seq_length-1;
		}
		else{
			t_ranges[0][0] = 0;
			t_ranges[nthreads-1][1] = seq_length-1;

			int t = 0;
			int count = 0;
			for(usize_t i=0; i< seq_length && t < nthreads-1; i++){
				if(n_mink_table[i] > -1){
					count++;
				}
				if(count == ratio){
					count = 0;
					t_ranges[t][1] = i;
					t_ranges[t+1][0] = i+1;
					t++;
				}
			}
		}

		std::cout<<"workload ratio "<<ratio<<", s_length "<<seq_length<<"\n";
		std::cout<<"ranges = ";
		for(int i=0; i<nthreads; i++){
			std::cout<<"["<<t_ranges[i][0]<<","<<t_ranges[i][1]<<"] ";
		}
		std::cout<<"\n";



		//double s_timers[nthreads][max_k+1][max_miss+1];
		double ***s_timers = new double**[nthreads];
		for(int i=0; i<nthreads; i++){
			s_timers[i] = new double*[n_max_k+1];
			for(int j=0; j<=n_max_k; j++){
				s_timers[i][j] = new double[n_max_miss+1];
				for(int k=0; k<=n_max_miss; k++){
					s_timers[i][j][k] = 0;
				}
			}
		}


		pthread_t *searcherPThreads = new pthread_t[nthreads];
		int rc;
		for(int i=0;i<nthreads;i++){
			EnlargerThread* st = new EnlargerThread(
					t_ranges[i][0],
					t_ranges[i][1],
					max_miss, n_max_miss,
					verse, reverse,
					f_ms, f_ds,
					max_k, n_mink_table, n_maxk_table,
					n_dcount, s_timers[i]
					);
			rc = pthread_create(&searcherPThreads[i], NULL, EnlargerThread::run, (void*)st);
			if(rc){
				printf("ERROR; return code from pthread_create() is %d\n", rc);
				exit(-1);
			}
		}
#ifdef DEBUG_BEMINDEX_H_
			std::cout<<"threads created\n";
#endif

		for(int i=0;i<nthreads;i++){
			rc = pthread_join(searcherPThreads[i], NULL);
			if(rc){
				printf("ERROR; return code from pthread_join() is %d\n", rc);
				exit(-1);
			}
		}
		std::cout<<end_time(start_t)<<"\n";

		write(n_o_file, min_k, n_max_k, n_max_miss, verse, reverse, f_ms, f_ds, n_mink_table, n_maxk_table, n_dcount, n_strictSA_size);


		std::cout<<"done "<<end_time(total_t)<<"\n";


		for(int i=min_k; i<=n_max_k; i++){
			for(int j=0; j<=n_max_miss; j++){
				double km_timer = 0;
				for(int k=0; k<nthreads; k++)
					km_timer += s_timers[k][i][j];
				std::cout<<i<<"\t"<<j<<"\t"<<km_timer<<"\n";
			}
		}


		for(int i=0; i<nthreads; i++){
			for(int j=0; j<=n_max_k; j++){
				delete [] s_timers[i][j];
			}
			delete[] s_timers[i];
		}
		delete [] s_timers;
//
		//delete t_ranges
		for(int i=0; i<nthreads; i++)
			delete [] t_ranges[i];
		delete [] t_ranges;

		delete [] searcherPThreads;



		delete [] n_mink_table;
		delete [] n_maxk_table;
		for(usize_t i=0; i<seq_length; i++)
			delete [] n_dcount[i];
		delete [] n_dcount;


	}




};

const usize_t BEMIndex::version = 0;

inline
usize_t nof(BEMIndex::iterator it){
	usize_t count = 0;
	while(it.next())
		count++;
	return count;
}


}


#endif /* BEMINDEX_H_ */
