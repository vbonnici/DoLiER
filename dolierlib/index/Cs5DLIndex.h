/*
 * Cs57DLIndex.h
 *
 *  Created on: Oct 8, 2013
 *      Author: vbonnici
 */

#ifndef CS5DLINDEX_H_
#define CS5DLINDEX_H_



#include "data_ts.h"
#include "DNA5Alphabet.h"
#include "CstyleIndex.h"


//#define DEBUG_CS5DLINDEX_H_


namespace dolierlib{


/*
 * The class to solve the document listing problem. It's suffix array based.
 */


class Cs5DLIndex{
public:
	dna5_t	*g_seq;			//indexed concatenation sequence
	usize_t g_length;		//indexed concatenation sequence length == SA length
	usize_t *g_SA;			//suffix array of the concatenation sequence

	usize_t *ss_lengths; 	//length of original sequences
	usize_t nof_seqs;		//number of original sequences
	usize_t *ss_ids;		//ss_ids[i] = ID of the original sequence
	usize_t *cum_lengths; 	//cumulative lengths of original sequences

	usize_t last_valid;//we don't want to search on last positions of SA where seq[SA[i]]=$

	bool freeOnExit;		//if true then delete data on deconstruction

	Cs5DLIndex(
			dna5_t	*g_seq,
			usize_t g_length,
			usize_t *g_SA,
			usize_t *ss_lengths,
			usize_t nof_seqs
			){
		this->g_seq = g_seq;
		this->g_length = g_length;
		this->g_SA = g_SA;
		this->ss_lengths = ss_lengths;
		this->nof_seqs = nof_seqs;

		build_ss_ids();
		search_last_valid();
		build_cum_lengths();

		freeOnExit = true;
	}

	Cs5DLIndex(const Cs5DLIndex &b){
		g_seq = b.g_seq;
		g_length = b.g_length;
		g_SA = b.g_SA;
		ss_lengths = b.ss_lengths;
		nof_seqs = b.nof_seqs;
		ss_ids = b.ss_ids;
		cum_lengths = b.cum_lengths;
		last_valid = b.last_valid;

		freeOnExit = false;
	}


	void build_ss_ids(){
		ss_ids = new usize_t[g_length];
		usize_t *ids= new usize_t[g_length];
		size_t i=0;
		for(usize_t d=0; d<nof_seqs; d++){
			for(usize_t j=0; j<ss_lengths[d]; j++){
				ids[i] = d;
				i++;
			}
		}

		for(i=0; i<g_length; i++){
			ss_ids[i] = ids[g_SA[i]];
		}

		delete [] ids;
	}

	void search_last_valid(){
		dna5_t ecode = DNA5Alphabet::UNDEF;
		for(last_valid = g_length-1; last_valid>=0; last_valid--){
			if(g_seq[g_SA[last_valid]] < ecode)
				break;
		}
	}

	void build_cum_lengths(){
		cum_lengths = new usize_t[nof_seqs];
		cum_lengths[0] = 0;
		for(usize_t i=1; i<nof_seqs; i++){
			cum_lengths[i] = cum_lengths[i-1] + ss_lengths[i-1];
		}
	}


	~Cs5DLIndex(){
		if(freeOnExit){
			delete [] ss_ids;
			delete [] cum_lengths;
		}
	}

#define Cs5DLIndex_codeAt(pos, shift) g_seq[g_SA[pos] + shift]



	usize_t find_ii_end;
	dna5_t find_cs;
	usize_t find_m;

	/**
	 * min_i( seq[SA[i]+shift] = c ), for i in [start,end]
	 * return i_end+1 if not found
	 */
	inline
	usize_t find_first(dna5_t c, usize_t i_start, usize_t i_end, usize_t shift){
//		if(i_start == i_end){
//			if(Cs5DLIndex_codeAt(i_start, shift) != c)
//				return i_end+1;
//			return i_start;
//		}
//std::cout<<"fin_first("<<i_start<<","<<i_end<<")\n";
//		usize_t find_ii_end = i_end;
//		dna5_t find_cs = CsDLIndex_codeAt(i_start, shift);
		find_ii_end = i_end;
		find_cs = Cs5DLIndex_codeAt(i_start, shift);

		if(find_cs > c){
			return i_end+1;
		}
		find_cs = Cs5DLIndex_codeAt(i_end, shift);
		if(find_cs < c){
			return i_end+1;
		}

//		usize_t find_m;
		while(true){
			if(i_start == i_end){
				break;
			}
			find_m = ((i_end - i_start)/2) + i_start;

			find_cs = Cs5DLIndex_codeAt(find_m, shift);
			if(c <= find_cs){
				i_end = find_m;
			}
			else{
				i_start = find_m+1;
			}

		}
		if(c != Cs5DLIndex_codeAt(i_start, shift)){
			return find_ii_end+1;
		}

		return i_start;

	}
	/**
	 * max_i( seq[SA[i]+shift] = c ), for i in [start,end]
	 * return i_end+1 if not found
	 */
	inline
	usize_t find_last(dna5_t c, usize_t i_start, usize_t i_end, usize_t shift){
//		if(i_start == i_end){
//			if(Cs5DLIndex_codeAt(i_start, shift) != c)
//				return i_end+1;
//			return i_start;
//		}
//		usize_t find_ii_end = i_end;
//		dna5_t find_cs = CsDLIndex_codeAt(i_start, shift);
		find_ii_end = i_end;
		find_cs = Cs5DLIndex_codeAt(i_start, shift);
		if(find_cs > c){
			return i_end+1;
		}
		find_cs = Cs5DLIndex_codeAt(i_end, shift);
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

			find_cs = Cs5DLIndex_codeAt(find_m, shift);
			if(c < find_cs){
				i_end = find_m-1;
			}
			else{
				i_start = find_m;
			}

		}
		if(c != Cs5DLIndex_codeAt(i_start, shift))
			return find_ii_end+1;
		return i_start;
	}


	/*
	 * calculated the reverse array D of the BWA heuristic
	 */
	inline
	void calculate_D(const dna5_t *pattern, const usize_t length, int* D){
#ifdef DEBUG_CS5DLINDEX_H_
		std::cout<<"Cs5DLIndex::calculate_D()\n";
#endif
		usize_t i, j=0;
		int z = 0;
		usize_t i_start = 0, i_end = g_length-1;
		dna5_t c;
		usize_t l = length-1;

		for(i=0; i<length; i++){
			c = pattern[i];
			i_start = find_first(c, i_start, i_end, i-j);


			if(i_start > i_end){
				i_start = 0;
				i_end = g_length-1;
				j = i + 1;
				//j = i;
				z = z+1;
			}
			else{
				i_end = find_last(c, i_start, i_end, i-j);
			}
			//D[i] = z
			D[l-i] = z;
		}
#ifdef DEBUG_CS5DLINDEX_H_
		std::cout<<"Cs5DLIndex::calculate_D<\n";
#endif

	}

	/*
	 * docs must be initialized (docs =new boolean[nof_seqs]) and/or cleared before.
	 * After, docs[i] = true if pattern appears in the i-th sequence, with at most mismatches
	 * D must be initialized out. (D = new ssize_t[length]).
	 */
	void dl_count_wD(const dna5_t* pattern, ssize_t* D, const usize_t length, const ssize_t mismatches, bool *docs){
#ifdef DEBUG_CS5DLINDEX_H_
		std::cout<<"Cs5DLIndex::dl_count_wD(";print_dna5(pattern, length);std::cout<<","<<length<<","<<mismatches<<")\n";
#endif
		calculate_D(pattern, length, D);
		_dl_count_wD(pattern, length, mismatches, D, docs);
#ifdef DEBUG_CS5DLINDEX_H_
		std::cout<<"Cs5DLIndex::dl_count_wD<\n";
#endif
	}



	/*
	 * docs must be initialized before (docs =new boolean[nof_seqs])  and/or cleared before.
	 * After, docs[i] = true if pattern appears in the i-th sequence, with at most mismatches
	 * D must be initialized and calculated out of this method.
	 *
	 * It's performs a bounded binary search backtracking algorithm.
	 */
	void _dl_count_wD(const dna5_t *pattern, const usize_t length, const size_t mismatches, const int *D, bool *idocs){
#ifdef DEBUG_CS5DLINDEX_H_
		std::cout<<"Cs5DLIndex::_dl_count_wD()\n";
#endif
		dna5_t* codes = new dna5_t[length];
		usize_t **ranges = new usize_t*[length];
		for(usize_t i=0; i<length; i++)
			ranges[i] = new usize_t[2];
		int* miss = new int[length];

		ssize_t si = 0;
		ssize_t last_si = length - 1;
		usize_t t_start, t_end;
		usize_t idoci;

		dna5_t ecode = DNA5Alphabet::UNDEF +1;

		codes[0] = Cs5DLIndex_codeAt(0, 0);
		ranges[0][0] = 0;
		//ranges[0][1] = find_first(DNA5Alphabet::UNDEF, 0, g_length -1, si) -1;
		ranges[0][1] = last_valid;
		miss[0] = mismatches;

#ifdef DEBUG_CS5DLINDEX_H_
		std::cout<<"codes[0] = "<<DNA5Alphabet::symbolFor(codes[0])<<"\n";
		std::cout<<"ranges[0][0] = "<<ranges[0][0]<<"\n";
		std::cout<<"ranges[0][1] = "<<ranges[0][1]<<"\n";
		std::cout<<"miss[0] = "<<miss[0]<<"\n";
		std::cout<<"\n";
#endif

		while(si >= 0){
			codes[si]++;

#ifdef DEBUG_CS5DLINDEX_H_
			if(si==3 || si==4){
			std::cout<<"si = "<<si<<"\n";
			std::cout<<"code = "<<(int)codes[si]<<"\t"<<DNA5Alphabet::symbolFor(codes[si]-1)<<"\n";
			std::cout<<"miss = "<<miss[si]<<",  D[si] = "<<D[si]<<"\n";
			std::cout<<"ranges = ["<<ranges[si][0]<<"]["<<ranges[si][1]<<"]\n";
			}
#endif

			if(codes[si]>DNA5Alphabet::UNDEF ||  miss[si]+1<D[si]){
				si--;
			}
			else{

				if(si == last_si){
					if(miss[si] > 0){
						if(Cs5DLIndex_codeAt(ranges[si][0],si) != DNA5Alphabet::UNDEF){
							t_end = find_first(DNA5Alphabet::UNDEF, ranges[si][0], ranges[si][1], si) -1;
#ifdef DEBUG_CS5DLINDEX_H_
						std::cout<<"@found["<<ranges[si][0]<<","<<t_end<<"] miss["<<miss[si]<<"]\n";
#endif
							for(idoci=ranges[si][0]; idoci<=t_end; idoci++){
								idocs[ss_ids[idoci]] = true;

//								if(ss_ids[idoci] == 12){
//									print_dna5(pattern, length);std::cout<<"\n";
////									std::cout<<"@@@@@@@@@@@found["<<ranges[si][0]<<","<<t_end<<"] miss["<<miss[si]<<"]\n";
//								}
							}
						}
					}
					else{
						t_start = find_first(pattern[si], ranges[si][0], ranges[si][1], si);
						if(t_start != ranges[si][1]+1){
							t_end = find_last(pattern[si], t_start, ranges[si][1], si);

#ifdef DEBUG_CS5DLINDEX_H_
						std::cout<<"@found["<<t_start<<","<<t_end<<"] miss["<<miss[si]<<"]\n";
#endif

							for(idoci=t_start; idoci<=t_end; idoci++){
								idocs[ss_ids[idoci]] = true;

//								if(ss_ids[idoci] == 12){
//									print_dna5(pattern, length);std::cout<<"\n";
////									std::cout<<"###########found["<<ranges[si][0]<<","<<t_end<<"] miss["<<miss[si]<<"]\n";
//								}
							}

						}
					}
					si--;
				}
				else{
					if(miss[si] == 0){
						t_start = ranges[si][0];
						if(Cs5DLIndex_codeAt(ranges[si][0],si) != pattern[si]){
							t_start = find_first(pattern[si], ranges[si][0], ranges[si][1], si);
						}
						if(t_start != ranges[si][1]+1){
							t_end = find_last(pattern[si], t_start, ranges[si][1], si);
							si++;
							codes[si] = Cs5DLIndex_codeAt(t_start, si);
							ranges[si][0] = t_start;
							ranges[si][1] = t_end;
							miss[si] = 0;
							codes[si-1] = DNA5Alphabet::UNDEF;
						}
						else{
							si--;
						}
					}
					else{
//						std::cout<<"--->\n";
						codes[si] = Cs5DLIndex_codeAt(ranges[si][0],si)+1;
						if(codes[si]-1 < DNA5Alphabet::UNDEF){

							t_end = find_last(codes[si]-1, ranges[si][0], ranges[si][1], si);

							if(t_end <= ranges[si][1]){
								codes[si+1] = 0x0;
								ranges[si+1][0] = ranges[si][0];
								ranges[si+1][1] = t_end;
								if(codes[si]-1 != pattern[si]){
									miss[si+1] = miss[si] -1;
								}
								else{
									miss[si+1]  = miss[si];
								}

								if(t_end != ranges[si][1])
									ranges[si][0] = t_end+1;
								else
									codes[si] = ecode;

								si++;
							}
						}
						else{
							si--;
						}
					}

				}

//				std::cout<<"---\n";
			}

#ifdef DEBUG_CS5DLINDEX_H_
			std::cout<<"---\n";
#endif
		}



		delete [] codes;
		for(usize_t i=0; i<length; i++)
			delete [] ranges[i];
		delete []  ranges;
		delete [] miss;

#ifdef DEBUG_CS5DLINDEX_H_
		std::cout<<"Cs5DLIndex::_dl_count_wD<\n";
#endif

	}









//	void pwm_count(const dna5_t* pattern, ssize_t* D, const usize_t length, const ssize_t mismatches, double **matrix){
//		calculate_D(pattern, length, D);
//		_pwm_count(pattern, length, mismatches, D, matrix);
//	}
//
//
//
//	void _pwm_count(const dna5_t *pattern, const usize_t length, const size_t mismatches, const int *D, double **matrix){
//#ifdef DEBUG_CS5DLINDEX_H_
//		std::cout<<"Cs5DLIndex::_dl_count_wD()\n";
//#endif
//		dna5_t* codes = new dna5_t[length];
//		usize_t **ranges = new usize_t*[length];
//		for(usize_t i=0; i<length; i++)
//			ranges[i] = new usize_t[2];
//		int* miss = new int[length];
//
//		ssize_t si = 0;
//		ssize_t last_si = length - 1;
//		usize_t t_start, t_end;
//		usize_t idoci;
//
//		dna5_t ecode = DNA5Alphabet::UNDEF +1;
//
//		codes[0] = Cs5DLIndex_codeAt(0, 0);
//		ranges[0][0] = 0;
//		//ranges[0][1] = find_first(DNA5Alphabet::UNDEF, 0, g_length -1, si) -1;
//		ranges[0][1] = last_valid;
//		miss[0] = mismatches;
//
//#ifdef DEBUG_CS5DLINDEX_H_
//		std::cout<<"codes[0] = "<<DNA5Alphabet::symbolFor(codes[0])<<"\n";
//		std::cout<<"ranges[0][0] = "<<ranges[0][0]<<"\n";
//		std::cout<<"ranges[0][1] = "<<ranges[0][1]<<"\n";
//		std::cout<<"miss[0] = "<<miss[0]<<"\n";
//		std::cout<<"\n";
//#endif
//
//		while(si >= 0){
//			codes[si]++;
//
//#ifdef DEBUG_CS5DLINDEX_H_
//			std::cout<<"si = "<<si<<"\n";
//			std::cout<<"code = "<<(int)codes[si]<<"\t"<<DNA5Alphabet::symbolFor(codes[si]-1)<<"\n";
//			std::cout<<"miss = "<<miss[si]<<",  D[si] = "<<D[si]<<"\n";
//			std::cout<<"ranges = ["<<ranges[si][0]<<"]["<<ranges[si][1]<<"]\n";
//#endif
//
//			if(codes[si]>DNA5Alphabet::UNDEF ||  miss[si]+1<D[si]){
//				si--;
//			}
//			else{
//
//				if(si == last_si){
//					if(miss[si] > 0){
//						if(Cs5DLIndex_codeAt(ranges[si][0],si) != DNA5Alphabet::UNDEF){
//							t_end = find_first(DNA5Alphabet::UNDEF, ranges[si][0], ranges[si][1], si) -1;
//#ifdef DEBUG_CS5DLINDEX_H_
//						std::cout<<"@found["<<ranges[si][0]<<","<<t_end<<"] miss["<<miss[si]<<"]\n";
//#endif
//							for(idoci=ranges[si][0]; idoci<=t_end; idoci++){
////								print_dna5(g_seq,g_SA[idoci],length);std::cout<<"\n";
////								std::cout<<length<<"\n";
//								for(usize_t i=0; i<length; i++){
////									std::cout<<i<<".";
////									std::cout<<DNA5Alphabet::symbolFor(g_seq[g_SA[idoci]+i]);
//									if( g_seq[g_SA[idoci]+i] < DNA5Alphabet::N){
////										std::cout<<(int)g_seq[g_SA[idoci+i]];
//										matrix[i][  g_seq[g_SA[idoci]+i] ]++;
//									}
//								}
////								std::cout<<"\n";
//							}
//						}
//					}
//					else{
//						t_start = find_first(pattern[si], ranges[si][0], ranges[si][1], si);
//						if(t_start != ranges[si][1]+1){
//							t_end = find_last(pattern[si], t_start, ranges[si][1], si);
//
//#ifdef DEBUG_CS5DLINDEX_H_
//						std::cout<<"@found["<<t_start<<","<<t_end<<"] miss["<<miss[si]<<"]\n";
//#endif
//
//							for(idoci=t_start; idoci<=t_end; idoci++){
//
////								print_dna5(g_seq,g_SA[idoci],length);std::cout<<"\n";
////								std::cout<<length<<"\n";
//								for(usize_t i=0; i<length; i++){
////									std::cout<<DNA5Alphabet::symbolFor(g_seq[g_SA[idoci]+i]);
////									std::cout<<i<<".";
//									if( g_seq[g_SA[idoci]+i] < DNA5Alphabet::N){
//										matrix[i][  g_seq[g_SA[idoci]+i] ]++;
////										std::cout<<(int)g_seq[g_SA[idoci+i]];
//									}
//								}
////								std::cout<<"\n";
//							}
//
//						}
//					}
//					si--;
//				}
//				else{
//					if(miss[si] == 0){
//						t_start = ranges[si][0];
//						if(Cs5DLIndex_codeAt(ranges[si][0],si) != pattern[si]){
//							t_start = find_first(pattern[si], ranges[si][0], ranges[si][1], si);
//						}
//						if(t_start != ranges[si][1]+1){
//							t_end = find_last(pattern[si], t_start, ranges[si][1], si);
//							si++;
//							codes[si] = Cs5DLIndex_codeAt(t_start, si);
//							ranges[si][0] = t_start;
//							ranges[si][1] = t_end;
//							miss[si] = 0;
//							codes[si-1] = DNA5Alphabet::UNDEF;
//						}
//						else{
//							si--;
//						}
//					}
//					else{
////						codes[si] = Cs5DLIndex_codeAt(ranges[si][0],si)+1;
////						t_end = find_last(codes[si]-1, ranges[si][0], ranges[si][1], si);
////
////						if(t_end <= ranges[si][1]){
////							codes[si+1] = 0x0;
////							ranges[si+1][0] = ranges[si][0];
////							ranges[si+1][1] = t_end;
////							if(codes[si]-1 != pattern[si]){
////								miss[si+1] = miss[si] -1;
////							}
////							else{
////								miss[si+1]  = miss[si];
////							}
////
////							if(t_end != ranges[si][1])
////								ranges[si][0] = t_end+1;//!!!!! genius
////							else
////								codes[si] = ecode;
////
////							si++;
////						}
//
//
//						codes[si] = Cs5DLIndex_codeAt(ranges[si][0],si)+1;
//						if(codes[si]-1 < DNA5Alphabet::UNDEF){
//
//							t_end = find_last(codes[si]-1, ranges[si][0], ranges[si][1], si);
//
//							if(t_end <= ranges[si][1]){
//								codes[si+1] = 0x0;
//								ranges[si+1][0] = ranges[si][0];
//								ranges[si+1][1] = t_end;
//								if(codes[si]-1 != pattern[si]){
//									miss[si+1] = miss[si] -1;
//								}
//								else{
//									miss[si+1]  = miss[si];
//								}
//
//								if(t_end != ranges[si][1])
//									ranges[si][0] = t_end+1;
//								else
//									codes[si] = ecode;
//
//								si++;
//							}
//						}
//						else{
//							si--;
//						}
//					}
//
//				}
//
////				std::cout<<"---\n";
//			}
//
//#ifdef DEBUG_CS5DLINDEX_H_
//			std::cout<<"---\n";
//#endif
//		}
//
//
//
//		delete [] codes;
//		for(usize_t i=0; i<length; i++)
//			delete [] ranges[i];
//		delete []  ranges;
//		delete [] miss;
//
//#ifdef DEBUG_CS5DLINDEX_H_
//		std::cout<<"Cs5DLIndex::_dl_count_wD<\n";
//#endif
//
//	}

};

}


#endif /* CS5DLINDEX_H_ */
