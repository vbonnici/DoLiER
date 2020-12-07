/*
 * SASEARCHER_H_.h
 *
 *  Created on: Nov 12, 2013
 *      Author: vbonnici
 */

#ifndef SASEARCHER_H_
#define SASEARCHER_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data_ts.h"

#include "CstyleIndex.h"
#include "DNA5Alphabet.h"
#include "main_common.h"



namespace dolierlib{


/*
 * A class to perform kmer searching over a suffix array
 * ...without any additional data structure (LCP or document listing support).
 *
 *
 *
 */


class SASearcher{
public:
	dna5_t *seq;  //indexed sequence
	usize_t *SA;  //suffix array of the sequence
	usize_t length; //sequence length == SA lenth

	usize_t last_valid; //we don't want to search on last positions of SA where seq[SA[i]]=$

public:

	SASearcher(dna5_t *seq, usize_t *SA, usize_t length){
		this->seq = seq;
		this->SA = SA;
		this->length = length;
		search_last_valid();
	}

#define SASearcher_codeAt(pos, shift) seq[SA[pos] + shift]
//	usize_t find_ii_end;
//	dna5_t find_cs;
//	usize_t find_m;


	/*
	 * find min_i( seq[SA[i]+shift] = c ), for i in [start,end]
	 * return end+1 if not found
	 */
	inline
	usize_t find_first(dna5_t c, usize_t i_start, usize_t i_end, usize_t shift){
		usize_t find_ii_end = i_end;
		dna5_t find_cs = SASearcher_codeAt(i_start, shift);
		usize_t find_m;

		if(find_cs > c){
			return i_end+1;
		}
		find_cs = SASearcher_codeAt(i_end, shift);
		if(find_cs < c){
			return i_end+1;
		}

//		usize_t find_m;
		while(true){
			if(i_start == i_end){
				break;
			}
			find_m = ((i_end - i_start)/2) + i_start;

			find_cs = SASearcher_codeAt(find_m, shift);
			if(c <= find_cs){
				i_end = find_m;
			}
			else{
				i_start = find_m+1;
			}

		}
		if(c != SASearcher_codeAt(i_start, shift)){
			return find_ii_end+1;
		}

		return i_start;

	}


	/*
	 * find max_i( seq[SA[i]+shift] = c ), for i in [start,end]
	 * return end+1 if not found
	 */
	inline
	usize_t find_last(dna5_t c, usize_t i_start, usize_t i_end, usize_t shift){
		usize_t find_ii_end = i_end;
		dna5_t find_cs = SASearcher_codeAt(i_start, shift);
		usize_t find_m;

		if(find_cs > c){
			return i_end+1;
		}
		find_cs = SASearcher_codeAt(i_end, shift);
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

			find_cs = SASearcher_codeAt(find_m, shift);
			if(c < find_cs){
				i_end = find_m-1;
			}
			else{
				i_start = find_m;
			}

		}
		if(c != SASearcher_codeAt(i_start, shift))
			return find_ii_end+1;
		return i_start;
	}


	/*
	 * return number of occurrences of the input kmer
	 */
	usize_t
	count(dna5_t *kmer, int k){
		usize_t c_start = 0;
		usize_t c_end = length-1;
		usize_t t_start = 0;
		usize_t t_end = 0;
		usize_t last_si = k-1;
		usize_t kk = static_cast<usize_t>(k);

		//it's a simple binary search over SA
		for(usize_t si=0; si<kk; si++){
			t_start = find_first(kmer[si], c_start, c_end, si);
			if(t_start != c_end+1){
				if(si != last_si){
					t_end = find_last(kmer[si], t_start, c_end, si);
					c_start = t_start;
					c_end = t_end;
				}
				else{
					t_end = find_last(kmer[si], t_start, c_end, si);
				}
			}
			else{
				return 0;
			}
		}

		return t_end - t_start +1;
	}

	/*
	 * return true if the sequence contains this kmer
	 */
	usize_t
	contains(dna5_t *kmer, int k){
		usize_t c_start = 0;
		usize_t c_end = length-1;
		usize_t t_start = 0;
		usize_t t_end = 0;
		usize_t last_si = k-1;
		usize_t kk = static_cast<usize_t>(k);

		//it's a simple binary search over SA
		for(usize_t si=0; si<kk; si++){
			t_start = find_first(kmer[si], c_start, c_end, si);
			if(t_start != c_end+1){
				if(si != last_si){
					t_end = find_last(kmer[si], t_start, c_end, si);
					c_start = t_start;
					c_end = t_end;
				}
			}
			else{
				return length;
			}
		}

		return t_start;
	}


	/*
	 * return true if the sequence contains this kmer. 
	 * The SA ranges are written in the "range" parameter, 
	 * such that range[0] = first index of the kmer
	 * and range[1] = last index + 1 of the kmer
	 */
	bool
	get_range(dna5_t *kmer, int k,usize_t *range){
		usize_t c_start = 0;
		usize_t c_end = length-1;
		usize_t t_start = 0;
		usize_t t_end = 0;
		usize_t last_si = k-1;
		usize_t kk = static_cast<usize_t>(k);

		//it's a simple binary search over SA
		for(usize_t si=0; si<kk; si++){
			t_start = find_first(kmer[si], c_start, c_end, si);
			if(t_start != c_end+1){
				//if(si != last_si){
					t_end = find_last(kmer[si], t_start, c_end, si);
					c_start = t_start;
					c_end = t_end;
				//}
			}
			else{
				return false;
			}
		}

		range[0] = c_start;
		range[1] = c_end+1;

		return true;
	}


/*
        public IELSAIterator find(B3Nucleotide[] q) {
            if(q.length<1)
                return null;
            
            int l = 0;
            int r = r_limit;
            int ll,lr, rl, rr;
            int m;
            int cq;

            for(int i=0; i<q.length; i++){
                    cq = cc(q,i);

                    if(cq < cc(sa[l] +i))
                            return null;
                    if(cq > cc(sa[r] +i))
                            return null;


                    ll = l;
                    lr = r;
                    while(true){
                            if(ll==lr)
                                    break;
                            m = ll + ((lr-ll)/2);
                            if(cq <= cc(sa[m] +i)){
                                    lr = m;
                            }
                            else if(cq > cc(sa[m] +i)){
                                    ll = m+1;
                            }
                    }
                    l = ll;
                    if(cq != cc(sa[l] +i))
                            return null;


                    rl = l;
                    rr = r;
                    while(true){
                            if(rl==rr)
                                    break;
                            m = rl + ((rr-rl)/2);
                            if((rr-rl)%2 != 0)
                                    m++;
                            if(cq < cc(sa[m] +i)){
                                    rr = m-1;
                            }
                            else if(cq >= cc(sa[m] +i)){
                                    rl = m;
                            }
                    }
                    r = rr;
                    if(cq != cc(sa[r] +i))
                            return null;

            }
            return new _Iterator(this, q.length, l, r+1);
        }
*/
	bool
	get_range_iter(dna5_t *q, usize_t k, usize_t *range){

			usize_t l = 0;
            usize_t r = length-1;
            usize_t ll,lr, rl, rr;
            usize_t m;
            dna5_t cq;

            for(usize_t i=0; i<k; i++){
                    cq = q[i];

					while(SA[l]+i >= length){
						l++;
					}
					while((r>=l) && (SA[r]+i>=length)){
						r--;
					}

                    if(cq < seq[SA[l]+i])
                            return false;
                    if(cq > seq[SA[r]+i])
                            return false;


                    ll = l;
                    lr = r;
                    while(true){
                            if(ll==lr)
                                    break;
                            m = ll + ((lr-ll)/2);
                            if(cq <= seq[SA[m] +i]){
                                    lr = m;
                            }
                            else if(cq > seq[SA[m] +i]){
                                    ll = m+1;
                            }
                    }
                    l = ll;
                    if(cq != seq[SA[l] +i])
                            return false;

                    rl = l;
                    rr = r;
                    while(true){
                            if(rl==rr)
                                    break;
                            m = rl + ((rr-rl)/2);
                            if((rr-rl)%2 != 0)
                                    m++;
                            if(cq < seq[SA[m] +i]){
                                    rr = m-1;
                            }
                            else if(cq >= seq[SA[m] +i]){
                                    rl = m;
                            }
                    }
                    r = rr;
                    if(cq != seq[SA[r] +i])
                            return false;

            }
			range[0] = l;
			range[1] = r;
			return true;
	}

	/*
	 * return true if the sequence contains the reverse complement of this kmer
	 */
	usize_t
	containsRC(dna5_t *kmer, int k){
		usize_t c_start = 0;
		usize_t c_end = length-1;
		usize_t t_start = 0;
		usize_t t_end = 0;
		usize_t last_si = k-1;
		usize_t kk = static_cast<usize_t>(k);

		//it's a simple binary search over SA
		for(usize_t si=kk-1; si>=0; si++){
			t_start = find_first(DNA5Alphabet::complement(kmer[si]), c_start, c_end, si);
			if(t_start != c_end+1){
				if(si != last_si){
					t_end = find_last(DNA5Alphabet::complement(kmer[si]), t_start, c_end, si);
					c_start = t_start;
					c_end = t_end;
				}
			}
			else{
				return length;
			}
		}

		return t_start;
	}

















	/*
	 * max_i( seq[SA[i]] < $)
	 */
	void search_last_valid(){
		dna5_t ecode = DNA5Alphabet::UNDEF;
		for(last_valid = length-1; last_valid>=0; last_valid--){
			if(seq[SA[last_valid]] < ecode)
				break;
		}
	}








	/*
	 * calculated the reverse array D of the BWA heuristic
	 */
	inline
	void calculate_D(const dna5_t *pattern, const usize_t length, int* D){
		usize_t i, j=0;
		int z = 0;
		usize_t i_start = 0, i_end = this->length-1;
		dna5_t c;
		usize_t l = length-1;

		for(i=0; i<length; i++){
			c = pattern[i];
			i_start = find_first(c, i_start, i_end, i-j);

			if(i_start > i_end){
				i_start = 0;
				i_end = this->length-1;
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
	}

//	usize_t
//	count(const dna5_t *kmer, const usize_t k, const usize_t m){
//
//	}

	/*
	 * return the number of occurrence of the kmer (pattern) allowing for mismatches.
	 * D must be initialized out. (D = new ssize_t[length])
	 */
	usize_t
	count(const dna5_t* pattern, ssize_t* D, const usize_t length, ssize_t mismatches){
		calculate_D(pattern, length, D);
		return count_wD(pattern, length, mismatches, D);
	}

	/*
	 * return the number of occurrence of the kmer (pattern) allowing for mismatches.
	 * D must be initialized and calculated out of this method.
	 *
	 * It's performs a bounded binary search backtracking algorithm.
	 */
	usize_t
	count_wD(const dna5_t *pattern, const usize_t length, const size_t mismatches, const ssize_t *D){


		usize_t ret = 0;


		dna5_t* codes = new dna5_t[length];
		usize_t **ranges = new usize_t*[length];
		for(usize_t i=0; i<length; i++)
			ranges[i] = new usize_t[2];
		ssize_t* miss = new ssize_t[length];

		ssize_t si = 0;
		ssize_t last_si = length - 1;
		usize_t t_start, t_end;

		dna5_t ecode = DNA5Alphabet::UNDEF +1;

		codes[0] = SASearcher_codeAt(0, 0);
		ranges[0][0] = 0;
		//ranges[0][1] = find_first(DNA5Alphabet::UNDEF, 0, g_length -1, si) -1;
		ranges[0][1] = last_valid;
		miss[0] = mismatches;

		while(si >= 0){

			codes[si]++;

			if(codes[si]>DNA5Alphabet::UNDEF ||  miss[si]+1<D[si]){
				si--;
			}
			else{

				if(si == last_si){
					if(miss[si] > 0){
						if(SASearcher_codeAt(ranges[si][0],si) != DNA5Alphabet::UNDEF){
							t_end = find_first(DNA5Alphabet::UNDEF, ranges[si][0], ranges[si][1], si) -1;
							ret += t_end - ranges[si][0] + 1;
//							for(idoci=ranges[si][0]; idoci<=t_end; idoci++){
//								idocs[ss_ids[idoci]] = true;
//							}
						}
					}
					else{
						t_start = find_first(pattern[si], ranges[si][0], ranges[si][1], si);
						if(t_start != ranges[si][1]+1){
							t_end = find_last(pattern[si], t_start, ranges[si][1], si);
							ret += t_end - t_start +1;
//							for(idoci=t_start; idoci<=t_end; idoci++){
//								idocs[ss_ids[idoci]] = true;
//
//							}
						}
					}
					si--;
				}
				else{
					if(miss[si] == 0){
						t_start = ranges[si][0];
						if(SASearcher_codeAt(ranges[si][0],si) != pattern[si]){
							t_start = find_first(pattern[si], ranges[si][0], ranges[si][1], si);
						}
						if(t_start != ranges[si][1]+1){
							t_end = find_last(pattern[si], t_start, ranges[si][1], si);
							si++;
							codes[si] = SASearcher_codeAt(t_start, si);
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
						codes[si] = SASearcher_codeAt(ranges[si][0],si)+1;
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
			}
		}

		delete [] codes;
		for(usize_t i=0; i<length; i++)
			delete [] ranges[i];
		delete []  ranges;
		delete [] miss;

		return ret;
	}

};


}


#endif /* SASEARCHER_H_ */
