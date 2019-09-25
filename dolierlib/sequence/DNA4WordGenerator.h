/*
 * WordGenerator.h
 *
 *  Created on: Nov 11, 2013
 *      Author: vbonnici
 */


/*
 * can be used to generate all the words at distance m (number of mismatches) of an input kmer.
 * It's a simple combinatorial algorithm with stop&restart.
 */

#ifndef DNA4WORDGENERATOR_H_
#define DNA4WORDGENERATOR_H_


namespace dolierlib{


class DNA4WordGenerator{

	int k;
	int m;
	dna5_t *dkmer;
	dna5_t *ckmer;
	int *mi;
	int si;

	dna5_t *s;

public:
	DNA4WordGenerator(dna5_t *base, int k, int m){
		this->m = m;
		this-> k = k;
		this->dkmer = base;
		ckmer = new dna5_t[k];
		mi = new int[k];

		si = 0;
		mi[0] = 0;
		ckmer[0] = DNA5Alphabet::A;

		s = new dna5_t[k];
	}


	~DNA4WordGenerator(){
		delete [] ckmer;
		delete [] mi;
		delete [] s;
	}

	dna5_t*
	word(){
		return s;
	}

	bool
	next(){
		if(si>=0){

			int miss = 0;
			while(si >= 0){
//				std::cout<<si<<"\t"<<((int)(ckmer[si] + 1))<<"\n";

				ckmer[si] = (ckmer[si] + 1);
				if(ckmer[si]>4){
					si--;
				}
				else{
					miss = (ckmer[si]-1 != dkmer[si]) ? 1 : 0;
					if(mi[si] + miss <= m){
						if(si == k-1){
							for(int i=0; i<k; i++){
								s[i] = ckmer[i]-1;
							}
							return true;
						}
						else{
							ckmer[si+1] = DNA5Alphabet::A;
							mi[si+1] = mi[si] + miss;
							si++;
						}
					}
				}
			}

			return false;

		}
		else
			return false;
	}

};


}

#endif /* WORDGENERATOR_H_ */
