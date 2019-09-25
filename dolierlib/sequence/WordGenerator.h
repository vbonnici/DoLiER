/*
 * WordGenerator.h
 *
 *  Created on: Nov 11, 2013
 *      Author: vbonnici
 */


/*
 * the same of DNA4WordGenerator but can input/output words in std::string type
 */


#ifndef WORDGENERATOR_H_
#define WORDGENERATOR_H_


namespace dolierlib{


class WordGenerator{

	int k;
	int m;
	dna5_t *dkmer;
	dna5_t *ckmer;
	int *mi;
	int si;

	std::string s;

public:
	WordGenerator(std::string base, int m){
		this->m = m;
		k = base.length();
		dkmer = to_dna5(base);
		ckmer = new dna5_t[k];
		mi = new int[k];

		si = 0;
		mi[0] = 0;
		ckmer[0] = DNA5Alphabet::A;

		s = "";
	}

	WordGenerator(dna5_t *base, int k, int m){
		this->m = m;
		this-> k = k;
		dkmer =  new dna5_t[k];
		for(int i=0; i<k; i++)
			dkmer[i] = base[i];
		ckmer = new dna5_t[k];
		mi = new int[k];

		si = 0;
		mi[0] = 0;
		ckmer[0] = DNA5Alphabet::A;

		s = "";
	}


	~WordGenerator(){
		delete [] dkmer;
		delete [] ckmer;
		delete [] mi;
	}

	std::string
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
							s = "";
							for(int i=0; i<k; i++){
								s = s + DNA5Alphabet::symbolFor(ckmer[i] - 1);
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
