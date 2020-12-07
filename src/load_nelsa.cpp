#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <new>


#include "timer.h"
#include "data_ts.h"

#include "FASTAReader.h"
#include "CstyleIndex.h"
#include "NSAIterator.h"
#include "SASearcher.h"


using namespace dolierlib;

int main(int argc, char* argv[]){
    std::cout<<"size max "<<std::numeric_limits<size_t>::max()<<"\n";
    std::cout<<"size max "<<std::numeric_limits<uint32_t>::max()<<"\n";

	TIMEHANDLE start_t = start_time();

    std::string ifasta = argv[1];
    std::string inelsa = argv[2];
    usize_t k = atoi(argv[3]);


	FASTAReader f_reader(ifasta);
    usize_t len = f_reader.get_huge_length();
	f_reader.close();
    std::cout<<"sequence length "<<len<<"\n";
    dna5_t *lseq = new dna5_t[ len ];
    f_reader.open(ifasta);
    bool readed = f_reader.read_huge(lseq, len);
    if(readed){
        std::cout<<"reading was ok\n";
    }
    f_reader.close();
    dna5_t *iseq = lseq;
    usize_t iseq_length = len;

    std::cout<<"done "<<end_time(start_t)<<"\n";




    
    usize_t length;
    usize_t *SA, *LCP, *NS;

    read_nelsa(inelsa, &length, SA, LCP, NS);
    std::cout<<"read_nelsa done "<<end_time(start_t)<<"\n";

    std::cout<<"ds length "<<length<<"\n";

    SASearcher searcher = SASearcher(iseq, SA, length);
    usize_t first_ncode = where_first_ncode(SA, NS, iseq_length);
    std::cout<<"first_ncode done "<<first_ncode<<"\n";
    std::cout<<"where_first_ncode done "<<end_time(start_t)<<"\n";



    char *it_syms = new char[k+1];
	it_syms[k] = '\0';
	dna5_t *kit_codes = new dna5_t[k];
    NSAIterator it = NSAIterator::begin(iseq,iseq_length,SA,LCP,NS,first_ncode,k);
    while(it.next()){
		it.get_kmer(kit_codes);
        print_dna5(kit_codes,k);std::cout<<" "<<it.multiplicity()<<" ["<<it.i_start<<","<<it.i_end<<"]\n";
    }

    std::cout<<"iteration done "<<end_time(start_t)<<"\n";

    for(usize_t i=0; i<100; i++){
		std::cout<<i<<"\t";
		for(usize_t j=0; j<20; j++){
			std::cout<<DNA5Alphabet::symbolFor(iseq[SA[i] + j]);
		}
		std::cout<<"\t"<<SA[i]<<"\t"<<LCP[i]<<"\t"<<NS[i]<<"\n";
	}


    delete [] SA;
    delete [] LCP;
    delete [] NS;
    delete [] iseq;

	exit(0);
}