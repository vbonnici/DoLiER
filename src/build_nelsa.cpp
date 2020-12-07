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

using namespace dolierlib;

int main(int argc, char* argv[]){
    std::cout<<"size max "<<std::numeric_limits<size_t>::max()<<"\n";
    std::cout<<"size max "<<std::numeric_limits<uint32_t>::max()<<"\n";

	TIMEHANDLE start_t = start_time();

    std::string ifasta = argv[1];
    std::string ofile = argv[2];



    std::string tseq;//string
	std::string iname;
	FASTAReader f_reader(ifasta);
    usize_t len = f_reader.get_huge_length();
	f_reader.close();
    std::cout<<len<<"\n";
    dna5_t *lseq = new dna5_t[ len ];
    f_reader.open(ifasta);
    bool readed = f_reader.read_huge(lseq, len);
    if(readed){
        std::cout<<"reading the sequence was fine\n";
    }
    f_reader.close();

    dna5_t *iseq = lseq;
    usize_t iseq_length = len;

    std::cout<<"sequence length:"<<iseq_length<<"\n";

    usize_t* SA = build_SA<dna5_t>(iseq, iseq_length);
    std::cout<<"SA done "<<end_time(start_t)<<"\n";
    std::cout<<SA[0]<<"\n";
    std::cout<<SA[1]<<"\n";
    std::cout<<SA[iseq_length-2]<<"\n";
    std::cout<<SA[iseq_length-1]<<"\n";
	usize_t* LCP = build_LCP<dna5_t>(iseq, iseq_length, SA);
    std::cout<<"LCP done "<<end_time(start_t)<<"\n";
	usize_t* NS = build_NS(iseq, iseq_length, SA);
    std::cout<<"NS done "<<end_time(start_t)<<"\n";

    write_nelsa(ofile, iseq_length, SA, LCP, NS);
    std::cout<<"writing done "<<end_time(start_t)<<"\n";

    delete [] iseq;
    delete [] SA;
    delete [] LCP;
    delete [] NS;

	exit(0);
}