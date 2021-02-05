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

#include "DNA5Alphabet.h"

using namespace dolierlib;

bool
is_unique_ext(
    NSAIterator &it,
    usize_t length, dna5_t *seq, 
    usize_t *SA, usize_t *invSA, usize_t *LCP, usize_t *NS,
    SASearcher &searcher
){
    //if(it.multiplicity() == 1) NOOOO!!!!! we need to check the k-1 suffix not the k-mer itself
    //    return true; 

    /*std::cout<<">";
    for(int i=0; i<it.k; i++){
        std::cout<<DNA5Alphabet::symbolFor(seq[SA[it.i_start]+i]);
    }
    std::cout<<" "<<it.i_start<<" "<<it.i_end<<"\n";*/
      
    bool ret = true;
    usize_t k = it.k-1;
    usize_t i = invSA[ SA[it.i_start]+1 ];

    usize_t s = i;
    while( (s>0) && (LCP[s]>=k) && (NS[s]>=k) ){
        s--;
    }

    usize_t e = i;
    while( (e<length-1) && (LCP[e+1]>=k) && (NS[e+1]>=k) ){
        e++;
    }

    /*std::cout<<s<<" - "<<e<<"\n";
    int si = s-1; if (si < 0) si = 0;
    int ei = e+1; if (ei > length-1) ei--;
    for(int i=si; i<=ei; i++){
        std::cout<<i<<" "<<SA[i]<<" ";
        for(int j=0; j<k+1; j++){
            std::cout<<DNA5Alphabet::symbolFor(seq[SA[i]+j]);
        }
        std::cout<<"\n";
    }*/

    for(i=s+1; i<=e; i++){
        if( (SA[i]+k<length) && (NS[i]>=k) && (SA[i-1]+k<length)  && (NS[i-1]>=k)){
            if( seq[SA[i]+k] != seq[SA[i-1]+k] ){
                //return false;
                ret = false;
                break;
            }
        }
    }
    /*if(ret)
        std::cout<<"TRUE\n";
    else
        std::cout<<"FALSE\n";
        */
    return ret;
}




bool 
is_uniquelly_extended(usize_t length, dna5_t *seq, usize_t *SA, usize_t *NS, usize_t k, usize_t start, usize_t end){
    if(end-start==1){
        return true;
    }
    for(usize_t i=start+1; i<end; i++){
        if( (SA[i]+k<length) && (NS[i]>=k) && (SA[i-1]+k<length)  && (NS[i-1]>=k)){
            if( seq[SA[i]+k-1] != seq[SA[i-1]+k-1] ){
                return false;
            }
        }
    }
    return true;
};


double
avg(std::list<usize_t> *l){
    double sum = 0.0;
    for(std::list<usize_t>::iterator it = l->begin(); it!=l->end(); it++){
        sum += ((double)(*it));
    }
    return ((double)sum) / ((double)l->size());
};


double
max(std::list<usize_t> *l){
    usize_t max = 0;
    for(std::list<usize_t>::iterator it = l->begin(); it!=l->end(); it++){
        if((*it) > max){
            max = *it;
        }
    }
    return max;
};

















usize_t
mhl(usize_t length, dna5_t *seq, usize_t *SA,usize_t *LCP, usize_t *NS, usize_t first_ncode){
    usize_t vmhl = 0;
    usize_t k = 1;
    while(vmhl == 0){
        //std::cout<<"searching for mhl at k "<<k<<"\n";
        NSAIterator it = NSAIterator::begin(seq,length,SA,LCP,NS,first_ncode,k);
        while(it.next()){
            if(it.multiplicity() == 1){
                vmhl = k;
                break;
            }
        }
        k++;
    }
    return vmhl;
};

usize_t
mrl(usize_t length, usize_t *LCP, usize_t *NS){
    usize_t m = 0;
    for(usize_t i=0; i<length; i++){
        if((LCP[i] > m) && (NS[i]>=LCP[i])){
            m = LCP[i];
        }
    }
    return m;
};




usize_t
mcl(usize_t length, dna5_t *seq, usize_t *SA,usize_t *LCP, usize_t *NS, usize_t first_ncode){
    
    usize_t m = 0;
    usize_t count;
    usize_t k = 1;

    while(m == 0){
        count = 0;
        NSAIterator it = NSAIterator::begin(seq,length,SA,LCP,NS,first_ncode,k);
        while(it.next()){
            count++;
        }

        if(count != pow(4,k)){
            m = k-1;
        }
        k++;
    }
    return m;

};

usize_t
mul(usize_t length, dna5_t *seq, 
    usize_t *SA, usize_t *invSA, usize_t *LCP, usize_t *NS 
    ,usize_t first_ncode
    , SASearcher &searcher
    ){

    usize_t k = 2;
	//dna5_t *kit_codes = new dna5_t[k];
    //usize_t *range = new usize_t[2];


    usize_t munique = 0;
    while(munique == 0){
        //delete [] kit_codes;

        //kit_codes = new dna5_t[k];

        //std::cout<<"searching for uniqueness at k = "<<k<<"\n";
        NSAIterator it = NSAIterator::begin(seq,length,SA,LCP,NS,first_ncode,k);
        while(it.next()){
		    /*it.get_kmer(kit_codes);
            if(searcher.get_range_iter(kit_codes+1, k-1,range)){
                if(is_uniquelly_extended(length, seq, SA, NS, k, range[0], range[1])){
                    munique = k;
                    break;
                }
            }*/
            if(is_unique_ext(it,length,seq,SA,invSA,LCP,NS, searcher)){
                munique = k;
                break;
            }
        }
        k++;
    }

    //delete [] kit_codes;
    //delete [] range;

    return munique;

};



int main(int argc, char* argv[]){
    std::cout<<"size max "<<std::numeric_limits<size_t>::max()<<"\n";
    std::cout<<"size max "<<std::numeric_limits<uint32_t>::max()<<"\n";

	TIMEHANDLE start_t = start_time();

    std::string ifasta = argv[1];
    usize_t min_k = atoi(argv[2]);
    usize_t max_k = atoi(argv[3]);

    std::string tseq;//string
	std::string iname;
	FASTAReader f_reader(ifasta);
	//f_reader.nextSeq(&iname, &tseq);
    usize_t len = f_reader.get_huge_length();
	f_reader.close();

    std::cout<<len<<"\n";

    dna5_t *lseq = new dna5_t[ len ];

    f_reader.open(ifasta);
    bool readed = f_reader.read_huge(lseq, len);
    if(readed){
        std::cout<<"reading was ok\n";
    }
    f_reader.close();

    dna5_t *iseq = lseq;
    usize_t iseq_length = len;

    std::cout<<"sequence length:"<<iseq_length<<"\n";

    usize_t *counts = new usize_t[256]; for(int i=0; i<256; i++)counts[i] = 0;
    usize_t sum = 0;

    for(usize_t i=0; i<iseq_length; i++){
        if(iseq[i]>=256){
            std::cout<<"OPS\n";
        }
        counts[ iseq[i] ]++;
    }
    for(usize_t i=0; i<256; i++){
        if(counts[i]>0)
            std::cout<<i<<" "<<counts[i]<<"\n";
        sum += counts[i];
    }
    std::cout<<"sum "<<sum<<"\n";
    std::cout<<"done "<<end_time(start_t)<<"\n";


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
    usize_t *invSA = new usize_t[iseq_length];
    for(usize_t i=0; i<iseq_length; i++){
        invSA[ SA[i] ] = i;
    }
    std::cout<<"invSA done "<<end_time(start_t)<<"\n";

    SASearcher searcher = SASearcher(iseq, SA, iseq_length);

    usize_t first_ncode = where_first_ncode(SA, NS, iseq_length);

    std::cout<<"lg_4(|C|) = "<<   log(iseq_length) / log(4) <<"\n";

    std::cout<<"mhl "<<mhl(iseq_length, iseq, SA, LCP, NS, first_ncode)<<"\n";
    std::cout<<"done "<<end_time(start_t)<<"\n";
    usize_t gmrl = mrl(iseq_length, LCP, NS);
    std::cout<<"mrl "<<mrl(iseq_length, LCP, NS)<<"\n";
    std::cout<<"done "<<end_time(start_t)<<"\n";
    std::cout<<"mcl "<<mcl(iseq_length, iseq, SA, LCP, NS, first_ncode)<<"\n";
    std::cout<<"done "<<end_time(start_t)<<"\n";
    usize_t k = 2;
    usize_t munique = mul(iseq_length, iseq, SA, invSA, LCP, NS, first_ncode, searcher);
    std::cout<<"minimal unique length "<<munique<<"\n";




    for(usize_t i=0; i<iseq_length; i++){
        if((LCP[i] == gmrl  && (NS[i]>=LCP[i]) ) ||(  i<iseq_length-1 && LCP[i+1] == gmrl && (NS[i+1]>=LCP[i])  )){
            std::cout<<i<<" "<<SA[i]<<" "<<LCP[i]<<" "<<NS[i]<<" : ";
            for(usize_t j = 0; j<gmrl; j++){
                std::cout << DNA5Alphabet::symbolFor( iseq[SA[i]+j] );
            }
            std::cout<<"|";
            for(usize_t j = 0; j<2; j++){
                std::cout << DNA5Alphabet::symbolFor( iseq[SA[i]+gmrl+j] );
            }
            std::cout<<"\n";
        }
    }


    //usize_t k = 2;
    k = 2;
    

    bool* starts = new bool[iseq_length];
    bool* covered = new bool[iseq_length];


    //for(k=min_k ; k<=max_k; k++){
    for(k=munique ; k<=gmrl+2; k++){
    //for(k=munique ; k<munique+1; k++){
        //std::cout<<"################################################################################\n";
        //std::cout<<k<<"\n";
        //std::cout<<"################################################################################\n";

        usize_t dk = 0;
        usize_t hk = 0;
        usize_t ucover = 0;
        usize_t nof_unique = 0;


        usize_t u_segments_nof = 0;
        usize_t u_segments_avg = 0;
        usize_t u_segments_max = 0;


        usize_t clength = 0;

        std::fill(starts, starts+iseq_length, false);
        std::fill(covered, covered+iseq_length, false);


        NSAIterator it = NSAIterator::begin(iseq,iseq_length,SA,LCP,NS,first_ncode,k);
        while(it.next()){
            dk++;
            if(is_unique_ext(it,iseq_length,iseq,SA,invSA,LCP,NS, searcher)){
                nof_unique++;

                for(usize_t posi=it.i_start; posi<it.i_end; posi++){
                    starts[SA[posi]] = true;
                    for(usize_t ki=0; ki<k; ki++){
                        covered[ SA[posi]+ki ] = true;
                    }
                }
            }
            if(it.multiplicity() == 1){
                hk++;
            }
        }
        for(usize_t i=0; i<iseq_length; i++){
            if(covered[i]){
                ucover++;
            }
        }

        usize_t p = 0;
        while(p<iseq_length){
            if(starts[p]){
                
                //NOOO!!! the k-mer after an expandable one is not expandale if the segment ends
                /*if((p<iseq_length-1) && !starts[p+1]){
                    std::cout<<"OPS "<<p<<" ";
                    for(int i=0; i<k; i++){
                        std::cout<<DNA5Alphabet::symbolFor(iseq[p+i]);
                    }
                    std::cout<<"\n";
                }*/

                //std::cout<<"> "<<p<<"\n";
                clength = 0;
                while(p<iseq_length && starts[p]){
                    p++;
                    clength++;
                }
                //clength += k-1;
                //the k-mer after an expandable one is not expandale if the segment ends
                clength += k;

                u_segments_nof++;
                u_segments_avg += clength;
                if(clength > u_segments_max){
                    u_segments_max = clength;
                }

                //std::cout<<p<<" "<<clength<<" "<<u_segments_nof<<" "<<u_segments_avg<<"\n";
            }
            p++;
        }



        std::cout<<k<<" "<<dk<<" "<<hk<<" "
            <<nof_unique<<" "<<((double)nof_unique)/((double)dk)
            <<" "
            <<ucover<<" "<<u_segments_nof<<" "<< (((double) u_segments_avg ) / ((double)u_segments_nof ))<<" "<<u_segments_max
            <<"\n";

    }


	std::cout<<"total time "<<end_time(start_t)<<"\n";


    delete [] iseq;
    delete [] SA;
    delete [] LCP;
    delete [] NS;
    delete [] covered;

	exit(0);
}