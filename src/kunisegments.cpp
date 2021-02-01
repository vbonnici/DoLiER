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




bool
is_unique_ext(
    NSAIterator &it,
    usize_t length, dna5_t *seq, 
    usize_t *SA, usize_t *invSA, usize_t *LCP, usize_t *NS,
    SASearcher &searcher
){
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

    //k++;
    //if(e-s == 1){
    //    return true;
    //}
    for(i=s+1; i<e; i++){
        if( (SA[i]+k<length) && (NS[i]>=k) && (SA[i-1]+k<length)  && (NS[i-1]>=k)){
            if( seq[SA[i]+k] != seq[SA[i-1]+k] ){
                //return false;
                ret = false;
                break;
            }
        }
    }

/*
    k = it.k;
    dna5_t *codes = new dna5_t[k];
    it.get_kmer(codes);
    usize_t *range = new usize_t[2];
    searcher.get_range_iter(codes+1, k-1,range);
    if(s!=range[0] ||  e!=range[1]){
        for(usize_t j=0; j<k; j++) std::cout<<codes[j]<<"|";

        std::cout<<" : "<<k<<" ";
        print_dna5(codes, k);
        std::cout<<" - ";
        print_dna5(codes+1, k-1);
        std::cout<<": ["<<s<<","<<e<<"] vs ["<<range[0]<<","<<range[1]<<"]\n";
    }
    delete [] codes;
    delete [] range;
*/
    //return true;
    return ret;
}




bool 
is_uniquelly_extended(usize_t length, dna5_t *seq, usize_t *SA, usize_t *NS, usize_t k, usize_t start, usize_t end){
    //std::cout<<"univocal "<<k<<" "<<start<<" "<<end<<"\n";

    if(end-start==1){
        return true;
    }
    for(usize_t i=start+1; i<end; i++){
        //std::cout<<"i:"<<i<<" "<<(SA[i]+k)<<" "<<(NS[i])<<" "<<(SA[i-1]) <<" "<<(NS[i-1])<<"\n";
        //std::cout<<(SA[i]+k<length)<<" "<<(NS[i]>=k)<<" "<<(SA[i-1]+k<length) <<" "<<(NS[i-1]>=k)<<"\n";
        if( (SA[i]+k<length) && (NS[i]>=k) && (SA[i-1]+k<length)  && (NS[i-1]>=k)){
            //std::cout<<k<<" "<<i<<" "<<SA[i-1]<<" "<<SA[i]<<" - "<<(seq[SA[i-1]+k-1])<<" "<<(seq[SA[i]+k-1])<<"\n";
            //print_dna5(seq, SA[i-1]-k, (2*k)+2);std::cout<<"\n";
            //print_dna5(seq, SA[i]-k, (2*k)+2);std::cout<<"\n";
            //std::cout<<k<<" "<<i<<" "<<DNA5Alphabet::symbolFor(seq[SA[i-1]+k])<<" "<<DNA5Alphabet::symbolFor(seq[SA[i]+k])<<"\n";

            if( seq[SA[i]+k-1] != seq[SA[i-1]+k-1] ){
                //std::cout<<"false\n";
                return false;
            }
        }
    }
    return true;
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


/*
void
print_uniques(
    dna5_t *iseq, usize_t iseq_length,
    SASearcher &searcher,
    usize_t *SA, usize_t *invSA, usize_t *LCP, usize_t *NS, 
    usize_t first_ncode, usize_t k
){
    char *it_syms = new char[k+1];
	it_syms[k] = '\0';
	dna5_t *kit_codes = new dna5_t[k];
    dna5_t *sit_codes = new dna5_t[k-1];
    usize_t *range = new usize_t[2];


    NSAIterator it = NSAIterator::begin(iseq,iseq_length,SA,LCP,NS,first_ncode,k);
    while(it.next()){
		//it.get_kmer(kit_codes);

        //if(searcher.get_range_iter(kit_codes+1, k-1,range)){
        //    if(is_uniquelly_extended(iseq_length, iseq, SA, NS, k, range[0], range[1])){
        if(is_unique_ext(it,iseq_length,iseq,SA,invSA,LCP,NS)){
                print_dna5(kit_codes, k);
                std::cout<<" "<<it.multiplicity()<<" ";
                print_dna5(kit_codes+1, k-1);
                std::cout<<" "<<(range[1] - range[0])<<" ["<<range[0]<<","<<range[1]<<"]";
                std::cout<<"\n";
            //}
        }
    }
}
*/





int main(int argc, char* argv[]){
    std::cout<<"size max "<<std::numeric_limits<size_t>::max()<<"\n";
    std::cout<<"size max "<<std::numeric_limits<uint32_t>::max()<<"\n";


    /*
    size_t ll = 3088269800;
    //dna5_t *llseq = new (std::nothrow) dna5_t[ ll ];
    size_t *llseq = new (std::nothrow) size_t[ ll ];

    if(!llseq){
        std::cout<<"no new\n";
    }else{
        std::cout<<"ok new\n";
    }

    for(usize_t i = 0; i<ll; i++){
        llseq[i] = i;
    }
    std::cout<<llseq[100]<<"\n";
    std::cout<<llseq[ll-1]<<"\n";

    exit(0);
    */

	TIMEHANDLE start_t = start_time();

    std::string ifasta = argv[1];
    usize_t min_k = atoi(argv[2]);
    usize_t max_k = atoi(argv[3]);
    std::string ianno = argv[4];

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

    /*
    std::cout<<"sequence length:"<<tseq.length()<<"\n";
    std::cout<<"reading done "<<end_time(start_t)<<"\n";

    usize_t iseq_length = tseq.length();

    dna5_t *iseq = new dna5_t[tseq.length()];
    for(usize_t i=0; i<iseq_length; i++){
        iseq[i] =  DNA5Alphabet::codeFor(tseq[i]);
    }
    //std::cout<<tseq<<"\n";
    tseq.clear();
    */

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

    bool *annos = new bool[iseq_length];
    for(usize_t i=0; i<iseq_length; i++){
        annos[i] = false;
    }

    std::ifstream file(ianno);
    std::string   line;
    std::string v;
    usize_t  val1;
    usize_t  val2;
    while(std::getline(file, line)){
        std::stringstream   linestream(line);
        std::string         data;
        std::getline(linestream, data, '\t');
        linestream >> v >> v >> val1 >> val2;
        for(usize_t i=val1; i<=val2; i++)
            annos[i] = true;

        //std::cout<<line<<"\n";
        //std::cout<<val1<<" "<<val2<<"\n";
    }

    usize_t annot = 0;
    for(usize_t i=0; i<iseq_length; i++){
       if( annos[i] == true ){
           annot++;
       }
    }

    std::cout<<"annotation cover "<<annot<<"\n";
    std::cout<<"annotations done "<<end_time(start_t)<<"\n";


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



    /*NSAIterator tit = NSAIterator::begin(iseq,iseq_length,SA,LCP,NS,first_ncode,max_k);
    while(tit.next()){
        if(is_unique_ext(tit,iseq_length,iseq,SA,invSA,LCP,NS, searcher)){
        }
    }
    exit(0);*/





    std::cout<<"lg_4(|C|) = "<<   log(iseq_length) / log(4) <<"\n";
    std::cout<<"mhl "<<mhl(iseq_length, iseq, SA, LCP, NS, first_ncode)<<"\n";
    std::cout<<"done "<<end_time(start_t)<<"\n";
    std::cout<<"mrl "<<mrl(iseq_length, LCP, NS)<<"\n";
    std::cout<<"done "<<end_time(start_t)<<"\n";
    std::cout<<"mcl "<<mcl(iseq_length, iseq, SA, LCP, NS, first_ncode)<<"\n";
    std::cout<<"done "<<end_time(start_t)<<"\n";
    //print_SA_LCP_N(iseq, iseq_length, SA, LCP, NS, 6);
    
    //std::cout<<"first N code at "<<first_ncode<<"\n";




    usize_t k = 2;
    
    //char *it_syms = new char[k+1];
	//it_syms[k] = '\0';
	//dna5_t *kit_codes = new dna5_t[k];
    //usize_t *range = new usize_t[2];


    /*

    NSAIterator it = NSAIterator::begin(iseq,iseq_length,SA,LCP,NS,first_ncode,k);
    while(it.next()){
		it.get_kmer(kit_codes);
        DNA5Alphabet::toSymbolArray(kit_codes, it_syms, k);
        std::cout<<it_syms<<" "<<it.i_start<<" "<<it.i_end<<"\n";

        for(int i=1; i<k; i++) sit_codes[i-1] = kit_codes[i];
        print_dna5(sit_codes,k-1);std::cout<<":";

        if(searcher.get_range(sit_codes, k-1,range)){
            std::cout<<range[0]<<" "<<range[1]<<"\nis it? "<<
            is_uniquelly_extended(iseq_length, iseq, SA, NS, k, range[0], range[1]) <<"\n";
        }
        
        std::cout<<"---\n";

    }

    exit(0);
    */


    
    //usize_t munique = mul(iseq_length, iseq, SA, LCP, NS, first_ncode, searcher);
    usize_t munique = mul(iseq_length, iseq, SA, invSA, LCP, NS, first_ncode, searcher);
    std::cout<<"minimal unique length "<<munique<<"\n";
    //exit(0);


    //print_uniques(iseq,iseq_length,searcher,SA, LCP, NS, first_ncode,munique);
    //exit(0);


    //munique = 9

    bool* covered = new bool[iseq_length];

    usize_t mm = munique;
    if(min_k != 0){
        mm = min_k;
    }

    //for(k=munique ; k<=max_k; k++){
    for(k=mm ; k<=max_k; k++){
        //delete [] it_syms;
        //delete [] kit_codes;
        //it_syms = new char[k+1]; it_syms[k] = '\0';
        //kit_codes = new dna5_t[k];


        usize_t dk = 0;
        usize_t ucover = 0;
        usize_t nof_unique = 0;
        usize_t ncover = 0;
        usize_t nof_nonunique = 0;

        usize_t a_ucover = 0;


        //std::list<usize_t> u_segments;  
        usize_t u_segments_nof = 0;
        usize_t u_segments_avg = 0;
        usize_t u_segments_max = 0;
        //std::list<usize_t> n_segments;  
        usize_t n_segments_nof = 0;
        usize_t n_segments_avg = 0;
        usize_t n_segments_max = 0;


        usize_t clength = 1;

        std::fill(covered, covered+iseq_length, false);

        NSAIterator it = NSAIterator::begin(iseq,iseq_length,SA,LCP,NS,first_ncode,k);
        while(it.next()){
            dk++;
            //it.get_kmer(kit_codes);
            //if(searcher.get_range_iter(kit_codes+1, k-1,range)){
                //if(is_uniquelly_extended(iseq_length, iseq, SA, NS, k, range[0], range[1])){
                if(is_unique_ext(it,iseq_length,iseq,SA,invSA,LCP,NS, searcher)){
                    nof_unique++;

                    for(usize_t posi=it.i_start; posi<it.i_end; posi++){
                        for(usize_t ki=0; ki<k; ki++){
                            covered[ SA[posi]+ki ] = true;
                        }
                    }
                }
            //}
        }

        for(usize_t i=0; i<iseq_length; i++){
            if(covered[i]){
                ucover++;
                if(annos[i]){
                    a_ucover++;
                }
            }
        }

        for(usize_t p=1; p<iseq_length; p++){
            if(covered[p] != covered[p-1]){
                if(!covered[p]){
                    //u_segments.push_back(clength);
                    u_segments_nof++;
                    u_segments_avg += clength;
                    if(clength > u_segments_max){
                        u_segments_max = clength;
                    }
                }
                clength = 1;
          
            }
            else{
                clength++;
            }
        }
        if(covered[iseq_length-1] == true){
            //u_segments.push_back(clength);
            u_segments_nof++;
            u_segments_avg += clength;
            if(clength > u_segments_max){
                u_segments_max = clength;
            }
        }



        std::fill(covered, covered+iseq_length, false);


        it = NSAIterator::begin(iseq,iseq_length,SA,LCP,NS,first_ncode,k);
        while(it.next()){
            //it.get_kmer(kit_codes);
            //if(searcher.get_range_iter(kit_codes+1, k-1,range)){
                //if(!is_uniquelly_extended(iseq_length, iseq, SA, NS, k, range[0], range[1])){
                if(is_unique_ext(it,iseq_length,iseq,SA,invSA,LCP,NS, searcher)){
                    nof_nonunique++;

                    for(usize_t posi=it.i_start; posi<it.i_end; posi++){
                        for(usize_t ki=0; ki<k; ki++){
                            covered[ SA[posi]+ki ] = true;
                        }
                    }
                }
            //}
        }

        for(usize_t i=0; i<iseq_length; i++){
            if(covered[i]) ncover++;
        }

        clength = 1;

        for(usize_t p=1; p<iseq_length; p++){
            if(covered[p] != covered[p-1]){
                if(!covered[p]){
                    //n_segments.push_back(clength);
                    n_segments_nof++;
                    n_segments_avg += clength;
                    if(clength > n_segments_max){
                        n_segments_max = clength;
                    }
                }
                clength = 1;
            }
            else{
                clength++;
            }
        }
        if(covered[iseq_length-1] == true){
            //n_segments.push_back(clength);
            n_segments_nof++;
            n_segments_avg += clength;
            if(clength > n_segments_max){
                n_segments_max = clength;
            }
        }


        std::cout<<k<<" "<<dk<<" "
            <<nof_unique<<" "<<((double)nof_unique)/((double)dk)<<" "
            <<ucover<<" "<<u_segments_nof++<<" "<< (((double) u_segments_avg ) / ((double)u_segments_nof ))<<" "<<u_segments_max<<" "
            <<" | "
            <<nof_nonunique<<" "<<ncover<<" "<<n_segments_nof<<" "<< (((double) n_segments_avg ) / ((double)n_segments_nof ))<<" "<<n_segments_max<<" "
            <<" | "
            <<a_ucover<<" "<< ((double)a_ucover)/((double)annot) <<" "
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