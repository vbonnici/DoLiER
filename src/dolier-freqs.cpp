/*
 * dolier-freqs.cpp
 *
 *  Created on: Nov 8, 2013
 *      Author: vbonnici
 */







#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include "data_ts.h"
#include "timer.h"
#include "trim.h"
#include "pars_t.h"

#include "DNA5Alphabet.h"
#include "FASTAReader.h"
#include "CstyleIndex.h"
#include "NSAIterator.h"
#include "Cs5DLIndex.h"
#include "main_common.h"
#include "SASearcher.h"
#include "FlatKmersIterator.h"

using namespace dolierlib;


//#define FREQS_DEBUG

char* scmd;
void pusage(){
	std::cout<<"Document listing problem. Output number of sequences where a kmer appears.\n";
	std::cout<<"Usage: "<<scmd<<" <seqs_file> [--kmer-list <kmers_file>] <o_file> <k> <m> [-r] [-p nthreads] [--kmers]\n";
	std::cout<<"\t<seqs_file> input FASTA file containing one or more sequences.\n";
	std::cout<<"\toptional [--kmer-list <kmers_file>] output only frequences of kmers listed in <kmers_file> (format one row one kmer).\n";
	std::cout<<"\t if no [--kmer-list] then output frequences of those kmers present in the set of input sequences.\n";
	std::cout<<"<t<o_file> output file.\n";
	std::cout<<"\t<k> kmer length. If [--kmer-list] then <k> will be ignored.\n";
	std::cout<<"\t<m> maximum number of allowed mismathes.\n";
	std::cout<<"\toptional [-r] do not output reverse complement. Given a kmer A and its reverse complement B, output A if A<=B otherwise output B even if B is not in input sequences. If [--kmer-list] then [-r] will be ignored.\n";
	std::cout<<"\t\tMoreover, the frequency of A is obtained considering all the occurrences of A plus all the occurrences of B.\n";
	std::cout<<"\toptional [-p nthreads] run in parallel using nthreads number of posix threads.\n";
	std::cout<<"\toptional [--kmers] also output kmers, not only their frequencies.\n";
	std::cout<<"If no [--kmer-list], frequencies are printed according to kmers lexicographic order.\n";
}

/*
Make uses of these theorems when dealing with reverse complements (-r).
- S a set of nucleotide sequences
- w a k-mer
- rc(w) the reverse complement of w
- rc(S) the reverse complements of S
- dl(w,S) the document listing of w in S, namely the set of sequences where w appears
- |dl(w,S)| the size of such set
We want to list every kmers that appears in S or rc(S).
We want to calculate
- frequency(w,S) = |dl(w,S)  U  dl(w,rc(S))|
- frequency(rc(w),S) = |dl(rc(w),S)  U  dl(rc(w),rc(S))|
We can proof that:
- frequency(w,S) = |dl(w,S)  U  dl(rc(w),S)|		so we don't need to store and index rc(S)
- frequency(w,S) = frequency(rc(w),S)				so we can calculate only one of them
 */


class mpars_t : public pars_t{
public:
	std::string iseqs;
	std::string klist_file;
	std::string ofile;
	int k;
	int m;
	bool opt_r;
	int nthreads;
	bool opt_kmers;

	mpars_t(int argc, char* argv[]) : pars_t(argc,argv){
		k = 0;
		m = 0;
		opt_r = false;
		nthreads = 1;
		opt_kmers = false;
	}
	~mpars_t(){}

	virtual void print(){
		std::cout<<"[ARG][iseqs]["<<iseqs<<"]\n";
		std::cout<<"[ARG][kmer_list]["<<klist_file<<"]\n";
		std::cout<<"[ARG][ofile]["<<ofile<<"]\n";
		std::cout<<"[ARG][k]["<<k<<"]\n";
		std::cout<<"[ARG][m]["<<m<<"]\n";
		std::cout<<"[ARG][r]["<<opt_r<<"]\n";
		std::cout<<"[ARG][nthreads]["<<nthreads<<"]\n";
		std::cout<<"[ARG][kmers]["<<opt_kmers<<"]\n";
	}

	virtual void usage(){
		pusage();
	}

	virtual void check(){
		print();
		std::cout<<"check...";
		iseqs = trim(iseqs);
		ofile = trim(ofile);

		if(		(iseqs.length() == 0)	||
				(ofile.length() == 0)	||
				(k < 1)					||
				(m < 0)					||
				(nthreads < 1)
		){
			usage();
			exit(1);
		}

	}

	virtual void parse(){
		std::string cmd;

		iseqs = next_string();

		cmd = next_string();
		if(cmd == "--kmer-list"){
			klist_file = next_string();
		}
		else{
			prev();
		}

		ofile = next_string();

		k = next_int();

		m = next_int();

		while(has_more()){
			cmd = next_string();
			if(cmd == "-r"){
				opt_r = true;
			}
			else if(cmd == "-p"){
				nthreads = next_int();
			}
			else if(cmd == "--kmers"){
				opt_kmers = true;
			}
			else{
				usage();
				exit(1);
			}
		}
		check();
	}

};






class Ikmer{
public:
	std::string kmer;
	usize_t value;
	Ikmer(std::string &kmer, usize_t freq){
		this->kmer = kmer;
		this->value = freq;
	}
	Ikmer(dna5_t *mer, int k, usize_t freq){
		this->kmer = to_string(mer,k);
		this->value = freq;
	}

	friend bool operator== (const Ikmer& lhs, const Ikmer& rhs){
		return lhs.kmer == rhs.kmer;
	}

	friend bool operator< (const Ikmer& lhs, const Ikmer& rhs) {
		return lhs.kmer < rhs.kmer;
	}
};






/*
 * Serial versions exist to avoid parallelization overhead in case of single thread.
 * They are also a basic example of how frequencies are computed.
 */
void run_from_list(mpars_t &pars);
void run_from_list_parallel(mpars_t &pars);
void run_serial(mpars_t &pars);
void run_parallel(mpars_t &pars);

int main(int argc, char* argv[]){
	TIMEHANDLE start_t = start_time();
	scmd = argv[0];

	mpars_t pars(argc, argv);
	pars.parse();


	if(pars.klist_file.length() >0){
		if(pars.nthreads == 1)
			run_from_list(pars);
		else
			run_from_list_parallel(pars);
	}
	else if(pars.nthreads == 1){
		run_serial(pars);
	}
	else{
		run_parallel(pars);
	}


	std::cout<<"total time "<<end_time(start_t)<<"\n";
	exit(0);
}











/*
 * get frequencies of a kmers list (from file) in serial
 */
void run_from_list(mpars_t &pars){
	DNA5MS_t f_ms(true);
		load_dna5ms(pars.iseqs, f_ms);
	CsFullyContainer f_ds(true);
		//build_fully_ds(f_ms, f_ds);
		build_strict_ds(f_ms, f_ds);
	Cs5DLIndex f_dlindex(f_ms.seq, f_ms.seq_length, f_ds.SA, f_ms.lengths, f_ms.nof_seqs);
	SASearcher sas(f_ms.seq, f_ds.SA, f_ms.seq_length);


	std::ofstream ofs;
	ofs.open((pars.ofile).c_str(), std::ios::out);
	if(!ofs.is_open() || ofs.bad()){
		std::cout<<"Error on opening output file : "<<pars.ofile<<" \n";
		exit(1);
	}


	int m = pars.m;


	usize_t f_dcount = 0;
	usize_t idoc_i;
	bool *f_idocs = new bool[f_ms.nof_seqs];


	usize_t total_kmers = 0;
	usize_t present_kmers = 0;

	double search_time = 0;
	double write_time = 0;
	TIMEHANDLE start_t = start_time();

	std::set<Ikmer> kmerlist;


	FlatKmersIterator it(pars.klist_file);
	std::string s;
	int k;
	while(it.next()){
		//it.get_kmer(it_codes);
		s = it.get_kmer();
		k = s.length();

		dna5_t *it_codes = to_dna5(s);
		dolierlib::ssize_t *D = new dolierlib::ssize_t[k];
//		char *it_syms = new char[k+1];
//		it_syms[k] = '\0';

		total_kmers++;


		memset(f_idocs, 0, f_ms.nof_seqs);

		start_t = start_time();
		f_dlindex.dl_count_wD(it_codes, D, k, m, f_idocs);
		search_time += end_time(start_t);

		if(pars.opt_r){
			dna5_t *it_codes_rev = new dna5_t[k];
			reverseComplement(it_codes, it_codes_rev, k);
			if(compare(it_codes_rev, it_codes, k) != 0){

				start_t = start_time();
				f_dlindex.dl_count_wD(it_codes_rev, D, k, m, f_idocs);
				search_time += end_time(start_t);
			}
			delete [] it_codes_rev;
		}

		f_dcount = 0;
		for(idoc_i=0; idoc_i<f_ms.nof_seqs; idoc_i++)
			if(f_idocs[idoc_i]) f_dcount++;


//		if(f_dcount > 0){
			present_kmers++;
			start_t = start_time();
			//DNA5Alphabet::toSymbolArray(it_codes, it_syms, k);
			if(pars.opt_kmers)
				ofs<<s<<"\t"<<f_dcount<<"\n";
			else
				ofs<<f_dcount<<"\n";
			write_time += end_time(start_t);
//		}


		delete [] it_codes;

	}

	std::cout<<"nof kmers "<<total_kmers<<"\n";
	std::cout<<"present kmers "<<present_kmers<<"\n";
	std::cout<<"write time "<<write_time<<"\n";

	ofs.flush();
	ofs.close();

	delete [] f_idocs;
}











class KmerFreq_t{
public:
	std::string kmer;
	usize_t freq;
	KmerFreq_t(std::string kmer, usize_t freq){
		this->kmer =kmer;
		this->freq = freq;
	}
	KmerFreq_t(const KmerFreq_t& e){
		this->kmer = e.kmer;
		this->freq = e.freq;
	}
	friend bool operator== (const KmerFreq_t& lhs, const KmerFreq_t& rhs){
		return lhs.kmer == rhs.kmer;
	}

	friend bool operator< (const KmerFreq_t& lhs, const KmerFreq_t& rhs) {
		return lhs.kmer < rhs.kmer;
	}
};


class kmerlist_thread_t{
public:
	int idx;
	int k;
	int m;
	bool reverse;
	DNA5MS_t &f_ms;
	CsFullyContainer &f_ds;
	usize_t **ranges;
	double *timers;
	std::vector<KmerFreq_t> &kfreqs_list;

	Cs5DLIndex *f_dlindex;

	kmerlist_thread_t(
			int _idx,
			int _k,
			int _m,
			bool _reverse,
			DNA5MS_t &_f_ms,
			CsFullyContainer &_f_ds,
			usize_t **_ranges,
			double *_timers,
			std::vector<KmerFreq_t> &_kfreqs_list
			)
	: idx(_idx), k(_k), m(_m), reverse(_reverse), f_ms(_f_ms), f_ds(_f_ds), ranges(_ranges), timers(_timers), kfreqs_list(_kfreqs_list)
	{
		f_dlindex = new Cs5DLIndex (f_ms.seq, f_ms.seq_length, f_ds.SA, f_ms.lengths, f_ms.nof_seqs);
	}
	~kmerlist_thread_t(){
		delete f_dlindex;
	}

	static void* run(void* argsptr){
		TIMEHANDLE start_t = start_time();

		kmerlist_thread_t* st = (kmerlist_thread_t*)argsptr;

		int idx = st->idx;
		int k = st->k;

		dna5_t *it_codes;// = new dna5_t[k];
		dna5_t *it_codes_rev = new dna5_t[k];
		dolierlib::ssize_t *D = new dolierlib::ssize_t[k];
		usize_t f_dcount = 0;
		usize_t idoc_i;
		bool* f_idocs = new bool[st->f_ms.nof_seqs];

		for(usize_t r=st->ranges[idx][0]; r<st->ranges[idx][1]; r++){
			//KmerFreq_t ele = st->kfreqs_list[r];
			//memcpy(it_codes,  , k * sizeof(dna5_t));
			it_codes = to_dna5(st->kfreqs_list[r].kmer);
			memset(f_idocs, 0, st->f_ms.nof_seqs * sizeof(bool));

			st->f_dlindex->dl_count_wD(it_codes, D, k, st->m, f_idocs);
			if(st->reverse){
				reverseComplement(it_codes, it_codes_rev, k);
				if(compare(it_codes_rev, it_codes, k) != 0){//added last time
					st->f_dlindex->dl_count_wD(it_codes_rev, D, k, st->m, f_idocs);
				}
			}

			f_dcount = 0;
			for(idoc_i=0; idoc_i<st->f_ms.nof_seqs; idoc_i++){
				if(f_idocs[idoc_i]) f_dcount++;
			}

			st->kfreqs_list[r].freq = f_dcount;

			delete [] it_codes;
		}

		//delete [] it_codes;
		delete [] it_codes_rev;
		delete [] D;
		delete [] f_idocs;

		st->timers[idx] = end_time(start_t);

		pthread_exit(NULL);
	}

};


/*
 * get frequencies of a kmers list (from file) in parallel
 */
void run_from_list_parallel(mpars_t &pars){
	std::vector<KmerFreq_t> kfreqs_list;

	FlatKmersIterator it(pars.klist_file);
	while(it.next()){
		kfreqs_list.push_back(KmerFreq_t(it.get_kmer(), 0));
	}
	int ratio = (int)floor(static_cast<double>(kfreqs_list.size()) /static_cast<double>(pars.nthreads));
	usize_t **t_ranges = new usize_t*[pars.nthreads];
	for(int i=0; i<pars.nthreads; i++)
		t_ranges[i] = new usize_t[2];

	if(pars.nthreads < 2){
		//never enter here
		t_ranges[0][0] = 0;
		t_ranges[0][1] = kfreqs_list.size();
	}
	else{
		t_ranges[0][0] = 0;
		t_ranges[0][1] = ratio;
		for(int t = 1; t<pars.nthreads; t++){
			t_ranges[t][0] = t_ranges[t-1][1];
			t_ranges[t][1] = t_ranges[t][0] + ratio;
		}
		t_ranges[pars.nthreads-1][1] = kfreqs_list.size();
	}

	std::cout<<"workload: ";
	for(int t=0; t<pars.nthreads; t++){
		std::cout<<"["<<t_ranges[t][0]<<","<<t_ranges[t][1]<<"] ";
	}
	std::cout<<"\n";


	DNA5MS_t f_ms(true);
		load_dna5ms(pars.iseqs, f_ms);
	CsFullyContainer f_ds(true);
		build_strict_ds(f_ms, f_ds);


	double *timers = new double[pars.nthreads];
	memset(timers, 0, pars.nthreads * sizeof(double));


	pthread_t *threads = new pthread_t[pars.nthreads];
	int rc;
	for(int i=0;i<pars.nthreads;i++){
		kmerlist_thread_t *st = new kmerlist_thread_t(i, pars.k, pars.m, pars.opt_r, f_ms, f_ds, t_ranges, timers, kfreqs_list);
		rc = pthread_create(&threads[i], NULL, kmerlist_thread_t::run, (void*)st);
		if(rc){
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(1);
		}
	}
	TIMEHANDLE start_t = start_time();
	for(int i=0;i<pars.nthreads;i++){
		rc = pthread_join(threads[i], NULL);
		if(rc){
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(1);
		}
	}
	std::cout<<"searching done "<<end_time(start_t)<<"\n";


	std::ofstream ofs;
	ofs.open((pars.ofile).c_str(), std::ios::out);
	if(!ofs.is_open() || ofs.bad()){
		std::cout<<"Error on opening output file : "<<pars.ofile<<" \n";
		exit(1);
	}
	double write_time = 0;
	if(pars.opt_kmers){
		for(std::vector<KmerFreq_t>::iterator IT = kfreqs_list.begin(); IT!=kfreqs_list.end(); IT++){
			ofs<<(*IT).kmer<<"\t"<<(*IT).freq<<"\n";
		}
	}
	else{
		for(std::vector<KmerFreq_t>::iterator IT = kfreqs_list.begin(); IT!=kfreqs_list.end(); IT++){
			ofs<<(*IT).freq<<"\n";
		}
	}
	ofs.flush();
	ofs.close();


	std::cout<<"timers:\t";
	for(int i=0; i<pars.nthreads; i++)
		std::cout<<timers[i]<<"\t";
	std::cout<<"\n";

	std::cout<<"write time "<<write_time<<"\n";

	delete [] t_ranges;
	delete [] threads;
}










































/*
 * get frequencies of a all distinct kmers present in the set of sequences, in serial
 */
void run_serial(mpars_t &pars){
	DNA5MS_t f_ms(true);
		load_dna5ms(pars.iseqs, f_ms);
	CsFullyContainer f_ds(true);
		build_fully_ds(f_ms, f_ds);
	Cs5DLIndex f_dlindex(f_ms.seq, f_ms.seq_length, f_ds.SA, f_ms.lengths, f_ms.nof_seqs);
	SASearcher sas(f_ms.seq, f_ds.SA, f_ms.seq_length);


	std::ofstream ofs;
	ofs.open((pars.ofile).c_str(), std::ios::out);
	if(!ofs.is_open() || ofs.bad()){
		std::cout<<"Error on opening output file : "<<pars.ofile<<" \n";
		exit(1);
	}


	int k = pars.k;
	int m = pars.m;

	char *it_syms = new char[k+1];
	it_syms[k] = '\0';
	dna5_t *it_codes = new dna5_t[k];
	dna5_t *it_codes_rev = new dna5_t[k];
	dolierlib::ssize_t *D = new dolierlib::ssize_t[k];

	usize_t f_dcount = 0;
	usize_t idoc_i;
	bool *f_idocs = new bool[f_ms.nof_seqs];

	bool *donot = NULL;
	if(pars.opt_r){
		donot = new bool[f_ms.seq_length];
		memset(donot, false, f_ms.seq_length);
	}
	usize_t donoti;

	usize_t total_kmers = 0;
	usize_t dis_donot = 0;

	double write_time = 0;
	TIMEHANDLE start_t = start_time();

	std::set<Ikmer> kmerlist;


	NSAIterator it = NSAIterator::begin(f_ms, f_ds, k);
	while(it.next()){
		it.get_kmer(it_codes);

//		std::cout<<">";
//		print_dna5(it_codes,k);

#ifdef  FREQS_DEBUG
		print_dna5(it_codes, k);
		if(!pars.opt_r || (donot[it.i_start]==false))std::cout<<".";
		std::cout<<"\n";
#endif

		total_kmers++;

		if((!pars.opt_r) || (donot[it.i_start]==false)){
//			std::cout<<" + ";

			memset(f_idocs, 0, f_ms.nof_seqs);

			donoti = f_ms.seq_length;

			f_dlindex.dl_count_wD(it_codes, D, k, m, f_idocs);
			if(pars.opt_r){
//				std::cout<<" . ";
				reverseComplement(it_codes, it_codes_rev, k);


				//if(compare(it_codes_rev, it_codes, k) > 0){
				if(compare(it_codes_rev, it_codes, k) != 0){
					f_dlindex.dl_count_wD(it_codes_rev, D, k, m, f_idocs);
					donoti = sas.contains(it_codes_rev, k);
					if(donoti < f_ms.seq_length){
						donot[donoti] = true;
						dis_donot++;
					}
				}
			}

			f_dcount = 0;
			for(idoc_i=0; idoc_i<f_ms.nof_seqs; idoc_i++)
				if(f_idocs[idoc_i]) f_dcount++;

			if(pars.opt_r){
				if(compare(it_codes_rev, it_codes, k)<0){
					kmerlist.insert( Ikmer(it_codes_rev, k, f_dcount) );
				}
				else{
					kmerlist.insert( Ikmer(it_codes, k, f_dcount) );
				}
			}
			else{
				start_t = start_time();
				if(pars.opt_kmers){
					DNA5Alphabet::toSymbolArray(it_codes, it_syms, k);
					ofs<<it_syms<<"\t"<<f_dcount<<"\n";
				}
				else
					ofs<<f_dcount<<"\n";
				write_time += end_time(start_t);
			}

		}
//		else{
//			std::cout<<" - ";
//		}

//		std::cout<<"\n";
	}

	if(pars.opt_r){
//		it = NSAIterator::begin(f_ms, f_ds, k);
//		while(it.next()){
//			it.get_kmer(it_codes);
//			DNA5Alphabet::toSymbolArray(it_codes, it_syms, k);
//			std::cout<<it_syms<<"\t"<<donot[it.i_start]<<"\n";
//		}


		start_t = start_time();
		if(pars.opt_kmers){
			for(std::set<Ikmer>::iterator IT = kmerlist.begin(); IT!=kmerlist.end(); IT++){
				ofs<< (*IT).kmer <<"\t"<< (*IT).value <<"\n";
			}
		}
		else{
			for(std::set<Ikmer>::iterator IT = kmerlist.begin(); IT!=kmerlist.end(); IT++){
				ofs<<(*IT).value <<"\n";
			}
		}
		write_time += end_time(start_t);
	}

	std::cout<<"nof kmers "<<total_kmers<<"\n";
	if(pars.opt_r)
		std::cout<<"rc discarded "<<dis_donot<<"\n";
	std::cout<<"write time "<<write_time<<"\n";

	ofs.flush();
	ofs.close();

	delete [] it_syms;
	delete [] it_codes;
	delete [] it_codes_rev;
	delete [] D;
	delete [] f_idocs;
	if(donot != NULL)
		delete [] donot;
}














class freq_thread_t{
public:
	int idx;
	int k;
	int m;
	bool reverse;
	DNA5MS_t &f_ms;
	CsFullyContainer &f_ds;
	Cs5DLIndex *f_dlindex;
	usize_t *strict_SA;
	usize_t *freqs;
	usize_t **ranges;
	double *timers;

	freq_thread_t(
			int _idx,
			int _k,
			int _m,
			bool _reverse,
			DNA5MS_t &_f_ms,
			CsFullyContainer &_f_ds,
			usize_t *_strict_SA,
			usize_t *_freqs,
			usize_t **_ranges,
			double *_timers
			)
	: idx(_idx), k(_k), m(_m), reverse(_reverse), f_ms(_f_ms), f_ds(_f_ds), strict_SA(_strict_SA), freqs(_freqs), ranges(_ranges), timers(_timers)
	{
		f_dlindex = new Cs5DLIndex (f_ms.seq, f_ms.seq_length, f_ds.SA, f_ms.lengths, f_ms.nof_seqs);
	}
	~freq_thread_t(){
		delete f_dlindex;
	}


	static void* run(void* argsptr){
		TIMEHANDLE start_t = start_time();

		freq_thread_t* st = (freq_thread_t*)argsptr;

		int idx = st->idx;
		int k = st->k;

		dna5_t *it_codes = new dna5_t[k];
		dna5_t *it_codes_rev = new dna5_t[k];
		dolierlib::ssize_t *D = new dolierlib::ssize_t[k];
		usize_t f_dcount = 0;
		usize_t idoc_i;
		bool* f_idocs = new bool[st->f_ms.nof_seqs];

		for(usize_t r = st->ranges[idx][0]; r < st->ranges[idx][1]; r++){
			memcpy(it_codes, &(st->f_ms.seq[ st->f_ds.SA[ st->strict_SA[r] ] ]) , k * sizeof(dna5_t));
			memset(f_idocs, 0, st->f_ms.nof_seqs * sizeof(bool));

			st->f_dlindex->dl_count_wD(it_codes, D, k, st->m, f_idocs);
			if(st->reverse){
				reverseComplement(it_codes, it_codes_rev, k);
				st->f_dlindex->dl_count_wD(it_codes_rev, D, k, st->m, f_idocs);
			}

			f_dcount = 0;
			for(idoc_i=0; idoc_i<st->f_ms.nof_seqs; idoc_i++){
				if(f_idocs[idoc_i]) f_dcount++;
			}

			st->freqs[r] = f_dcount;

		}

		delete [] it_codes;
		delete [] it_codes_rev;
		delete [] D;
		delete [] f_idocs;

		st->timers[idx] = end_time(start_t);

		pthread_exit(NULL);
	}
};



/*
 * get frequencies of a all distinct kmers present in the set of sequences, in parallel
 */
void run_parallel(mpars_t &pars){
	DNA5MS_t f_ms(true);
		load_dna5ms(pars.iseqs, f_ms);
	CsFullyContainer f_ds(true);
		build_fully_ds(f_ms, f_ds);
	Cs5DLIndex f_dlindex(f_ms.seq, f_ms.seq_length, f_ds.SA, f_ms.lengths, f_ms.nof_seqs);
	SASearcher sas(f_ms.seq, f_ds.SA, f_ms.seq_length);

	int k = pars.k;

	NSAIterator it = NSAIterator::begin(f_ms, f_ds, k);

//	usize_t nof_kmers = nof(it);
//	std::cout<<"nof kmers: "<<nof_kmers<<"\n";

	char *it_syms = new char[k+1];
	it_syms[k] = '\0';
	dna5_t *it_codes = new dna5_t[k];
	dna5_t *it_codes_rev = new dna5_t[k];
	usize_t nof_kmers = 0;
	usize_t real_kmers = 0;

	bool *donot = NULL;
	if(pars.opt_r){
		donot = new bool[f_ms.seq_length];
		memset(donot, false, f_ms.seq_length * sizeof(bool));
	}
	usize_t donoti;

	if(!pars.opt_r){
		nof_kmers = nof(it);
		real_kmers = nof_kmers;
	}
	else{
		while(it.next()){
			nof_kmers++;
			it.get_kmer(it_codes);
			if(donot[it.i_start] == false){
				real_kmers++;

				reverseComplement(it_codes, it_codes_rev, k);
				if(compare(it_codes_rev, it_codes, k) > 0){
					donoti = sas.contains(it_codes_rev, k);
					if(donoti < f_ms.seq_length){
						donot[donoti] = true;
					}
				}
			}
		}
		it = NSAIterator::begin(f_ms, f_ds, k);
	}

	std::cout<<"nof kmers "<<nof_kmers<<"\n";
	if(pars.opt_r)
		std::cout<<"rc discared "<<(nof_kmers - real_kmers)<<"\n";

//	it = NSAIterator::begin(f_ms, f_ds, k);
//	while(it.next()){
//		it.get_kmer(it_codes);
//		DNA5Alphabet::toSymbolArray(it_codes, it_syms, k);
//		std::cout<<it_syms<<"\t"<<donot[it.i_start]<<"\n";
//	}

	usize_t *strict_SA = new usize_t[real_kmers];
	usize_t *freqs = new usize_t[real_kmers];
	int i = 0;
	it = NSAIterator::begin(f_ms, f_ds, k);
	while(it.next()){
		if(!pars.opt_r || (donot[it.i_start] == false)){
			strict_SA[i] = it.i_start;
			i++;
		}
	}

	int nthreads = pars.nthreads;

	int ratio = (int)floor(static_cast<double>(real_kmers) /static_cast<double>(nthreads));
	usize_t **t_ranges = new usize_t*[nthreads];
	for(int i=0; i<nthreads; i++)
		t_ranges[i] = new usize_t[2];


	if(nthreads < 2){
		//never enter here
		t_ranges[0][0] = 0;
		t_ranges[0][1] = real_kmers;
	}
	else{
		t_ranges[0][0] = 0;
		t_ranges[0][1] = ratio;
		for(int t = 1; t<nthreads; t++){
			t_ranges[t][0] = t_ranges[t-1][1];
			t_ranges[t][1] = t_ranges[t][0] + ratio;
		}
		t_ranges[nthreads-1][1] = real_kmers;
	}

	std::cout<<"workload: ";
	for(int t=0; t<nthreads; t++){
		std::cout<<"["<<t_ranges[t][0]<<","<<t_ranges[t][1]<<"] ";
	}
	std::cout<<"\n";


	double *timers = new double[nthreads];
	memset(timers, 0, nthreads * sizeof(double));


	pthread_t *threads = new pthread_t[nthreads];
	int rc;
	for(int i=0;i<nthreads;i++){
		freq_thread_t *st = new freq_thread_t(i, k, pars.m, pars.opt_r, f_ms, f_ds, strict_SA, freqs, t_ranges, timers);
		rc = pthread_create(&threads[i], NULL, freq_thread_t::run, (void*)st);
		if(rc){
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(1);
		}
	}
	TIMEHANDLE start_t = start_time();
	for(int i=0;i<nthreads;i++){
		rc = pthread_join(threads[i], NULL);
		if(rc){
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(1);
		}
	}

	std::cout<<"searching done "<<end_time(start_t)<<"\n";

	std::ofstream ofs;
	ofs.open((pars.ofile).c_str(), std::ios::out);
	if(!ofs.is_open() || ofs.bad()){
		std::cout<<"Error on opening output file : "<<pars.ofile<<" \n";
		exit(1);
	}




	double write_time = 0;
	double sort_time  =0;

	std::set<Ikmer> kmerlist;

//	char *it_syms = new char[k+1];
//	it_syms[k] = '\0';
//	dna5_t *it_codes = new dna5_t[k];
	it = NSAIterator::begin(f_ms, f_ds, k);
	usize_t r = 0;
	while(it.next()){
		if(!pars.opt_r || (donot[it.i_start] == false)){
//			start_t = start_time();
			it.get_kmer(it_codes);
//			DNA5Alphabet::toSymbolArray(it_codes, it_syms, k);
//			ofs<<it_syms<<"\t"<<freqs[r]<<"\n";
//			r++;
//			write_time += end_time(start_t);

			if(pars.opt_r){
				start_t = start_time();
				reverseComplement(it_codes, it_codes_rev, k);
				if(compare(it_codes_rev, it_codes, k)<0){
					kmerlist.insert( Ikmer(it_codes_rev, k, freqs[r]) );
				}
				else{
					kmerlist.insert( Ikmer(it_codes, k, freqs[r]) );
				}
				sort_time += end_time(start_t);
			}
			else{
				start_t = start_time();
				if(pars.opt_kmers){
					DNA5Alphabet::toSymbolArray(it_codes, it_syms, k);
					ofs<<it_syms<<"\t"<<freqs[r]<<"\n";
				}
				else
					ofs<<freqs[r]<<"\n";
				write_time += end_time(start_t);
			}
			r++;
		}
	}
	if(pars.opt_r){
		start_t = start_time();
		if(pars.opt_kmers){
			for(std::set<Ikmer>::iterator IT = kmerlist.begin(); IT!=kmerlist.end(); IT++){
				ofs<< (*IT).kmer <<"\t"<< (*IT).value <<"\n";
			}
		}
		else{
			for(std::set<Ikmer>::iterator IT = kmerlist.begin(); IT!=kmerlist.end(); IT++){
				ofs<<(*IT).value <<"\n";
			}
		}
		write_time += end_time(start_t);
	}

	ofs.flush();
	ofs.close();

	delete [] strict_SA;
	delete [] freqs;
	delete [] t_ranges;
	delete [] threads;
	delete [] it_syms;
	delete [] it_codes;
	delete [] it_codes_rev;
	if(donot != NULL)
		delete [] donot;


	std::cout<<"timers:\t";
	for(int i=0; i<nthreads; i++)
		std::cout<<timers[i]<<"\t";
	std::cout<<"\n";
	std::cout<<"write time "<<write_time<<"\n";
	if(pars.opt_r)
		std::cout<<"sort time "<<sort_time<<"\n";
}
