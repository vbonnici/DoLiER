/*
 * FASTA_reader.h
 *
 *  Created on: Sep 26, 2013
 *      Author: vbonnici
 */




/*
 * Read FASTA file.
 */

#ifndef FASTAREADER_H_
#define FASTAREADER_H_


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include <bitset>
#include <vector>

#include "DNA5Alphabet.h"

#include "trim.h"



namespace dolierlib{

class FASTAReader{


	std::string _fname;
	std::ifstream ifs;

	std::string _line;
	std::string _seq;

public:
	FASTAReader(){
	}

	FASTAReader(const std::string& fname){
		open(fname);
	}

	void open(const std::string& fname) {
		_fname = fname;
		ifs.open(fname.c_str());
		if(!ifs.is_open() || ifs.bad()){
			std::cout<<"Error on opening input file : "<<fname<<" \n";
			exit(1);
		}
		_line = "";
	}

	void close(){
		if(ifs.is_open())
			ifs.close();
		_line = "";
	}


	/*
	 * read next sequence in the file
	 */
	bool nextSeq(
			std::string* name, //output pointer
			std::string* seq //output pointer
			){
		if(_line == ""){
			if(!getline(ifs, _line))
				return false;
		}
		*name = *(new std::string(trim(_line).substr(1)));

		bool readed = false;
		*seq = *(new std::string(""));
		while(getline(ifs, _line) && _line[0]!='>'){
			trim(_line);
			if(_line.length() > 0){
				readed = true;
				*seq += _line;
			}
		}

		if(ifs.eof() || ifs.bad()){
			_line = "";
		}

		return readed;
	}

	usize_t 
	get_huge_length(){
		if(!getline(ifs, _line))
			return 0;

		usize_t len = 0;

		while(getline(ifs, _line) && _line[0]!='>'){
			trim(_line);
			len += _line.length();
		}

		return len;
	}

	bool
	read_huge(dna5_t *seq, usize_t seq_length){
		if(!getline(ifs, _line))
			return false;

		usize_t i=0, j=0;
		while(getline(ifs, _line) && _line[0]!='>'){
			trim(_line);
			for(j=0; j<_line.length(); j++, i++){
				seq[i] =  DNA5Alphabet::codeFor(_line[j]);
			}
		}
		return i==seq_length;
	};


	/*
	 * read all the sequences and put them in the output vector
	 */
	void readAll(std::vector<std::string>& sequences){

		std::string name;
		std::string seq;

		while(nextSeq(&name,&seq)){
			sequences.push_back(toUpperCase(seq));
		}
	}
//
//	void readAll(std::vector<char*>& sequences){
//		std::string name;
//		std::string seq;
//		while(nextSeq(&name,&seq)){
//			char* cc = new char[seq.length()];
//			seq.copy(cc, seq.length(),0);
//			sequences.push_back(cc);
//		}
//	}

//	void readAll(std::vector<DNA5Sequence*>& sequences){
//		std::string name;
//		std::string seq;
//		while(nextSeq(&name,&seq)){
//			sequences.push_back(new DNA5Sequence(seq));
//		}
//	}


	/*
	 * read all the sequences in dna5_t output format
	 */
	void readAll(
			std::vector<dna5_t*>& sequences, //output sequences
			std::vector<usize_t>& lengths,   //output sequences lengths
			bool appendUndef //if true then append $ to the end of each sequence
			){
		std::string name;
		std::string seq;
		dna5_t *codes;

		int i,l;
		while(nextSeq(&name,&seq)){
			if(appendUndef){
				l = seq.length();
				codes = new dna5_t[l+1];
				for(i=0;i<l;i++){
					codes[i] = DNA5Alphabet::codeFor(seq[i]);
				}
				codes[l] = DNA5Alphabet::UNDEF;

				sequences.push_back(codes);
				lengths.push_back(l+1);
			}
			else{
				l = seq.length();
				codes = new dna5_t[l];
				for(i=0;i<l;i++){
					codes[i] = DNA5Alphabet::codeFor(seq[i]);
				}

				sequences.push_back(codes);
				lengths.push_back(l);
			}
		}
	}


};

}


#endif /* FASTAREADER_H_ */
