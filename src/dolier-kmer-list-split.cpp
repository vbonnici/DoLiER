/*
 * dolier-list-split.cpp
 *
 *  Created on: Nov 14, 2013
 *      Author: vbonnici
 */



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdint.h>


void pusage(){
	std::cout<<"Usage: cmd <ifile> <oprefix> <n>\n";
}

int main(int argc, char* argv[]){
	if(argc != 4){
		pusage();
		exit(1);
	}

	std::string ifile = argv[1];
	std::string oprefix = argv[2];
	double n = atof(argv[3]);

	if(n < 2){
		pusage();
		exit(1);
	}


	double total = 0;

	std::string s;

	std::ifstream ifs;
	ifs.open((ifile).c_str(), std::ios::in);
	if(!ifs.is_open() || ifs.bad()){
		std::cout<<"Error on opening input file : "<<ifile<<" \n";
		exit(1);
	}
	while(ifs >> s){
		total++;
	}
	ifs.close();


	double each = total / n;

	std::cout<<"packing on "<<each<<"\n";


	ifs.open((ifile).c_str(), std::ios::in);
	if(!ifs.is_open() || ifs.bad()){
		std::cout<<"Error on opening input file : "<<ifile<<" \n";
		exit(1);
	}



	for(double i=0; i<n-1; i++){
		std::ofstream ofs;
		std::stringstream ss;
		ss<<oprefix<<i;
		ofs.open((ss.str()).c_str(), std::ios::out);
		if(!ofs.is_open() || ofs.bad()){
			std::cout<<"Error on opening output file : "<<oprefix<<" \n";
			exit(1);
		}

		for(double j=0; j<each; j++){
			ifs >> s;
			ofs << s <<"\n";
		}

		ofs.flush();
		ofs.close();
	}


	std::ofstream ofs;
	std::stringstream ss;
	ss<<oprefix<<(n-1);
	ofs.open((ss.str()).c_str(), std::ios::out);
	if(!ofs.is_open() || ofs.bad()){
		std::cout<<"Error on opening output file : "<<(ss.str())<<" \n";
		exit(1);
	}

	while(ifs >> s){
		ofs << s <<"\n";
	}

	ofs.flush();
	ofs.close();



	ifs.close();
}
