/*
 * pars_t.h
 *
 *  Created on: Dec 17, 2013
 *      Author: vbonnici
 */

#ifndef PARS_T_H_
#define PARS_T_H_


#include "data_ts.h"


/*
 * An abstract class for command line arguments parsing.
 * See any executable source, using it, in src for an example.
 */

namespace dolierlib{


class pars_t{
private:
	int a_i;
	int argc;
	char** argv;

public:

	pars_t(int argc, char* argv[]){
		this->argc = argc;
		this->argv = argv;
		a_i = 1;
	}
	virtual ~pars_t(){}

protected:
	void assert_ai(){
		if(a_i<1 || a_i >= argc){
			usage();
			exit(1);
		}
	}
	bool has_more(){
		return a_i < argc;
	}
	void prev(){
		a_i--;
	}


	void parse_int(int *n, char* str){
		std::istringstream iss(str);
		if(!(iss>>*n)){
			usage();
			exit(1);
		}
	}
	void parse_float(float *n, char* str){
		std::istringstream iss(str);
		if(!(iss>>*n)){
			usage();
			exit(1);
		}
	}
	void parse_double(double *n, char* str){
		std::istringstream iss(str);
		if(!(iss>>*n)){
			usage();
			exit(1);
		}
	}
	void parse_size_t(size_t *n, char* str){
		std::istringstream iss(str);
		if(!(iss>>*n)){
			usage();
			exit(1);
		}
	}

	void parse_usize_t(usize_t *n, char* str){
		std::istringstream iss(str);
		if(!(iss>>*n)){
			usage();
			exit(1);
		}
	}

	std::string next_string(){
		assert_ai();
		std::string s = argv[a_i];
		a_i++;
		return s;
	}
	int next_int(){
		assert_ai();
		int r = 0;
		parse_int(&r, argv[a_i]);
		a_i++;
		return r;
	}
	double next_double(){
		assert_ai();
		double r = 0;
		parse_double(&r, argv[a_i]);
		a_i++;
		return r;
	}
	float next_float(){
		assert_ai();
		float r = 0;
		parse_float(&r, argv[a_i]);
		a_i++;
		return r;
	}
	size_t next_size_t(){
		assert_ai();
		size_t r = 0;
		parse_size_t(&r, argv[a_i]);
		a_i++;
		return r;
	}

	usize_t next_usize_t(){
		assert_ai();
		usize_t r = 0;
		parse_usize_t(&r, argv[a_i]);
		a_i++;
		return r;
	}

public:
//	bool is_one_of(std::string p, std::string ... ss){
//		for(std::string s : ss){
//			if(s == p){
//				return true;
//			}
//		}
//		return false;
//	}
//	void check_one_of(std::string p, std::string ... ss){
//		if(!is_one_of(p, ss)){
//			usage();
//			exit(1);
//		}
//	}

public:
	virtual void check() =0;
	virtual void parse() =0;
	virtual void usage() {std::cout<<"Wrong parameters!";};
	virtual void print() =0;
};

}



#endif /* PARS_T_H_ */
