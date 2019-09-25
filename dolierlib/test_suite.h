/*
 * test_suite.h
 *
 *  Created on: Oct 16, 2013
 *      Author: vbonnici
 */

#ifndef TEST_SUITE_H_
#define TEST_SUITE_H_


/*
 * a set of functions for binomial tests, or other distribution
 */

//#define DEBUG_TEST_SUITE_H_


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <vector>
#include <limits>

#include "data_ts.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


namespace dolierlib{


/*
 * Pietro's binomial test
 */
inline
double dolier_binomial(const usize_t f_dcount, const usize_t f_total, const usize_t b_dcount, const usize_t b_total){
	double  p = (double)b_dcount / (double)b_total;
	double result = 0;
	for(usize_t i = f_dcount; i<=f_total; i++){
		result += gsl_ran_binomial_pdf(static_cast<unsigned int>(i), p, static_cast<unsigned int>(f_total));
	}
	return result;
}


/*
 * a more efficient binomial test, thanks to gsl
 */
inline
double dolier_binomial_Q(const usize_t f_dcount, const usize_t f_total, const usize_t b_dcount, const usize_t b_total){
	double  p = (double)b_dcount / (double)b_total;
	//return gsl_cdf_binomial_Q(static_cast<unsigned int>(f_dcount), p, static_cast<unsigned int>(f_total));
	return gsl_ran_binomial_pdf(static_cast<unsigned int>(f_dcount), p, static_cast<unsigned int>(f_total))
			+ gsl_cdf_binomial_Q(static_cast<unsigned int>(f_dcount), p, static_cast<unsigned int>(f_total));
}


/*
 * Pietro's binomial test with Bon Ferroni correction
 */
inline
double dolier_binomial_bonf(const usize_t f_dcount, const usize_t f_total, const usize_t b_dcount, const usize_t b_total, const usize_t nof_kmers){
	double p = (double)b_dcount / (double)b_total;
	double result = 0;
	for(usize_t i = f_dcount; i<=f_total; i++){
		result += gsl_ran_binomial_pdf(static_cast<unsigned int>(i), p, static_cast<unsigned int>(f_total));
	}
	return double(nof_kmers) * result;
}


/*
 * a more efficient binomial test with Bon Ferroni correction, thanks to gsl
 */
inline
double dolier_binomial_bonf_Q(const usize_t f_dcount, const usize_t f_total, const usize_t b_dcount, const usize_t b_total, const usize_t nof_kmers){
	double  p = (double)b_dcount / (double)b_total;
	//return gsl_cdf_binomial_Q(static_cast<unsigned int>(f_dcount), p, static_cast<unsigned int>(f_total)) * nof_kmers;
	return (
			gsl_ran_binomial_pdf(static_cast<unsigned int>(f_dcount), p, static_cast<unsigned int>(f_total))
			+
			gsl_cdf_binomial_Q(static_cast<unsigned int>(f_dcount), p, static_cast<unsigned int>(f_total))
			)
			* nof_kmers;
}























/*
 * other distributions tried during tests
 */




inline
double dolier_poisson_Q(const usize_t f_dcount, const usize_t f_total, const usize_t b_dcount, const usize_t b_total){
	double  p = (double)b_dcount / (double)b_total;
	double mu = p *  static_cast<double>(f_total);
	return gsl_ran_poisson_pdf(static_cast<unsigned int>(f_dcount),mu)
			+ gsl_cdf_poisson_Q(static_cast<unsigned int>(f_dcount), mu);
//	return gsl_ran_binomial_pdf(static_cast<unsigned int>(f_dcount), p, static_cast<unsigned int>(f_total))
//			+ gsl_cdf_binomial_Q(static_cast<unsigned int>(f_dcount), p, static_cast<unsigned int>(f_total));
}

inline
double dolier_poisson_bonf_Q(const usize_t f_dcount, const usize_t f_total, const usize_t b_dcount, const usize_t b_total, const usize_t nof_kmers){
	double  p = (double)b_dcount / (double)b_total;
	double mu = p *  static_cast<double>(f_total);
	if(mu==0) std::cout<<"#mu=0#\t"<<f_dcount<<"\t"<<f_total<<"\t"<<b_dcount<<"\t"<<b_total<<"\t"<<p<<"\n";
	return 		(gsl_ran_poisson_pdf(static_cast<unsigned int>(f_dcount),mu)
				+ gsl_cdf_poisson_Q(static_cast<unsigned int>(f_dcount), mu))
			* nof_kmers;
//	return (
//			gsl_ran_binomial_pdf(static_cast<unsigned int>(f_dcount), p, static_cast<unsigned int>(f_total))
//			+
//			gsl_cdf_binomial_Q(static_cast<unsigned int>(f_dcount), p, static_cast<unsigned int>(f_total))
//			)
//			* nof_kmers;
}



inline
double dolier_poisson_mu1_Q(const usize_t f_dcount, const usize_t f_total, const usize_t b_dcount, const usize_t b_total){
	return gsl_ran_poisson_pdf(static_cast<unsigned int>(f_dcount),1)
			+ gsl_cdf_poisson_Q(static_cast<unsigned int>(f_dcount), 1);
}
inline
double dolier_poisson_mu1_bonf_Q(const usize_t f_dcount, const usize_t f_total, const usize_t b_dcount, const usize_t b_total, const usize_t nof_kmers){
	return 		(gsl_ran_poisson_pdf(static_cast<unsigned int>(f_dcount),1)
				+ gsl_cdf_poisson_Q(static_cast<unsigned int>(f_dcount), 1))
			* nof_kmers;
}





const double SQRT_PI  = sqrt(3.141592653589793238462);

inline
double dolier_normal(double x, double sigma, double mu){
	return (1/(sigma * SQRT_PI))*(exp(-(((x-mu)*(x-mu)) / (2*sigma*sigma))));
}


inline
double dolier_approx_binomial(const usize_t f_dcount, const usize_t f_total, const usize_t b_dcount, const usize_t b_total){
	double  p = (double)b_dcount / (double)b_total;
	double result = 0;

	double np = static_cast<double>(f_total) * p;
	double npp = np * (1-p);

	for(usize_t i = f_dcount; i<=f_total; i++){
		result += dolier_normal(static_cast<double>(i), np, npp);
	}
	return result;
}




inline
double dolier_binomial(unsigned int x, double p, unsigned int total){
	return gsl_ran_binomial_pdf(x, p, total);
}

inline
double dolier_normall(unsigned int x, double p, unsigned int total){
	double np = static_cast<double>(total) * p;
	double npp = np * (1-p);
	return dolier_normal(static_cast<double>(x), np, npp);
}



}



#endif /* TEST_SUITE_H_ */
