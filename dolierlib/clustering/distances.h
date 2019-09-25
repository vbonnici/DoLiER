/*
 * distances.h
 *
 *  Created on: Feb 7, 2014
 *      Author: vbonnici
 */

#ifndef DISTANCES_H_
#define DISTANCES_H_


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include "data_ts.h"
#include "trim.h"
#include "DNA5Alphabet.h"
#include "CstyleIndex.h"



namespace dolierlib{



//class VectorDistance{
//public:
//	virtual ~VectorDistance(){};
//
//	virtual double distance(double *a, double *b, double *weights, size_t dim) =0;
//};



double dist_tanimoto(double *a, double *b, double *weights, size_t dim){
	double min = 0;
	double max = 0;
	for(size_t i=0; i<dim; i++){
		if(a[i] < b[i]){
			min += a[i] * weights[i];
			max += b[i] * weights[i];
		}else if(a[i] > b[i]){
			min += b[i] * weights[i];
			max += a[i] * weights[i];
		}
		else{
			min += a[i] * weights[i];
			max += b[i] * weights[i];
		}
	}
	return min/max;
}


double dist_euclidean(double *a, double *b, size_t dim){
	double dist = 0;
	for(size_t i=0; i<dim; i++){
		dist += (a[i]-b[i]) * (a[i]-b[i]);
	}
	return sqrt(dist);
}




}



#endif /* DISTANCES_H_ */
