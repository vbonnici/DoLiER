/*
 * SWM11_SQcomparator.h
 *
 *  Created on: Feb 4, 2014
 *      Author: vbonnici
 */

#ifndef SWM11_SQCOMPARATOR_H_
#define SWM11_SQCOMPARATOR_H_


#include "SWM11_Point.h"
#include "SWM11_point_comparator.h"

using namespace dolierlib::SWM11;

namespace dolierlib{
namespace SWM11{




class SWM11_SQcomparator : public SWM11_point_comparator{
public:
	SWM11_SQcomparator(){}

	virtual double distance(double a[], double b[], int dim){
		double dist = 0;
		for(int i=0; i<dim; i++){
			dist += (a[i]-b[i]) * (a[i]-b[i]);
		}
		return exp(dist);
	}

	virtual double distance(SWM11_Point *a, SWM11_Point *b){
		double dist = 0;

		int dim = 0;
		a->dimension >= b->dimension ? dim = b->dimension : dim = a->dimension;


		int i;
		double it;
		for(i=0; i<dim; i++){
			it = a->coordinates[ i ] - b->coordinates[i];
			it = it * it;
			dist += it;
		}
		for(i=dim; i<a->dimension; i++)
			it += a->coordinates[i] * a->coordinates[i];
		for(i=dim; i<b->dimension; i++)
			it += b->coordinates[i] * b->coordinates[i];


		std::cout<<"dist("<<a->id<<","<<b->id<<") = "<<dist<<" => "<<exp(dist)<<"\n";

		//return exp(dist);
		return dist;
	}

	virtual bool equals(SWM11_Point *a, SWM11_Point *b){
		if(a->dimension != b->dimension)
			return false;
		for(int i=0; i<a->dimension; i++){
			if(a->coordinates[i] != b->coordinates[i])
				return false;
		}
		return true;
	}
};

}
}



#endif /* SWM11_SQCOMPARATOR_H_ */
