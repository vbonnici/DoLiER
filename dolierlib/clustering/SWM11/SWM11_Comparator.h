/*
 * SMW11_comparator.h
 *
 *  Created on: Feb 5, 2014
 *      Author: vbonnici
 */

#ifndef SMW11_COMPARATOR_H_
#define SMW11_COMPARATOR_H_


#include "SWM11_Point.h"

using namespace dolierlib::SWM11;

namespace dolierlib{
namespace SWM11{



enum ComparatorType{SQ};


//template <ComparatorType T>
//class C{
//
//};
//
//template<>
//class C<SQ>{
//
//};



template <ComparatorType T>
class SWM11_comparator{
public:
	double distance(SWM11_Point *a, SWM11_Point *b){
		return 0;
	}

	bool equals(SWM11_Point *a, SWM11_Point *b){
		return true;
	}
//	bool equals(SWM11_Point *a, SWM11_Point *b){
//		return distance(a,b) == 0;
//	}


	double distSQToClosest(SWM11_Point * p, SWM11_Point * ctrs, int numCtrs)
	{
			return distSQToClosest(p, ctrs, 0, numCtrs);
	}


	double distSQToClosest(SWM11_Point * p, SWM11_Point * ctrs, int s, int e)
	{
		double closest = distance( p, ctrs + s);
		double dis;
		int i;
		for (i = s+1; i < e; i++)
		{
			dis = distance(p, ctrs+i);
			if ( dis < closest )
				closest = dis;
		}
		return closest;
	}

	int findClosest(SWM11_Point * data, int n, SWM11_Point * a)
	{
		int indexOfClosest = 0;
		double distToClosest = distance(a, data+0);
		int i;
		double dis;
		for(i=1; i < n; i++)
		{
			dis = distance(a, data+i);
			if ( dis < distToClosest)
			{
				distToClosest = dis;
				indexOfClosest = i;
			}
		}
		return indexOfClosest;
	}

	bool contains(SWM11_Point * set, int size, SWM11_Point * p)
	{
		int i;
		for (i=0; i < size; i++)
		{
			if ( equals(set + i, p) )
				return true;
		}
		return false;
	}

};


template <>
class SWM11_comparator<SQ>{
public:
	double distance(SWM11_Point *a, SWM11_Point *b){
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

		return exp(dist);
	}

	bool equals(SWM11_Point *a, SWM11_Point *b){
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

#endif /* SMW11_COMPARATOR_H_ */
