/*
 * SWM11_point_comparator.h
 *
 *  Created on: Feb 4, 2014
 *      Author: vbonnici
 */

#ifndef SWM11_POINT_COMPARATOR_H_
#define SWM11_POINT_COMPARATOR_H_



#include "SWM11_Point.h"

using namespace dolierlib::SWM11;

namespace dolierlib{
namespace SWM11{

class SWM11_point_comparator{
public:
	virtual ~SWM11_point_comparator(){};


	virtual double distance(double a[], double b[], int dim) =0;

	virtual double distance(SWM11_Point *a, SWM11_Point *b) =0;

	virtual bool equals(SWM11_Point *a, SWM11_Point *b) =0;
//	bool equals(SWM11_point *a, SWM11_point *b){
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

}
}



#endif /* SWM11_POINT_COMPARATOR_H_ */
