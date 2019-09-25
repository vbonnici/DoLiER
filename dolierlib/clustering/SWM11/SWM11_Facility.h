/*
 * SWM11_Facility.h
 *
 *  Created on: Feb 5, 2014
 *      Author: vbonnici
 */

#ifndef SWM11_FACILITY_H_
#define SWM11_FACILITY_H_

#include "SWM11_Point.h"
#include "SWM11_frandom.h"
#include "SWM11_point_comparator.h"

#include <vector>

using namespace dolierlib::SWM11;

namespace dolierlib{
namespace SWM11{



class SWM11_Facility{
public:
	SWM11_Point *p;
	double *totals;

	std::vector<SWM11_Point*> samples;
	//SWM11_Point **samples;
	//int numSamples;
	int maxSamples;


	SWM11_Facility(int dim, int _maxSamples){
		p = new SWM11_Point(dim);
		totals = (double*)calloc(dim, sizeof(double));
		maxSamples = _maxSamples;
	}


//	SWM11_Facility(const SWM11_Facility &f){
//		p = f.p;
//		totals = f.totals;
//		samples = f.samples;
//	}

	~SWM11_Facility(){
		free(totals);
	}


	void
	copy(SWM11_Facility *f){
		samples.clear();
		p->copy(f->p);
		int i;
		for(i=0; i<p->dimension; i++){
			totals[i] = f->totals[i];
		}
		for(std::vector<SWM11_Point*>::iterator IT=f->samples.begin(); IT!=f->samples.end(); IT++){
			samples.push_back(new SWM11_Point(*(*IT)));
		}
	}


	void assign(SWM11_Point *p){
		p->weight++;
		int i;
		for (i=0; i < p->dimension; i++){
			totals[i] += p->coordinates[i];
		}

		if(samples.size() < maxSamples){
			//copy
			samples.push_back(p);
		}
		else{
			for(i=0; i<samples.size(); i++){
				if(SWM11_frandom::prob(1.0/(maxSamples+1))){
					samples[i] = p;
					break;
				}
			}
		}

	}


	int numToKeepFromA(int n1, int n2){
		int num = 0;
		int sum = n1+n2;
		int i;
		for (i=0; i < maxSamples; i++){
			if( SWM11_frandom::prob( 1.0* n1 / sum ) ){
				num++;
				n1--;
				sum--;
			}
			else{
				n2--;
				sum--;
			}
		}
		return num;
	}

	void mergeSamples(SWM11_Facility *b){
	  int i, n1, n2;
	  n1 = p->weight -1;
	  n2 = b->p->weight -1;



	  // Note:  If n1 >= MAX_SAMPLES, samples full.
	  // 			Not full otherwise.
	  // 			Same with n2.


	  if( n1 + n2 <= maxSamples ){
	      // add all samples from b to a.
	      for(i=0; i < n2; i++){
		  // move b->samples[i] to a->samples[i+ n1]
		  //copy( a->samples[i+n1], b->samples[i] );
	    	  samples.push_back(b->samples[i]);
		}
	      //a->numSamples = n1 + n2;
	      // update a's total.
	      return;	// done
	    }

	  // else, we have at least MAX_SAMPLES total.
	  // Need to determine how many from a to keep
	  // (and, similarly, how many from b to keep)
	  int keepFromA = numToKeepFromA(n1, n2);
	  int keepFromB = maxSamples - keepFromA;
	  for(i=0; i < keepFromB; i++)
	    {
	      // move b->samples[i] to a->samples[i + keepFromA].
	      //copy( a->samples[i + keepFromA], b->samples[i] );
		  if(i+keepFromA < samples.size()){
			  samples[i + keepFromA] = b->samples[i];
		  }
		  else{
			  samples.push_back(b->samples[i]);
		  }
		  //samples[i + keepFromA] = b->samples[i];
	    }

	  // NOTE:  Do *NOT* delete anything in b.
	  // 		  That is reserved for main (or equivalent)

	}


	void assign(SWM11_Facility * b){
		mergeSamples(b);
		p->weight += b->p->weight;
		int i;
		for (i=0; i < p->dimension; i++){
			totals[i] += b->totals[i];
		}
	}

	void setToCoM(){
		int i;
		for (i=0; i < p->dimension;  i++){
			p->coordinates[i] = totals[i] / p->weight;
		}
	}







//	template<ComparatorType T>
	static
	SWM11_Facility*
	findNearest(SWM11_point_comparator &comp, SWM11_Point * p, SWM11_Facility *fs[], int n, double & distSQ)
	{
		int closest = 0;
		int i;
		double dist;
		distSQ = comp.distance(p, fs[0]->p);
		for (i=1; i < n; i++)
		{
			dist = comp.distance(p, fs[i]->p);
			if( dist < distSQ)
			{
				distSQ = dist;
				closest = i;
			}
		}
		return fs[closest];
	}


//	template<ComparatorType T>
	static
	int findNN(SWM11_point_comparator &comp, SWM11_Facility * fs[], int n, int x)
	{
		double dis, dis2;
		int j, nn;
		if( x == 0 )
		{
			dis = comp.distance( fs[0]->p, fs[1]->p );
			nn = 1;
		}
		else
		{
			dis = comp.distance( fs[0]->p, fs[x]->p );
			nn = 0;
		}
		for (j=1; j < n; j++)
		{
			if ( j != x )
			{
				dis2 = comp.distance( fs[x]->p, fs[j]->p);
				if( dis2 < dis )
				{
					nn = j;
					dis = dis2;
				}
			}
		}
		return nn;
	}


	static
	void ballKMeans(SWM11_point_comparator &comp, SWM11_Facility * f, double dist)
	{
		// recenter f at the CoM of all sampled points within the given distance.

		int numInBall = 1;
		int dim = f->p->dimension;
		double * totals = (double*) calloc(dim, sizeof(double));
		int i, j;
		for (i=0; i < dim; i++)
		{
			totals[i] = f->p->coordinates[i];
		}

		for (i=0; i < f->samples.size(); i++)
		{
			// if f->samples[i] within dist of f->p
			// add it to the balled coordinates.
			if( comp.distance( f->samples[i], f->p ) <= dist )
			{
				for (j=0; j < dim; j++)
				{
					totals[j] += f->samples[i]->coordinates[j];
				}
				numInBall++;
			}
		}

		for (i=0; i < dim; i++)
		{
			f->p->coordinates[i] = totals[i] / numInBall;
		}

		free(totals);
	}



	static
	void ballKMeansPoints(SWM11_point_comparator &comp, SWM11_Facility * f, double dist)
	{
		// recenter f at the CoM of all sampled points within the given distance.

		int numInBall = 1;
		int dim = f->p->dimension;
		double * totals = (double*) calloc(dim, sizeof(double));
		int i, j;
		for (i=0; i < dim; i++)
		{
			totals[i] = f->p->coordinates[i];
		}

		for (i=0; i < f->samples.size(); i++)
		{
			// if f->samples[i] within dist of f->p
			// add it to the balled coordinates.
			if( comp.distance( f->samples[i], f->p ) <= dist )
			{
				for (j=0; j < dim; j++)
				{
					totals[j] += f->samples[i]->coordinates[j];
				}
				numInBall++;
			}
		}

		for (i=0; i < dim; i++)
		{
			f->p->coordinates[i] = totals[i] / numInBall;
		}

		free(totals);
	}

};


}
}


#endif /* SWM11_FACILITY_H_ */
