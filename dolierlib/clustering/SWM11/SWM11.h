/*
 * SWM11.h
 *
 *  Created on: Feb 5, 2014
 *      Author: vbonnici
 */

#ifndef SWM11_H_
#define SWM11_H_

#include "SWM11_Point.h"
#include "SWM11_point_comparator.h"
#include "SWM11_Facility.h"
#include "SWM11_batch.h"


#include "data_ts.h"
#include "timer.h"
#include "trim.h"
#include "pars_t.h"
#include "dimers.h"
#include "DNA5Alphabet.h"
#include "FASTAReader.h"
#include "CstyleIndex.h"



using namespace dolierlib::SWM11;
using namespace std;

namespace dolierlib{
namespace SWM11{

const bool RESTART_ON_PT = false;
const double E = 2.718281828;
const double ALPHA = 2.0;
const double COFL = 3 * ALPHA + 2 * ( E / (E - 1) );
const double KOFL = 6 * ALPHA + 1;
const double BETA = 2 * ALPHA * ALPHA * COFL + 2 * ALPHA;

class SWM11_ApproximateFacility{
public:
	int facilityNumber;
	double projection;

	SWM11_ApproximateFacility(){
		facilityNumber = -1;
		projection = 0.0;
	}
};



class SWM11_algorithm{
public:

	double **vectors;
	int nof_vectors;
	int dimension;
	int k;
	int maxFacilities;
	int maxSamples;

	SWM11_Facility **facilities;
	SWM11_ApproximateFacility * approxFacs;
	double *projectedPoint;

	SWM11_point_comparator *comp;

	SWM11_algorithm(double **vectors, int nof_vectors, int dimension, int k, int maxFacilities, int maxSamples, SWM11_point_comparator *comp){
		this->vectors = vectors;
		this->nof_vectors = nof_vectors;
		this->dimension = dimension;
		this->k = k;
		this->maxFacilities = maxFacilities;
		this->maxSamples = maxSamples;
		this->comp = comp;

		facilities = NULL;
		approxFacs = NULL;
		projectedPoint = NULL;
//		facilities = new SWM11_Facility*[maxFacilities];
//		approxFacs = new SWM11_ApproximateFacility[maxFacilities];
//		projectedPoint = new double [ dimension ];
	}
	~SWM11_algorithm(){
		if(approxFacs != NULL)
		delete [] approxFacs;
		if(projectedPoint != NULL)
		delete [] projectedPoint;

		if(facilities != NULL)
		for(int i=0; i<maxFacilities; i++)
			delete facilities[i];
		delete [] facilities;
	}


	void run(std::vector<dna5_t*> f_sequences,std::vector<usize_t> f_lengths){
		//TODO delete previous run
		facilities = new SWM11_Facility*[maxFacilities];
		approxFacs = new SWM11_ApproximateFacility[maxFacilities];
		projectedPoint = new double [ dimension ];

		int pit = 0;
		int i;

		SWM11_Point **points = new SWM11_Point*[nof_vectors];
		for(i=0; i<nof_vectors; i++){
			points[i] = new SWM11_Point(i, vectors[i], dimension);
		}

		for (i=0; i < dimension; i++){
			projectedPoint[i] = SWM11_frandom::rand();
		}
		for(i=0; i<maxFacilities; i++){
			facilities[i] = new SWM11_Facility(dimension, maxSamples);
		}
		int numFacilities = 0;

		double delta, f;
		double logn = log(nof_vectors)/log(2);
		f = 1 / (k*(1+logn));


		std::cout<<"_run_1\n";

		SWM11_Point *x1 = new SWM11_Point(pit, vectors[pit], dimension);
		pit++;
		makeFacility(x1, facilities, approxFacs, numFacilities, dimension, projectedPoint);

		double fakedst;
		int numOldFacs;
		SWM11_Facility *y;
		while(pit < nof_vectors){

			SWM11_Point *x = new SWM11_Point(pit, vectors[pit], dimension);
			pit++;

			y = findApproxNearest(comp, x, facilities, approxFacs, numFacilities, delta, projectedPoint, dimension );

//			y = SWM11_Facility::findNearest(*comp, x, facilities, numFacilities, fakedst);
//			std::cout<<"approx_find(";print_dna5(f_sequences[x->id], f_lengths[x->id]);std::cout<<"->";
//			print_dna5(f_sequences[y->p->id], f_lengths[y->p->id]);std::cout<<"\n";


			if(SWM11_frandom::prob(delta/f)){
				makeFacility(x, facilities, approxFacs, numFacilities, dimension, projectedPoint);
				while( numFacilities == maxFacilities ){
					f = f * BETA;
					numOldFacs = numFacilities;
					numFacilities = 1;
					// recenter fac 0 at CoM:
					facilities[0]->setToCoM();
					// reproject it.
					approxFacs[0].projection = project(facilities[0]->p->coordinates, projectedPoint, dimension );
					approxFacs[0].facilityNumber = 0;


					for(i=1; i < numOldFacs; i++)
					{
						facilities[ i ]->setToCoM();

						y = findApproxNearest(comp,  facilities[i]->p, facilities, approxFacs, numFacilities, delta, projectedPoint, dimension );

//						y = SWM11_Facility::findNearest(*comp, facilities[i]->p, facilities, numFacilities, fakedst);
//						std::cout<<"approx_find(";print_dna5(f_sequences[facilities[i]->p->id], f_lengths[facilities[i]->p->id]);std::cout<<"->";
//						print_dna5(f_sequences[y->p->id], f_lengths[y->p->id]);std::cout<<"\n";


						if( SWM11_frandom::prob(delta/f) )
						{
							// make this a facility
							// Do so by moving the memory to the earlier spot.
//							copyFacility( facilities[numFacilities], facilities[i] );
							facilities[numFacilities]->copy(facilities[i]);
							approxFacs[numFacilities].projection = project(facilities[numFacilities]->p->coordinates, projectedPoint, dimension );
							approxFacs[numFacilities].facilityNumber = numFacilities;
							insert(approxFacs, numFacilities);

							numFacilities++;
						}
						else
						{
							facilities[i]->assign(y);
						}
					}
				}
			}
			else{
				y->assign(x);
			}
		}






		int j;
		if( k < numFacilities )
		{
			std::cout<<"OOOOOOOOOOO\n";
			int * finals = SWM11_batch::ORSSSeed(comp, facilities, numFacilities, k);

			SWM11_batch::LocalSearch(comp, facilities, numFacilities, k, finals);

			for (i=0; i < min(numFacilities, k); i++)
			{
				// swap facilities[i] with facilities[ finals[i] ]
				y = facilities[ finals[i] ];
				facilities[ finals[i] ] = facilities[i];
				facilities[i] = y;
			}

			double * dists = new double[k];

			// Assign all other facilities to their nearest of the final k.
			// Then facilities[0..k-1] can be used as the final k-means.
			for (i=0; i < k; i++)
			{
				j = SWM11_Facility::findNN(*comp, facilities, k, i);
				dists[i] = comp->distance( facilities[i]->p, facilities[j]->p );
				SWM11_Facility::ballKMeans(*comp, facilities[i], dists[i] / 9 );
			}
		}

	}
private:

	void makeFacility(
			SWM11_Point * x,
			SWM11_Facility *fs[],
			SWM11_ApproximateFacility *afs,
			int &numFacilities,
			int dim,
			double *proj)
	{
		// remember to copy x, because we will change it later.
		fs[ numFacilities ]->p->copy(x);

		fs[ numFacilities ]->p->weight = 1.0;
		int i;
		for (i=0; i < x->dimension; i++)
		{
			fs[ numFacilities]->totals[i] = x->coordinates[i];
		}

		afs[numFacilities].facilityNumber = numFacilities;
		afs[numFacilities].projection = project(x->coordinates, proj, dim );

		insert( afs, numFacilities );

		numFacilities++;
	}

	double project(double x[], double y[], int d)
	{
		double sum = 0.0;
		int i;
		for (i=0; i < d; i++)
		{
			sum += ( x[i] * y[i] );
		}
		return sum;
	}

	void insert(SWM11_ApproximateFacility * afs, int n)
	{
		int i = n;
		int tempint;
		double tempdbl;
		while( i > 0 && afs[i-1].projection > afs[i].projection )
		{
			// swap i-1, i within afs.
			tempint = afs[i-1].facilityNumber;
			afs[i-1].facilityNumber = afs[i].facilityNumber;
			afs[i].facilityNumber = tempint;

			tempdbl = afs[i-1].projection;
			afs[i-1].projection = afs[i].projection;
			afs[i].projection = tempdbl;

			i--;
		}
	}



	SWM11_Facility * findApproxNearest(
			SWM11_point_comparator *comp,
			SWM11_Point * p,
			SWM11_Facility *fs[],
			SWM11_ApproximateFacility * afs,
			int n, double & distSQ, double * proj, int dim)
	{
		double projection = project(p->coordinates, proj, dim);
		int loc = binarySearch(afs, n, projection);

		if( loc == -1 )
		{
			// nearest is fs[ afs[0].facilityNumber ].
			distSQ = comp->distance( p, fs[ afs[0].facilityNumber]->p );
			return fs[ afs[0].facilityNumber ];
		}
		else if ( loc == n-1)
		{
			// nearest is fs[ afs[n-1].facilityNumber ].
			distSQ = comp->distance( p, fs[ afs[n-1].facilityNumber]->p );
			return fs[ afs[n-1].facilityNumber ];
		}
	//	else:
		// either fs[ afs[loc].facilityNumber ]
		// or 	  fs[ afs[loc+1].facilityNumber ]
		distSQ = comp->distance( p, fs[ afs[loc].facilityNumber]->p);
		double dis = comp->distance( p, fs[ afs[loc+1].facilityNumber]->p );
		if( distSQ <= dis )
		{
			return fs[afs[loc].facilityNumber];
		}
		else
		{
			distSQ = dis;
			return fs[ afs[loc+1].facilityNumber ];
		}
	}


	int binarySearch(SWM11_ApproximateFacility * afs, int n, double target)
	{
		if( target < afs[0].projection )
			return -1;
		if ( target > afs[ n-1 ].projection )
			return n-1;
		int low = 0;
		int high = n-1;
		int mid;
		// Invariant:  nums[low] <= target <= nums[high]
		while( high - low > 1 )
		{
			mid = (high + low) / 2;
			if ( target <= afs[mid ].projection )
			{
				high = mid;
			}
			if ( afs[mid].projection <= target )
			{
				low = mid;
			}
		}
		return low;
	}
};



}
}


#endif /* SWM11_H_ */
