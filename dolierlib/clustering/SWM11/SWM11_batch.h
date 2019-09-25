/*
 * SWM11_batch.h
 *
 *  Created on: Feb 5, 2014
 *      Author: vbonnici
 */

#ifndef SWM11_BATCH_H_
#define SWM11_BATCH_H_

#include <math.h>
#include "SWM11_Facility.h"
#include "SWM11_Point.h"
#include "SWM11_point_comparator.h"

using namespace dolierlib::SWM11;
using namespace std;

namespace dolierlib{
namespace SWM11{


class SWM11_batch{
public:

static
int * ORSSSeed(SWM11_point_comparator *comp, SWM11_Facility * fs[], int n, int k)
{
	int dim = fs[0]->p->dimension;
	int * ctrs = new int[k];
	int i;
	double * probs  = new double[ n ];
	double * dist_sq = new double[ n ];
	// compute these among first two.
	fillFirstTwo(comp, fs, n, ctrs, dim, dist_sq);
	for (i=2; i < min(n,k); i++)
	{
		// select ith point:
		determineProbabilities(dist_sq, fs, probs, n);
		ctrs[i] = selectWithProb(probs, n);

		updateDist(comp, dist_sq, n, ctrs[i], fs);

	}
	delete [] dist_sq;
	delete [] probs;
	return ctrs;
}



static
void fillFirstTwo(SWM11_point_comparator *comp, SWM11_Facility * fs[], int n, int * ctrs, int dim, double dist_sq[])
{
	double OPT_1 = 0; // optimal 1-means cost.
	double * com = new double[ dim ]; // weighted com of points
	double * probs = new double[ n ];
	double sum = 0.0;

	double total_weight = 0.0;
	int i, j;
	double weight;
	double ctr_to_c1_sq;

	for (i=0; i < n; i++)
	{
		weight = fs[i]->p->weight;
		total_weight += weight;
		for(j=0; j < dim; j++)
		{
			com[j] += ( weight * fs[i]->p->coordinates[j] );
		}
	}
	for (j=0; j < dim; j++)
	{
		com[j] = com[j] / total_weight;
	}
	for (i=0; i < n; i++)
	{
		dist_sq[i] = comp->distance( com, fs[i]->p->coordinates, dim);
		OPT_1 += dist_sq[i] * fs[i]->p->weight;
	}
	for (i=0; i < n; i++)
	{
		probs[i] = (OPT_1 + total_weight* dist_sq[i] )/ (2 * total_weight * OPT_1 );
		sum += probs[i];
	}

	for (i=0; i < n; i++)
	{
		probs[i] = probs[i] / sum;
	}

	// Select first with probability equal to probs[i].
	ctrs[0] = selectWithProb(probs, n);
	if( ctrs[0] == n)
	{
		cerr << "Error in ORSS Seed selection: first point." << endl;
		exit(1);
	}
	ctr_to_c1_sq = comp->distance(com, fs[ ctrs[0] ]->p->coordinates, dim);
	SWM11_Facility * y = fs[ ctrs[0] ]; // first center

	// refill dist_sq and probs for point two:
	sum = 0.0;
	probs[ ctrs[0] ] = 0;
	for (i=0; i < n; i++)
	{
		if( i != ctrs[0] )
		{
			dist_sq[i] = comp->distance( fs[i]->p, y->p ) * fs[i]->p->weight;
			probs[i] = dist_sq[i] / ( OPT_1 + total_weight * ctr_to_c1_sq );
			sum += probs[i];
		}
	}
	for (i=0; i < n; i++)
	{
		probs[i] = probs[i] / sum;
	}

	ctrs[1] = selectWithProb(probs, n);

	dist_sq[ ctrs[1] ] = 0;


	if( ctrs[1] == n)
	{
		cerr << "Error in ORSS Seed selection: second point." << endl;
		exit(1);
	}

	for (i=0; i < n; i++)
	{
		dist_sq[i] = min( comp->distance(fs[i]->p, fs[ctrs[0]]->p),
						  comp->distance(fs[i]->p, fs[ctrs[1]]->p ) );
	}

	delete [] com;
	delete [] probs;
}



static
int selectWithProb(double probs[], int n)
{
	double cutoff = SWM11_frandom::rand(); // is [0,1]
	double cur_pt = 0.0;
	int i;
	for (i=0; i < n; i++)
	{
		cur_pt += probs[i];
		if ( cur_pt >= cutoff )
		{
			return i;
		}
	}

	cout << "Going to return n = " << n << endl;
	cout << "cutoff = " << cutoff << endl;
	double sum = 0.0;
	for (i=0; i < n; i++)
	{
		sum += probs[i];
		cout << " prob[ " << i << " ] = " << probs[i] << endl;
	}
	cout << "sum = " << sum << endl;

	return n;
}


static
void updateDist(SWM11_point_comparator *comp, double dist_sq[], int n, int chosen, SWM11_Facility * fs[])
{
	int i;
	for (i=0; i < n; i++)
	{
		dist_sq[i] = min( dist_sq[i], comp->distance( fs[i]->p, fs[ chosen ]->p ));
	}


}

static
void determineProbabilities(double dist_sq[], SWM11_Facility * fs[], double probs[], int n)
{
	double sum = 0;
	int i;
	for (i=0; i < n; i++)
	{
		probs[i] = dist_sq[i] * fs[i]->p->weight;
		sum += probs[i];
	}
	//cerr << "\tSum in detProb = " << sum << endl;
	for (i=0; i < n; i++)
	{
		probs[i] = probs[i] / sum;
	}
}

static
void LocalSearch(SWM11_point_comparator *comp, SWM11_Facility * fs[], int n, int k, int start[])
{
	int i, ii, j, was;
	bool changed;
	double bestCost = kMeansCost(comp, fs, n, k, start);
	double cost;
	int * best = new int[k];
	do
	{
		changed = false;
		for (i=0; i < n; i++)
		{
			if( ! contains( start, k, i ) )
			{
				for (j=0; j < k; j++)
				{
					was = start[j];
					start[j] = i;
					cost = kMeansCost(comp, fs, n, k, start);
					if( cost < bestCost )
					{

						// copy start into best
						for (ii = 0; ii < k; ii++)
						{
							best[ii] = start[ii];
						}
						bestCost = cost;
						changed = true;
					}
					start[j] = was;
				}
			}
		}
		if( changed )
		{
			for (i=0; i < k; i++)
			{
				start[i] = best[i];
			}
		}
	}while( changed );

	delete [] best;
}

static
bool contains(int a[], int n, int x)
{
	// is there an i s.t. a[i] == x?
	int i;
	for (i=0; i < n; i++)
	{
		if ( a[i] == x ) return true;
	}

	return false;
}

static
double distNearest(SWM11_point_comparator *comp, SWM11_Facility *fs[], int x, int k, int start[])
{
	double dis, dis2;
	int i;
	dis = comp->distance( fs[x]->p, fs[ start[0] ]->p );
	for (i=1; i < k; i++)
	{
		dis2 = comp->distance( fs[x]->p, fs[ start[i] ]->p );
		if( dis2 < dis)
		{
			dis = dis2;
		}
	}
	return dis;
}

static
double kMeansCost(SWM11_point_comparator *comp, SWM11_Facility * fs[], int n, int k, int start[])
{
	double sum = 0.0;
	int i;
	for (i=0; i < n; i++)
	{
		sum += distNearest(comp, fs, i, k, start);
	}
	return sum;
}


};


}
}


#endif /* SWM11_BATCH_H_ */
