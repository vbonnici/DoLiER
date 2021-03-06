/*
 * seq1.h
 *
 *  Created on: Feb 5, 2014
 *      Author: vbonnici
 */



/*
 * An inefficient version of seq1
 */


#ifndef SEQ2_H_
#define SEQ2_H_


#include <stdlib.h>
#include <vector>
#include <limits>

#include "distances.h"


namespace dolierlib{
namespace seq1{



class seq2_algorithm{
public:

	double **vectors;
	size_t nof_vectors;
	size_t dim;
	double *weights;
	double distThr;

	double maxTDist;
	size_t noCluster;


	std::vector< std::set<size_t> > clusters;
	std::vector<double*> centroids;



	seq2_algorithm(double **vectors, size_t nof_vectors, size_t dim, double *weights, double distThr){
		this->vectors = vectors;
		this->nof_vectors = nof_vectors;
		this->dim = dim;
		this->weights = weights;
		this->distThr = distThr;

		this->maxTDist = std::numeric_limits<double>::max();
		this->noCluster = std::numeric_limits<size_t>::max();

	}
	~seq2_algorithm(){
		doDelete();
	}

	void doDelete(){
		for(std::vector<double*>::iterator IT = centroids.begin(); IT!=centroids.end(); IT++){
			delete [] (*IT);
		}
	}

	//void run(){
	void run(
			std::vector<dna5_t*> &f_sequences,
			std::vector<usize_t> &f_lengths){
		doDelete();
		clusters.clear();
		centroids.clear();

		if(nof_vectors > 0){
			newCluster(0);
		}

		size_t toCluster;
		for(size_t vit = 1; vit < nof_vectors; vit++){
			toCluster = closestCluster(vit);

//			print_dna5(f_sequences[vit], f_lengths[vit]); std::cout<<"\t";
//			if(toCluster == noCluster) std::cout<<clusters.size()<<"\t#\n";
//			else std::cout<<toCluster<<"\n";

			if(toCluster == noCluster){
				newCluster(vit);
			}
			else{
				addToCluster(vit, toCluster);
			}
		}
	}


	double* copy(const double *src, size_t dim){
		double *dest = new double[dim];
		for(size_t i=0; i<dim; i++)
			dest[i] = src[i];
		return dest;
	}

	void newCluster(size_t from_i){
		std::set<size_t> elements;
		elements.insert(from_i);
		clusters.push_back(elements);
		centroids.push_back(copy(vectors[from_i], dim));
	}
	void addToCluster(size_t ielement, size_t icluster){
		clusters[icluster].insert(ielement);
		//mean(centroids[icluster], vectors[ielement], dim);
		//mean0(centroids[icluster], vectors[ielement], dim);
		for(size_t i=0; i<dim; i++){
			centroids[icluster][i] = 0;
		}
		for(std::set<size_t>::iterator IT=clusters[icluster].begin(); IT!=clusters[icluster].end(); IT++){
			for(size_t i=0; i<dim; i++){
				centroids[icluster][i] += vectors[*IT][i];
			}
		}
		double ss = static_cast<double>(clusters[icluster].size());
		for(size_t i=0; i<dim; i++){
			centroids[icluster][i] /= ss;
		}
	}

	void mean(double *a, const double *b, size_t dim){
		for(size_t i=0; i<dim; i++){
			a[i] = (a[i] + b[i])/2;
		}
	}

//	void mean0(double *a, const double *b, size_t dim){
//		for(size_t i=0; i<dim; i++){
//			if(a[i]!=0 && b[i]!=0)
//				a[i] = (a[i] + b[i])/2;
//			else{
//				if(a[i]==0)
//					a[i] = b[i];
//			}
//		}
//	}

	size_t closestCluster(size_t ielement){
		size_t ret = noCluster;
		double minDist = maxTDist;

		double *ivector = vectors[ielement];
		double dist;

		for(size_t i=0; i<centroids.size(); i++){
//			dist = distance(ivector, centroids[i], weights, dim);
//			if(tanimoto(ivector,centroids[i],dim) < distThr  &&  dist < minDist){
//				minDist = dist;
//				ret = i;
//			}
			dist = dist_tanimoto(ivector,centroids[i], weights, dim);
//			std::cout<<"\t\t"<<i<<"\t"<<dist<<"\n";

			if(dist >= distThr  &&  dist < minDist){
				minDist = dist;
				ret = i;
			}
		}


//		if(ret != noCluster){
//			size_t count = 0;
//			for(size_t i=0; i<centroids.size(); i++){
//				dist = tanimoto(ivector,centroids[i], weights, dim);
//				if(dist == minDist)
//					count++;
//			}
//			if(count > 1)
//				std::cout<<"#\t"<<count<<"\n";
//		}

		return ret;
	}

//	double distance(double *a, double *b, double *weights, size_t dim){
//		double dist = 0;
//
//		for(size_t i=0; i<dim; i++){
//			dist += ((a[i]-b[i]) * (a[i]-b[i])) /weights[i];
//		}
//		return sqrt(dist);
//	}


//	double tanimoto(double *a, double *b, double *weights, size_t dim){
//		double min = 0;
//		double max = 0;
//		for(size_t i=0; i<dim; i++){
//			if(a[i] < b[i]){
//				min += a[i] * weights[i];
//				max += b[i] * weights[i];
//			}else if(a[i] > b[i]){
//				min += b[i] * weights[i];
//				max += a[i] * weights[i];
//			}
//			else{
//				min += a[i] * weights[i];
//				max += b[i] * weights[i];
//			}
//		}
//		return min/max;
//	}

//	double tanimoto(){
//
//	}





	void get_metoids_by_mean(std::vector<size_t> &metoids){
		metoids.clear();

		size_t el;
		double *vec;

		size_t i,j;
		double count;
		size_t choose;
		double minDist, dist;

		double *mean = new double[dim];

		for(i=0; i<clusters.size(); i++){

			choose = *(clusters[i].begin());

			if(clusters[i].size() > 1){
				for(j=0; j<dim; j++)
					mean[j] = 0;

				count = static_cast<double>(clusters[i].size());
				for(std::set<size_t>::iterator IT=clusters[i].begin(); IT!=clusters[i].end(); IT++){
					el = (*IT);
					vec = vectors[el];
					for(j=0; j<dim; j++){
						mean[j] += vec[j];
					}
				}
				for(j=0; j<dim; j++){
					mean[j] /= count;
				}


				minDist = std::numeric_limits<double>::max();
				for(std::set<size_t>::iterator IT=clusters[i].begin(); IT!=clusters[i].end(); IT++){
					el = (*IT);
					vec = vectors[el];
					dist = dist_tanimoto(vec,mean, weights, dim);
					//dist = dist_euclidean(vec,mean, dim);
					if(dist < minDist){
						choose = el;
						minDist = dist;
					}
				}
			}

			metoids.push_back(choose);
		}

		delete [] mean;
	}

	void get_metoids_by_mindist(std::vector<size_t> &metoids){
		metoids.clear();



		size_t c;
		size_t choose;
		double minDist, dist;

		for(c=0; c<clusters.size(); c++){

			choose = *(clusters[c].begin());;

			if(clusters[c].size()>1){
				minDist = std::numeric_limits<double>::max();
				for(std::set<size_t>::iterator iIT=clusters[c].begin(); iIT!=clusters[c].end(); iIT++){
					dist = 0;
					for(std::set<size_t>::iterator jIT=clusters[c].begin(); jIT!=clusters[c].end(); jIT++){
						dist += dist_tanimoto(vectors[*iIT], vectors[*jIT], weights, dim);
						//dist += dist_euclidean(vectors[*iIT], vectors[*jIT], dim);
					}
					if(dist < minDist){
						choose = *iIT;
						minDist = dist;
					}
				}
			}

			metoids.push_back(choose);
		}
	}





};


}
}





#endif /* SEQ2_H_ */
