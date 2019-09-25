/*
 * fvector.h
 *
 *  Created on: Feb 4, 2014
 *      Author: vbonnici
 */

#ifndef FVECTOR_H_
#define FVECTOR_H_


#include <limits>
#include <vector>
#include "data_ts.h"

#include <gsl/gsl_rng.h>



using namespace dolierlib;

namespace dolierlib{



/*
 * shuffle sequences... remember to swap also their lengths
 */
void shuffle(
		std::vector<dna5_t*> &f_sequences,
		std::vector<usize_t> &f_lengths,
		int ntimes)
{
	srand( (unsigned)time(NULL) );
	size_t f,s;
	//dna5_t *tmp_s; usize_t tmp_l;
	for(int i=0; i<ntimes; i++){
		f = static_cast<size_t>(ceil((rand()/(RAND_MAX +1.0))*(static_cast<double>(f_sequences.size() -1))));
		s = static_cast<size_t>(ceil((rand()/(RAND_MAX +1.0))*(static_cast<double>(f_sequences.size() -1))));
//		std::cout<<f<<"\t"<<s<<"\n";
		std::swap(f_sequences[f], f_sequences[s]);
		std::swap(f_lengths[f], f_lengths[s]);
	}
}



/*
 * swap randomly the first leading_possequences with a forward position
 */
void leading_shuffle(
		std::vector<dna5_t*> &f_sequences,
		std::vector<usize_t> &f_lengths,
		size_t leading_pos)
{
	srand( (unsigned)time(NULL) );


	size_t f;
	for(size_t i=0; i<leading_pos; i++){
		f = static_cast<size_t>(ceil((rand()/(RAND_MAX))*(static_cast<double>(f_sequences.size()- leading_pos -1)))) + leading_pos;
		//f = static_cast<size_t>((gsl_rng_uniform (r) * static_cast<double>(f_sequences.size() - leading_pos))) + leading_pos
//		std::cout<<i<<"\t"<<f<<"\n";

//		print_dna5(f_sequences[i], f_lengths[i]);std::cout<<"\t";print_dna5(f_sequences[f], f_lengths[f]);std::cout<<"\n";

		std::swap(f_sequences[i], f_sequences[f]);
		std::swap(f_lengths[i], f_lengths[f]);

//		print_dna5(f_sequences[i], f_lengths[i]);std::cout<<"\t";print_dna5(f_sequences[f], f_lengths[f]);std::cout<<"\n";
	}

}


/*
 * Get kmers spectra.
 * Only non-zero columns are reported.
 * Let X to be the length of the longest input sequences, then theoretically we have sum_{i=1}^{i<=X}(4^i) features (possible kmers)
 * but some of these theoretical kmers do not appear in the input sequences, so they are excluded from the output spectra.
 * Moreover, if sequences are all of the same length, the working X is X-1.
 *
 */
void fvector(
		std::vector<dna5_t*> &f_sequences,
		std::vector<usize_t> &f_lengths,
		double ***vectors,  //output spectra, it's pointer
		size_t *vector_length, //length of spectra, it's a pointer
		double **weights
		)
{

	DNA5MS_t f_ms(true);
	concat(f_ms, f_sequences, f_lengths, false);
	f_ms.areUndefTeminated = true;

	CsFullyContainer f_ds(true);
		build_fully_ds(f_ms, f_ds);
		//build_strict_ds(f_ms, f_ds);
	Cs5DLIndex f_dlindex(f_ms.seq, f_ms.seq_length, f_ds.SA, f_ms.lengths, f_ms.nof_seqs);
	SASearcher sas(f_ms.seq, f_ds.SA, f_ms.seq_length);


	size_t max_length = 0;
	for(size_t i=0; i<f_lengths.size(); i++){
		if(f_lengths[i] > max_length)
			max_length = f_lengths[i];
	}



	size_t nof_features = 0;
	size_t *nof_kfeatures = new size_t[max_length+1];
	for(size_t i=1; i<=max_length; i++){
		NSAIterator it = NSAIterator::begin(f_ms, f_ds, i);
		nof_kfeatures[i] = nof(it);
		nof_features += nof_kfeatures[i];

		std::cout<<i<<"\t"<<nof_kfeatures[i]<<"\t"<<f_sequences.size()<<"\n";

		//if sequences are all of the same length, the working X is X-1.
		if(nof_kfeatures[i] == f_sequences.size()  &&  i<=max_length){
			if(i>1){
				std::cout<<"#\n";
				max_length = i-1;
				nof_features -= nof_kfeatures[i];
			}
			else{
				max_length = i;
			}
			break;
		}
	}


	//double **table = new double*[f_sequences.size()];
	(*vectors) =  new double*[f_sequences.size()];
	for(size_t i=0; i<f_sequences.size(); i++){
		(*vectors)[i] =  new double[nof_features];
		for(size_t j=0; j<nof_features; j++)
			(*vectors)[i][j] = 0;
	}
	(*weights) =  new double[nof_features];

	size_t c_nf = 0;
	size_t ii;
	for(size_t k=1; k<=max_length; k++){
		NSAIterator it = NSAIterator::begin(f_ms, f_ds, k);
		dna5_t *kmer  = new dna5_t[k];
		while(it.next()){
			it.get_kmer(kmer);
//			std::cout<<c_nf<<"\t";print_dna5(kmer,k); std::cout<<"\n";

			for(ii=it.i_start; ii<it.i_end; ii++){
				(*vectors)[f_dlindex.ss_ids[ii]][c_nf]++;
				(*weights)[c_nf] = k;
			}
			c_nf++;
		}
	}

	*vector_length = nof_features;
}



/*
 * normalize vectors by column, from [0,max] to [min,max]/(max-min)
 */
void normalize(double **matrix, size_t size1, size_t size2)
{
	for(size_t i=0; i<size2; i++){
		double min = std::numeric_limits<double>::max();
		double max = std::numeric_limits<double>::min();

		for(size_t j=0; j<size1; j++){
			if(matrix[j][i]< min)
				min = matrix[j][i];
			if(matrix[j][i] > max)
				max = matrix[j][i];
		}

		double f = max - min;
		for(size_t j=0; j<size1; j++){
			matrix[j][i] = (matrix[j][i] - min) / f;
		}
	}
}


/*
 * normalize vectors by column, from [0,max] to [0,max]/(max)
 */
void normalize_by_sum(double **matrix, size_t size1, size_t size2)
{
	for(size_t i=0; i<size2; i++){
		double sum = 0;
		for(size_t j=0; j<size1; j++){
			sum += matrix[j][i];
		}
		for(size_t j=0; j<size1; j++){
			matrix[j][i] /= sum;
		}
	}
}



/*
 * Keep only those columns having at least one element >= min_value.
 * Unselected columns are putted at the end of the matrix.
 */
void keep_only(
		double **matrix, //input matrix
		size_t size1, //number of rows
		size_t size2, //number of columns
		double min_value, //min value for selection
		size_t *new_size2 //result number of columns
		)
{
	size_t to_size = size2;
	size_t to_pos = 0;
	for(size_t j=0; j<size2; j++){
		bool keep = false;
		for(size_t i=0; i<size1; i++){
			if(matrix[i][j] >= min_value){
				keep = true;
				break;
			}
		}

		if(keep){
			if(j != to_pos){
				for(size_t i=0; i<size1; i++){
					matrix[i][to_pos] = matrix[i][j];
				}
			}
			to_pos++;
		}
		else{
			to_size--;
		}
	}
	*new_size2 = to_size;
}


/*
 * Same behaviour of keep_only, but it swap the weights matrix, too.
 */
void keep_only(double **matrix, size_t size1, size_t size2, double min_value, size_t *new_size2, double *weights)
{
	size_t to_size = size2;
	size_t to_pos = 0;
	for(size_t j=0; j<size2; j++){
		bool keep = false;
		for(size_t i=0; i<size1; i++){
			if(matrix[i][j] >= min_value){
				keep = true;
				break;
			}
		}

		if(keep){
			if(j != to_pos){
				for(size_t i=0; i<size1; i++){
					matrix[i][to_pos] = matrix[i][j];
				}
				weights[to_pos] = weights[j];
			}
			to_pos++;
		}
		else{
			to_size--;
		}
	}
	*new_size2 = to_size;
}




/*
 * Theoretical spectra are vectors if length sum_{i=1}^{i<=X}(4^i).
 * Thus, weights[i] correspond to the specific length i of the releated theoretical kmer.
 *
 */
void get_fvectors_weights(double *weights, size_t length, std::vector<dna5_t*> &f_sequences, std::vector<usize_t> &f_lengths){
//	size_t base = 0;
//	size_t ebase = 0;
//	size_t cbase = 0;
//	for(size_t i=0;i<length; i++){
//		if(cbase == ebase){
//			base++;
//			ebase = static_cast<size_t>(pow(4, base));
//			cbase = 0;
//		}
//		weights[i] = base;
//		cbase++;
//	}

	DNA5MS_t f_ms(true);
	concat(f_ms, f_sequences, f_lengths, false);
	f_ms.areUndefTeminated = true;

	CsFullyContainer f_ds(true);
		build_fully_ds(f_ms, f_ds);
		//build_strict_ds(f_ms, f_ds);
	Cs5DLIndex f_dlindex(f_ms.seq, f_ms.seq_length, f_ds.SA, f_ms.lengths, f_ms.nof_seqs);
	SASearcher sas(f_ms.seq, f_ds.SA, f_ms.seq_length);


	size_t max_length = 0;
	for(size_t i=0; i<f_lengths.size(); i++){
		if(f_lengths[i] > max_length)
			max_length = f_lengths[i];
	}



	size_t nof_features = 0;
	size_t *nof_kfeatures = new size_t[max_length+1];
	for(size_t i=1; i<=max_length; i++){
		NSAIterator it = NSAIterator::begin(f_ms, f_ds, i);
		nof_kfeatures[i] = nof(it);
		nof_features += nof_kfeatures[i];

		std::cout<<i<<"\t"<<nof_kfeatures[i]<<"\t"<<f_sequences.size()<<"\n";

		//if sequences are all of the same length, the working X is X-1.
		if(nof_kfeatures[i] == f_sequences.size()  &&  i<=max_length){
			if(i>1){
				std::cout<<"#\n";
				max_length = i-1;
				nof_features -= nof_kfeatures[i];
			}
			else{
				max_length = i;
			}
			break;
		}
	}


	size_t c_nf = 0;
	size_t ii;
	for(size_t k=1; k<=max_length; k++){
		NSAIterator it = NSAIterator::begin(f_ms, f_ds, k);
		dna5_t *kmer  = new dna5_t[k];
		while(it.next()){
			it.get_kmer(kmer);
//			std::cout<<c_nf<<"\t";print_dna5(kmer,k); std::cout<<"\n";
			for(ii=it.i_start; ii<it.i_end; ii++){
				weights[c_nf] = k;
			}
			c_nf++;
		}
	}
}



void print_matrix(double **matrix, size_t size1, size_t size2){
	for(size_t i=0; i<size1; i++){
		for(size_t j=0; j<size2; j++){
			std::cout<<matrix[i][j]<<" ";
		}
		std::cout<<"\n";
	}
}



}


#endif /* FVECTOR_H_ */
