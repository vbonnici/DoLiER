/*
 * pwm-supp.h
 *
 *  Created on: Dec 9, 2013
 *      Author: vbonnici
 */


/*
 * a collection of procedures to deal with PWM
 */


#ifndef PWM_SUPP_H_
#define PWM_SUPP_H_




#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <math.h>

#include <vector>
#include <limits>

#include "data_ts.h"


namespace dolierlib{



/*
 * load PWM from file and return it as a double matrix.
 * The matrix is a  gt_k(rows) * 4(columns) matrix.
 */
double**  get_gtpwm(
		std::string &gt_pwm_file,//input file
		int* gt_k //PWM length, output pointer
		){
	std::string line;
	std::string field;
	int k = 0;

	std::ifstream ifs;
	ifs.open(gt_pwm_file.c_str());
	if(!ifs.good()){
		std::cout<<"ERROR on opening: "<<gt_pwm_file<<"\n";
		exit(-1);
	}
	while(getline(ifs,line)) k++;

	k--;
	ifs.close();

	*gt_k = k;

	double **matrix = new double*[k];
	for(int i=0; i<k; i++){
		matrix[i] = new double[4];
		memset(matrix[i],0,4*sizeof(double));
	}


	ifs.open(gt_pwm_file.c_str());
	if(!ifs.good()){
		std::cout<<"ERROR on opening: "<<gt_pwm_file<<"\n";
		exit(-1);
	}
	getline(ifs, line);//skip first line
	for(int i=0; i<k; i++){
		getline(ifs, line);

		std::stringstream sss(line.c_str());
		for(int j=0; j<4; j++){
			sss>>field;
			matrix[i][j] = atof(field.c_str());
		}
	}
	ifs.close();

	return  matrix;
}


void print_matrix(std::ostream &os, double **matrix, int rows){
	for(int i=0; i<rows; i++){
		for(int j=0; j<4; j++){
			os<<matrix[i][j]<<"\t";
		}
		os<<"\n";
	}
}

void print_matrix(double **matrix, int rows){
	for(int i=0; i<rows; i++){
		for(int j=0; j<4; j++){
			std::cout<<matrix[i][j]<<"\t";
		}
		std::cout<<"\n";
	}
}

/*
 * normalize matrix by rows.
 */
void normalize(double **matrix, int rows){
	double max = matrix[0][0];
	double min = matrix[0][1];
	for(int i=0; i<rows; i++){
		for(int j=0; j<4; j++){
			if(matrix[i][j] > max)
				max = matrix[i][j];
			if(matrix[i][j] < min)
				min = matrix[i][j];
		}
	}
	for(int i=0; i<rows; i++){
		for(int j=0; j<4; j++){
			matrix[i][j] = (matrix[i][j] - min) / (max - min);
		}
	}
}


/*
 * mormalize matrix by pseudocounts
 */
void pseudocount(double **matrix, int rows, double eps){
	bool zero = false;
	for(int i=0; i<rows; i++){
		for(int j=0; j<4; j++){
			if(matrix[i][j] == 0){
				zero = true;
				break;
			}
		}
	}
	if(zero){
		for(int i=0; i<rows; i++){
			double sum = 0;
			for(int j=0; j<4; j++){
				matrix[i][j] += eps;
				sum += matrix[i][j];
			}
			for(int j=0; j<4; j++){
				matrix[i][j] = matrix[i][j] / sum;
			}
		}
	}
}


//double best_score(double **matrix, int rows, int k){
//	double bscore =   0;
//	double *max = new double[rows];
//	memset(max, 0, rows * sizeof(double));
//	for(int i=0; i<rows; i++){
//		for(int j=0; j<4; j++){
//			if(matrix[i][j] > max[i])
//				max[i] = matrix[i][j];
//		}
//	}
//	for(int i=0; i <= rows-k; i++){
//		double score = 0;
//		for(int j=0; j<k; j++)
//			score += max[i+j];
//		if(score > bscore)
//			bscore = score;
//	}
//	delete [] max;
//	return bscore;
//}
//
//double best_sum_score(double **matrix, int rows, int k){
//	double bscore =   0;
//	double *max = new double[rows];
//	memset(max, 0, rows * sizeof(double));
//	for(int i=0; i<rows; i++){
//		for(int j=0; j<4; j++){
//			if(matrix[i][j] > max[i])
//				max[i] = matrix[i][j];
//		}
//	}
//	for(int i=0; i <= (rows>=k ? rows-k : 0); i++){
//		double score = 0;
//		for(int j=0; j<  (rows>=k ? k : rows); j++)
//			score += max[i+j];
//		if(score > bscore)
//			bscore = score;
//	}
//	delete [] max;
//	return bscore;
//}


/*
 * get the best theoretical score of length k that can be obtained over this matrix.
 */
double best_prob_score(double **matrix, int rows, int k){
	double bscore =   0;
	double *max = new double[rows];
	memset(max, 0, rows * sizeof(double));
	for(int i=0; i<rows; i++){
		for(int j=0; j<4; j++){
			if(matrix[i][j] > max[i])
				max[i] = matrix[i][j];
		}
	}
	for(int i=0; i <= (rows>=k ? rows-k : 0); i++){
		double score = 1;
		for(int j=0; j<(rows>=k ? k : rows); j++)
			score *= max[i+j];
		if(score > bscore)
			bscore = score;
	}
	delete [] max;
	return bscore;
}

//double sum_score(dna5_t *kmer, int kmer_k, double **matrix, int matrix_k){
//	double sum = 0;
//	if(kmer_k == matrix_k){
//		for(int i=0; i<kmer_k; i++){
//			sum += matrix[i][ kmer[i] ];
//		}
//	}
//	else if(kmer_k < matrix_k){
//		double csum = 0;
//
//		for(int i=0; i<= matrix_k - kmer_k; i++){
//
//			csum = 0;
//			for(int j=0; j<kmer_k; j++){
//				csum += matrix[j+i][ kmer[j] ];
//			}
//			if(csum > sum)
//				sum = csum;
//		}
//
//	}
//	else{
//		double csum = 0;
//		for(int i=0; i <= kmer_k - matrix_k; i++){
//			csum = 0;
//			for(int j=0; j<matrix_k; j++){
//				csum += matrix[j][ kmer[i+j] ];
//			}
//			if(csum > sum)
//				sum = csum;
//		}
//	}
//	return sum;
//}


/*
 * get the score of a kmer (of length kmer_k) over the pwm matrix (of length matrix_k)
 */
double prob_score(dna5_t *kmer, int kmer_k, double **matrix, int matrix_k){
	double sum = 1;
	if(kmer_k == matrix_k){
		//kmer length == matrix length, it's simply the product of
		for(int i=0; i<kmer_k; i++){
			sum *= matrix[i][ kmer[i] ]; //it works because A=0,C=1;G=2,T=3  both for dna5_t and matrix columns
		}
	}
	else if(kmer_k < matrix_k){
		//the kmer is shorter than the matrix. Get the best score shifting the kmer over the matrix.
		sum = 0;
		double csum = 1;

		for(int i=0; i<= matrix_k - kmer_k; i++){

			csum = 1;
			for(int j=0; j<kmer_k; j++){
				csum *= matrix[j+i][ kmer[j] ];
			}
			if(csum > sum)
				sum = csum;
		}

	}
	else{
		//the kmer is longer than the matrix. Get the best score shifting the matrix over the kmer.
		sum = 0;
		double csum = 1;
		for(int i=0; i <= kmer_k - matrix_k; i++){
			csum = 1;
			for(int j=0; j<matrix_k; j++){
				csum *= matrix[j][ kmer[i+j] ];
			}
			if(csum > sum)
				sum = csum;
		}
	}
	return sum;
}


}

#endif /* PWM_SUPP_H_ */
