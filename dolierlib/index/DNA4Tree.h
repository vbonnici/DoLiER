/*
 * DNA4Tree.h
 *
 *  Created on: Nov 14, 2013
 *      Author: vbonnici
 */

#ifndef DNA4TREE_H_
#define DNA4TREE_H_


#include <stdint.h>

#include "data_ts.h"
#include "DNA5Alphabet.h"

namespace dolierlib{


/*
 * An implementation of a prefix tree to store kmers over the 2bit alphabet {A,C,G,T}.
 * Each node has 4 childs.
 */


class DNA4Tree{

public:

	class DNA4TreeNode{
	public:
		DNA4TreeNode **childs;
		DNA4TreeNode(){
			childs = new DNA4TreeNode*[4];

			childs[0] = NULL;
			childs[1] = NULL;
			childs[2] = NULL;
			childs[3] = NULL;
		}
		~DNA4TreeNode(){
			delete [] childs;
		}
	};


	DNA4TreeNode *root;

	DNA4Tree(){
		root = new DNA4TreeNode();
	}

	bool insert(dna5_t *kmer, int k){
		DNA4TreeNode *n = root;
		for(int i=0; i<k-1; i++){
			if(n->childs[kmer[i]] == NULL){
				n->childs[kmer[i]] =  new DNA4TreeNode();
			}
			n = n->childs[kmer[i]];
		}
		if(n->childs[kmer[k-1]] == NULL){
			n->childs[kmer[k-1]] =  new DNA4TreeNode();
			return false;
		}
		return true;
	}

	bool exists(dna5_t *kmer, int k){
		DNA4TreeNode *n = root;
		for(int i=0; i<k; i++){
			n = n->childs[kmer[i]];
			if(n == NULL)
				return false;
		}
		return true;
	}
};


}


#endif /* DNA4TREE_H_ */
