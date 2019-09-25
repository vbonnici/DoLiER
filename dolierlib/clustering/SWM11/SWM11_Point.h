/*
 * SWM11_point.h
 *
 *  Created on: Feb 4, 2014
 *      Author: vbonnici
 */

#ifndef SWM11_POINT_H_
#define SWM11_POINT_H_



namespace dolierlib{
namespace SWM11{


class SWM11_Point{
public:
	int id;
	double *coordinates;
	int dimension;
	double weight;

	SWM11_Point(int _id, double *_coordinates, int _dimension)
	 : id(_id), dimension(_dimension){
		coordinates = (double*)malloc(dimension * sizeof(double));
		for(int i=0; i<dimension; i++){
			coordinates[i] = _coordinates[i];
		}
		weight = 1;
	}
	SWM11_Point(int _dimension){
		id = -1;
		dimension = _dimension;
		coordinates = (double*)calloc(dimension, sizeof(double));
		weight = 1.0;
	}
//	SWM11_Point(){
//		id = -1;
//		coordinates = NULL;
//		dimension = 0;
//		weight = 1.0;
//	}

	SWM11_Point(const SWM11_Point &p){
		id = p.id;
		dimension = p.dimension;
		weight = 1.0;
		coordinates = new double[dimension];
		for(int i=0; i<dimension; i++)
			coordinates[i] = p.coordinates[i];
	}

	void
	copy(SWM11_Point *p){
//		if(coordinates != NULL){
//			free(coordinates);
//		}
		id = p->id;
		dimension = p->dimension;
		//unsafe
		for(int i=0; i<dimension; i++){
			coordinates[i] = p->coordinates[i];
		}
		weight = p->weight;

	}


	~SWM11_Point(){
		if(coordinates != NULL)
			free(coordinates);
	}


};


}
}



#endif /* SWM11_POINT_H_ */
