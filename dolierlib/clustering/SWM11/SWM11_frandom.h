/*
 * SWM11_frandom.h
 *
 *  Created on: Feb 5, 2014
 *      Author: vbonnici
 */

#ifndef SWM11_FRANDOM_H_
#define SWM11_FRANDOM_H_




using namespace dolierlib::SWM11;

namespace dolierlib{
namespace SWM11{

class SWM11_frandom{
private:

	SWM11_frandom(){
		srand( (unsigned)time(NULL) );
	};
	~SWM11_frandom(){
	}

public:
	static SWM11_frandom& instance()
	   {
	      static SWM11_frandom INSTANCE;
	      return INSTANCE;
	   }

	static double random(){
		return rand()/(double(RAND_MAX)+1);
	}
	static bool prob(double f){
		return rand() < f;
	}
	static int random(int low, int high){
		return low + ( random() * static_cast<double>(high - low) );
	}
};



}
}



#endif /* SWM11_FRANDOM_H_ */
