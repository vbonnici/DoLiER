/*
 * data_ts.h
 *
 *  Created on: Sep 26, 2013
 *      Author: vbonnici
 */


/*
 * core data types definition
 */

#ifndef DATA_TS_H_
#define DATA_TS_H_

#include <stdint.h>
#include "limits.h"

namespace dolierlib{


typedef uint32_t 	_ui32;
typedef int32_t 	_si32;

typedef uint8_t		_ui8;
typedef uint16_t	_ui16;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
typedef char		symbol_t;
typedef _ui32		code_t;
typedef _ui8		dna5_t;
typedef _ui32		usize_t;
typedef _si32		ssize_t;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


typedef uint64_t 	_ui64;

typedef _ui8		dna7_t;
typedef _ui16		code7_t;

typedef _ui16		code4_t;

}



#endif /* DATA_TS_H_ */
