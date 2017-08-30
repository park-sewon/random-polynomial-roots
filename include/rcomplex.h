/*
the header file defines a countable subset of complex number which is a pair of rationals
*/
#ifndef RCOMPLEX_H
#define RCOMPLEX_H

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <utility>
using namespace iRRAM;

class R_COMPLEX
{
	private:
		RATIONAL re;
		RATIONAL im;
	public:
		R_COMPLEX();
		R_COMPLEX(RATIONAL, RATIONAL);
		R_COMPLEX(RATIONAL);
		~R_COMPLEX();
		RATIONAL real() const;
		RATIONAL imag() const;

		friend R_COMPLEX operator +(const R_COMPLEX& lhs, const R_COMPLEX& rhs)
		{
			return R_COMPLEX(lhs.real()+rhs.real(), lhs.imag()+rhs.imag());
		}
		friend R_COMPLEX operator -(const R_COMPLEX& lhs, const R_COMPLEX& rhs)
		{
			return R_COMPLEX(lhs.real()-rhs.real(), lhs.imag()-rhs.imag());
		}
};
REAL abs(R_COMPLEX );
void print(R_COMPLEX);

#endif
