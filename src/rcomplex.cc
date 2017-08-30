#include "../include/rcomplex.h"
#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <stdexcept>

#include <stack>
#include <queue>
#include <utility>
#include <vector>

using namespace iRRAM;




/*
 * Complex number with rational real and rational imag
 *
 */
R_COMPLEX::R_COMPLEX()
{
	re = 0;
	im = 0;
}
R_COMPLEX::R_COMPLEX(RATIONAL r, RATIONAL c)
{
	re = r;
	im = c;
}

R_COMPLEX::R_COMPLEX(RATIONAL r)
{
	re = r;
	im = RATIONAL(0);
}

R_COMPLEX::~R_COMPLEX()
{
}
RATIONAL R_COMPLEX::real() const
{
	return re;
}
RATIONAL R_COMPLEX::imag() const
{
	return im;
}
REAL abs(R_COMPLEX A)
{
	return sqrt(A.real()*A.real() + A.imag()*A.imag());
}
void print(R_COMPLEX R)
{
	cout << "("<<REAL(R.real())<<", "<<REAL(R.imag())<<")";
}