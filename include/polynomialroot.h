/*
the header file defines programs for approximating roots of a complex polynomial
*/
#ifndef POLYNOMIALROOT_H
#define POLYNOMIALROOT_H

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "polynomial.h"
#include "rcomplex.h"


#include <utility>
using namespace iRRAM;

namespace iRRAM{template <> struct is_continuous<std::vector< std::pair<REAL, int> > > : public std::true_type{};}

namespace iRRAM{
inline sizetype geterror( const std::vector<std::pair<REAL, int> > & l){	sizetype error, lerror;
	l[0].first.geterror(error);
	for (unsigned int i=0; i< l.size(); i++)
	{
		lerror = geterror(l[i].first);
		sizetype_max(error,error,lerror);
	}
	return error;}}


namespace iRRAM{
inline void seterror(std::vector<std::pair<REAL, int> > & l, sizetype & error)
{
	for (unsigned int i=0; i< l.size(); i++)
		l[i].first.seterror(error);
}}

std::vector<R_COMPLEX >
root_approximation_naive(int , POLYNOMIAL );

std::vector<R_COMPLEX >
root_approximation_newton(int , POLYNOMIAL );

#endif
