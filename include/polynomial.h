/*
the header file defines polynomial represented as a list of its complex coefficients
*/
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <utility>
using namespace iRRAM;

class POLYNOMIAL
{
	public:
		std::vector<COMPLEX> coef;
		int degree;

		POLYNOMIAL(int, COMPLEX *);
		POLYNOMIAL(int, std::vector<COMPLEX>);
		POLYNOMIAL(COMPLEX);
		POLYNOMIAL();
		~POLYNOMIAL();
		COMPLEX operator () (const COMPLEX&);
};

POLYNOMIAL operator + (const POLYNOMIAL&, const POLYNOMIAL&);
POLYNOMIAL operator - (const POLYNOMIAL&, const POLYNOMIAL&);
POLYNOMIAL operator * (const POLYNOMIAL&, const POLYNOMIAL&);


/*
POLYNOMIAL operator + (const POLYNOMIAL&, const REAL&);
POLYNOMIAL operator + (const REAL&, const POLYNOMIAL&);
POLYNOMIAL operator - (const POLYNOMIAL&);
POLYNOMIAL operator - (const POLYNOMIAL&, const REAL&);
POLYNOMIAL operator - (const REAL&, const POLYNOMIAL&);
POLYNOMIAL operator * (const POLYNOMIAL&, const REAL&);
POLYNOMIAL operator * (const REAL&, const POLYNOMIAL&);
POLYNOMIAL operator / (const POLYNOMIAL&, const REAL&);
*/

namespace iRRAM{template <> struct is_continuous<POLYNOMIAL> : public std::true_type{};}
// void geterror(const POLYNOMIAL &);
int geterror_exp(const POLYNOMIAL & );

// Q = deriv(P,k) := Q^{(k)}
POLYNOMIAL deriv(POLYNOMIAL, int);
// y = evaluate(P, x) := P(x)
COMPLEX evaluate(POLYNOMIAL, COMPLEX);
// Taylor coefficients: a_k = CoefAt(P, k, z) := P^{(k)}(z)/k!
COMPLEX CoefAt(POLYNOMIAL, int, COMPLEX);
// Translation and dilation by g(x) = f(ax + b)
// a : REAL, b : COMPLEX 
POLYNOMIAL translation(POLYNOMIAL , REAL , COMPLEX );





#endif
