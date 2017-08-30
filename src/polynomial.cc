#include "../include/polynomial.h"
#include "../include/utilities.h"

#include "iRRAM/lib.h"
#include "iRRAM/core.h"

#include <stack>
#include <queue>
#include <utility>
#include <vector>
#include <cmath>
using namespace iRRAM;


// void REALMATRIX::geterror (sizetype& error) const
// {
//   int i;
//   sizetype lerror; 
//   ELEMENT((*this),0,0).geterror(error); 
//   for (i=0;i<(*this).maxcolumn;i++) 
//   for (j=0;j<(*this).maxrow;j++) {
//       ELEMENT((*this),i,j).geterror(lerror);
//       sizetype_max(error,error,lerror);
//   }
// }

// should be sup norm induced metric
int geterror_exp(const POLYNOMIAL & P)
{
	sizetype error; 
  	P.coef[0].geterror(error);
  	int ex = error.exponent;
  	int tmp;

	for(int i=1; i< P.degree; i++)
	{
		P.coef[0].geterror(error);
		tmp = error.exponent;
	  	if (tmp >ex)
	  		ex = tmp;
	}

	return tmp;
}



POLYNOMIAL::POLYNOMIAL(int d, COMPLEX* c)
{
	std::vector<COMPLEX> co;
	co.reserve(d+1);
	for (int i = 0; i<d+1; i++)
	{	
		co.push_back(c[i]);
	}
	degree = d;
	coef = co;
}


POLYNOMIAL::POLYNOMIAL(int d, std::vector<COMPLEX> v)
{
	std::vector<COMPLEX> co;
	co.reserve(d+1);
	for (int i = 0; i<d+1; i++)
	{	
		co.push_back(v[i]);
	}
	degree = d;
	coef = co;
}

POLYNOMIAL::POLYNOMIAL(COMPLEX x)
{
	std::vector<COMPLEX> co;
	co.reserve(1);
	co.push_back(x);
	degree = 0;
	coef = co;
}



POLYNOMIAL::POLYNOMIAL()
{
}

POLYNOMIAL::~POLYNOMIAL()
{
}


COMPLEX POLYNOMIAL::operator () (const COMPLEX& z)
{
	COMPLEX fx = COMPLEX(0,0);
	for(int i=0; i < degree + 1; i++)
		fx = fx + power(z,i) * coef[i];
	return fx;
}

POLYNOMIAL operator + (const POLYNOMIAL& P , const POLYNOMIAL& Q)
{
	int d;
	if ( P.degree > Q.degree)
		d =  P.degree;
	else
		d = Q.degree;

	COMPLEX c[d+1];
	for(int i=0; i< d+1; i++)
	{
		if(i < P.degree +1 && i < Q.degree+1)
			c[i] = P.coef[i] + Q.coef[i];
		else if (i > P.degree)
			c[i] = Q.coef[i];
		else if (i > Q.degree)
			c[i] = P.coef[i];
	}

	POLYNOMIAL R = POLYNOMIAL(d, c);
	return R;
}

POLYNOMIAL operator - (const POLYNOMIAL& P)
{

	COMPLEX c[P.degree+1];
	for(int i=0; i< P.degree+1; i++)
		c[i] =  - P.coef[i];

	POLYNOMIAL R = POLYNOMIAL(P.degree, c);
	return R;
}
POLYNOMIAL operator - (const POLYNOMIAL& P, const POLYNOMIAL& Q)
{
	int d;
	int p = P.degree;
	int q = Q.degree;

	if (p > q)
		d = p;
	else
		d = q;

	COMPLEX c[d+1];

	for(int i=0; i< d+1; i++)
	{
		if(i < p +1 && i < q+1)
			c[i] = P.coef[i] - Q.coef[i];
		else if (i > p)
			c[i] = COMPLEX(0,0) - Q.coef[i];
		else if (i > q)
			c[i] = P.coef[i];
	}

	POLYNOMIAL R = POLYNOMIAL(d, c);
	return R;
}


POLYNOMIAL operator * (const POLYNOMIAL& P, const POLYNOMIAL& Q)
{
	int d = P.degree + Q.degree;
	int p = P.degree;
	int q = Q.degree;
	COMPLEX c[d+1];

	for (int i=0; i < p + 1; i++)
	{
		for(int j=0; j < q +1; j++)
		{
			c[i+j] = c[i+j] +  P.coef[i] * Q.coef[j];
		}
	}
	POLYNOMIAL R = POLYNOMIAL(d, c);
	
	return R;
}


COMPLEX evaluate(POLYNOMIAL P, COMPLEX z)
{
	COMPLEX fx = COMPLEX(0,0);
	for(int i=0;i<P.degree + 1;i++)
		fx = fx + power(z,i) * P.coef[i];
	return fx;

}

POLYNOMIAL deriv(POLYNOMIAL P, int k)
{
	int d;
	if (k == 0)
		return POLYNOMIAL(P.degree, P.coef);

	if (k > P.degree)
	{
		d = 0;
		COMPLEX C[1];
		C[0] = COMPLEX(0,0);
		return POLYNOMIAL(d, C);

	}
	else
	{
		d = P.degree - k;
		INTEGER tmp;
		COMPLEX C[d+1];
		for(int i=0; i< d+1; i++)
		{
			tmp = 1;
			for (int j=0; j< k; j++)
			{
				tmp = tmp * (j+i+1);
			}
			C[i] = P.coef[k+i] * REAL(tmp);
		}
	
		return POLYNOMIAL(d, C);
	}
}

// kth taylor coefficient of P at z
COMPLEX CoefAt(POLYNOMIAL P, int k, COMPLEX z)
{
	COMPLEX tmp[P.degree+1];
	for(int i=0; i<P.degree + 1; i++)
		tmp[i] = P.coef[i];

	for(int j=0; j<k; j++)
	{
		for(int i=0;i<P.degree-j;i++)
		{
			tmp[i] = tmp[i+1] * REAL(i+1);
		}
	}

	COMPLEX fx = COMPLEX(0,0);
	for(int i=0;i < P.degree + 1 - k; i++)
	{
		fx = fx + power(z,i) * tmp[i];
	}
	return  fx / REAL(factorial(k));
}	

// translation and dilation by g(x) = f(a + bx)
POLYNOMIAL translation(POLYNOMIAL P, REAL a, COMPLEX b)
{
	COMPLEX C[P.degree +1];
	for(int i=0; i<P.degree + 1; i++)
	{
		C[i] = power(a, i) * deriv(P,i)(b) / REAL(factorial(i));
	}

	return POLYNOMIAL(P.degree, C);
}



/*

POLYNOMIAL operator + (const POLYNOMIAL& P , const REAL& a)
{
	REAL c[P.degree+1];
	for(int i=0; i< P.degree+1; i++)
	{
		c[i] = P.coef[i];
	}
	c[0] += a;
	return POLYNOMIAL(P.degree, c);
}
POLYNOMIAL operator + (const REAL& a , const POLYNOMIAL& P)
{
	return P + a;
}


POLYNOMIAL operator - (const POLYNOMIAL& P, const REAL& c)
{
	return P + REAL(-c);
}
POLYNOMIAL operator - (const REAL& c, const POLYNOMIAL& P)
{
	return (-P) + c;
}

POLYNOMIAL operator * (const POLYNOMIAL& P , const REAL& a)
{
	REAL c[P.degree+1];
	for(int i=0; i< P.degree+1; i++)
	{
		c[i] = P.coef[i] * a;
	}
	return POLYNOMIAL(P.degree, c);
}
POLYNOMIAL operator * (const REAL& c, const POLYNOMIAL& P)
{
	return P*c;
}
POLYNOMIAL operator / (const POLYNOMIAL& P, const REAL& c)
{
	return P*REAL(1/c);
}

*/
