#include "../include/polynomialroot.h"
#include "../include/complexplane.h"
#include "../include/utilities.h"

#include "iRRAM/lib.h"
#include "iRRAM/core.h"

#include <stack>
#include <queue>
#include <utility>
#include <vector>
#include <cmath>
#include <time.h>

using namespace iRRAM;


// < relation to use pre-impelented sorting function and prioirty queue.
/*
struct less_than_key
{
    inline bool operator() (const std::pair<R_OPENDISC, int>& C1, 
    	const std::pair<R_OPENDISC, int>& C2)
	{
		return (C1.first.center < C2.first.center);
	}
};

struct IC_LESS
{
    inline bool operator() (const INTERVALCOMPONENT& C1, 
    	const INTERVALCOMPONENT& C2)
	{
		return (C1.lower + C1.upper < C2.lower + C2.upper);
	}
};
*/

struct disc_wider
{
    inline bool operator() (const std::pair<R_OPENDISC, int>& C1, 
    	const std::pair<R_OPENDISC, int>& C2)
	{
		return (C1.first.radius < C2.first.radius);
	}
};

struct box_wider
{
    inline bool operator() (const std::pair<R_CLOSEDBOX, int>& C1, 
    	const std::pair<R_CLOSEDBOX, int>& C2)
	{
		return (C1.first.width < C2.first.width);
	}
};

struct component_wider
{
    inline bool operator() (const COMPONENT& C1, 
    	const COMPONENT& C2)
	{
		return (C1.Wc() < C2.Wc());
	}
};


/*

std::vector<std::pair<RATIONALINTERVAL, int> > 
casting(std::vector<INTERVALCOMPONENT> P)
{

	std::vector<std::pair<RATIONALINTERVAL, int> > Q;

	for (int i =0; i<(int) P.size(); i++)
	{	
		Q.push_back(std::pair<RATIONALINTERVAL, int> (RATIONALINTERVAL(P[i].Mc(), P[i].Wc()/2), P[i].kc));
	}
	return Q;
} 
*/


/******************
GRAEFFE'S ITERATION
*******************/

POLYNOMIAL polyE(POLYNOMIAL P)
{
	double td = P.degree / 2;
	int d = std::floor(td);

	COMPLEX C[d+1];

	for(int i = 0; i<d+1; i++)
	{
		C[i] = P.coef[i*2];
	}
	return POLYNOMIAL(td, C);
}

POLYNOMIAL polyO(POLYNOMIAL P)
{
	double td = (P.degree-1) / 2;
	int d = std::floor(td);

	COMPLEX C[d+1];

	for(int i = 0; i<d+1; i++)
	{
		C[i] = P.coef[1+i*2];
	}
	return POLYNOMIAL(td, C);
}

POLYNOMIAL Giterate(POLYNOMIAL P, int N)
{

	POLYNOMIAL Q = POLYNOMIAL(P.degree, P.coef);
	COMPLEX c[2];
	c[0] = COMPLEX(0,0); c[1] = COMPLEX(1,0);
	int k = pow(-1, P.degree);
	POLYNOMIAL T = POLYNOMIAL(1, c);
	if(k > 0)
		for(int i=1; i<N+1; i++)
		{
			Q = (polyE(Q) * polyE(Q)) - (T * (polyO(Q) * polyO(Q)));
		}
	else
		for(int i=0; i<N; i++)
		{
			Q = (T*(polyO(Q) * polyO(Q))) - (polyE(Q) * polyE(Q));
		}
	return Q;
}




/*
COMPUTABLE PREDICATES FOR ROOT CONTAINMENT OF AN OPEN DISC 
*/

bool unit_disc_test(POLYNOMIAL P, int k)
{
	REAL LHS;
	REAL RHS = 0;
	for (int i=0; i<P.degree+1; i++)
		if(i!=k)
			RHS += abs(P.coef[i]);
	LHS = abs(P.coef[k]);
	int ans = choose(LHS>RHS, LHS<RHS, LHS*2 < RHS*3 && 2*RHS < 3*LHS) == 1;
	return ans;
}

// Soft pallet test: \tilde{\mathcal{T}}_k(P,D)
bool softTTest(POLYNOMIAL P, int k, R_OPENDISC D)
{
	COMPLEX m = COMPLEX(D.center.real(), D.center.imag());
	REAL r = D.radius;
	REAL LHS = abs(CoefAt(P,k,m)) * power(r,k);
	REAL RHS = 0;
	bool ans;
	for(int i=0;i<P.degree+1;i++)
	{
		if(i!=k)
		{
			RHS += abs(CoefAt(P,i,m)) * power(r,i);
		}
	}
	ans = choose(LHS>RHS, LHS<RHS, LHS*2 < RHS*3 && 2*RHS < 3*LHS) == 1;
	return ans;
}

// Soft pallet test on Graeffe iteration: \tilde{\mathcal{T}}_k(G_{\log(1+\log d) + 5}(P_D),D(0,1))
// on a interval with rational endpoints
bool softGTest(POLYNOMIAL P, int k, R_OPENDISC D)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, D.radius, COMPLEX(D.center.real(), D.center.imag())), N);
	// R_OPENDISC U(R_COMPLEX(0, 0), 1);
	// return softTTest(Q, k, U);
	return unit_disc_test(Q,k);
}


// Soft pallet test on Graeffe iteration: \tilde{\mathcal{T}}_k(G_{\log(1+\log d) + 5}(P_D),D(0,1))
// on a interval with real endpoints
bool softGTest(POLYNOMIAL P, int k, OPENDISC D)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, D.radius, D.center), N);
	// R_OPENDISC U(R_COMPLEX(0, 0), 1);
	// return softTTest(Q, k, U);
	return unit_disc_test(Q,k);

}


// First i such that T_i(G_N(f_{m,r}), 0, 1) /\ T_i(G_N(f_{m,3r}), 0, 1)holds 
//    0 if T_0(G_N(f_{m,r}), 0, 1) 
//   -1 if - T_i(G_N(f_{m,r}), 0, 1) for all i
// -i-1 if T_i(G_N(f_{m,r}), 0, 1) /\  -T_i(G_N(f_{m,3r}), 0, 1)
int softGStarThree(POLYNOMIAL P, R_OPENDISC B, int k)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, B.radius, COMPLEX(B.center.real(), B.center.imag())), N);
	POLYNOMIAL R = Giterate(translation(P, 3*B.radius, COMPLEX(B.center.real(), B.center.imag())), N);

	// R_OPENDISC U(R_COMPLEX(0, 0), 1);

	// for(int i =0;i<k+1;i++)
	// {
	// 	if(softTTest(Q, i, U))
	// 	{
	// 		if (i == 0)
	// 			return 0;
	// 		if(softTTest(R, i, U))
	// 			return i;
	// 		else
	// 			return -i-1;
	// 	}
	// }

	for(int i =0;i<k+1;i++)
	{
		if(unit_disc_test(Q, i))
		{
			if (i == 0)
				return 0;
			if(unit_disc_test(R, i))
				return i;
			else
				return -i-1;
		}
	}



	return - 1;
}

// First i such that \tilde{T}_i(G_N(f_{m,r}), 0, 1) holds 
//  -1 if all fail
int softGStar(POLYNOMIAL P, R_OPENDISC B, int k)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, B.radius, COMPLEX(B.center.real(), B.center.imag())), N);
	R_OPENDISC U(R_COMPLEX(0, 0), 1);

	for(int i =0;i<k+1;i++)
		// if(softTTest(Q, i, U))
		if(unit_disc_test(Q, i))

			return i;
	return - 1;
}






// Return a list of (D_i, d_i) such that each D_i is disjoint (1,3)--isolating for positive number of roots.
// Implementation based on simple naive subdivision algorithm with Graeffe's iteration
std::vector<R_COMPLEX >
root_approximation_naive(int prec, POLYNOMIAL P)
{
	int n = P.degree;
	RATIONAL min = 1;


/*
COMPUTING ROOT BOUND D
*/

	R_OPENDISC D = R_OPENDISC();
	while(!softGTest(P, n, D))
	{
		D = D.multiply(2);
	}
	// cout <<D.radius<<"\n";

	REAL epsilon = power(2, prec);
    std::priority_queue<std::pair<R_CLOSEDBOX, int> , std::vector<std::pair<R_CLOSEDBOX, int> >, box_wider> Q0;
    Q0.push(std::pair<R_CLOSEDBOX, int>(R_CLOSEDBOX(R_COMPLEX(0,0), 2*D.radius), P.degree));

	int k;


	int nroot=0;
	int f = 0;
	int ff = 0;

	std::vector<std::pair<R_OPENDISC, int> > solutionDiscs;
	solutionDiscs.reserve(P.degree);

	cout << D.radius <<"\n";

/*
MAIN LOOP 
*/
	while(!Q0.empty())
	{
		std::pair<R_CLOSEDBOX, int> fp = Q0.top();	
		Q0.pop();
		
		R_CLOSEDBOX a = fp.first;
		R_OPENDISC b = R_OPENDISC(a);


		int max_multiplicity = fp.second;

		ff = 0;
		for (int i = 0; i < solutionDiscs.size(); i++)
		{
			if (intersect(b, solutionDiscs[i].first))
			{
				ff = 1; break;
			}
		}

		cout<< a.id<<"\n";
		if (ff == 0) 
		{
			k = softGStarThree(P, b, max_multiplicity);

			if (k==-1)   // failed
			{	
				for(int si=1;si<5;si++)
					Q0.push(std::pair<R_CLOSEDBOX, int>(a.subdivide(si), max_multiplicity));

			}
			else if (k < 0)	// T(3D) failed but T(D) succeeded with i = -k - 1
			{
				for(int si=1;si<5;si++)
					Q0.push(std::pair<R_CLOSEDBOX, int>(a.subdivide(si), -k - 1));

			}
			else if (k != 0) // D, 3D has k > 0 roots
			{
				if(choose(b.radius < epsilon, b.radius > epsilon / 2) == 2)
				{
					for(int si=1;si<5;si++)
						Q0.push(std::pair<R_CLOSEDBOX, int>(a.subdivide(si), k));
				}
				else
				{
					solutionDiscs.push_back(std::pair<R_OPENDISC, int> (b, k));
					nroot+=k;
				}
			}
		}
	}

	std::vector<R_COMPLEX > ansvec;
	for (int i=0; i<(int) solutionDiscs.size(); i++)
	{
		for (int j=0; j<solutionDiscs[i].second; j++)
			ansvec.push_back(solutionDiscs[i].first.center);
	}
	return ansvec;
}


/*
CYAP2016 NEWTON METHOD
*/



COMPONENT
Newton(POLYNOMIAL P, COMPONENT C, REAL epsilon)
{

	REAL L, R;
	COMPLEX Xprime; // pick a point which is C.wc/2 away from the boundary of the commoment
	COMPONENT D;
	R_CLOSEDBOX ans;

	RATIONAL new_left, new_below, new_width;

	R_COMPLEX R_Xc = C.right_most_box.center + R_COMPLEX(C.wc()/2, 0);
	COMPLEX Xc = COMPLEX(R_Xc.real(), R_Xc.imag());
	RATIONAL r = C.Wc() / 2;

	L = 4 * r * abs(deriv(P, 1)(Xc));
	R = abs(P(Xc));

	if( choose(L>R, L<R, 2*L< 3*R && 2*R < 3*L) != 2) // when f'(Xc) != 0
	{

		Xprime = Xc - REAL(C.kc) * P(Xc)/ deriv(P,1)(Xc);
		if (is_in(C, Xprime))
		{

			OPENDISC II = OPENDISC(Xprime, REAL( C.wc()) / 8 / C.Nc);

			if( softGTest(P, C.kc, II))
			{

				new_width = C.wc() / (2*C.Nc);

				INTEGER index_x = ( (real(Xprime) - C.left_most ) / new_width ).as_INTEGER();
				INTEGER index_y = ( (imag(Xprime) - C.below_most) / new_width ).as_INTEGER();

				if (index_x == 0)
					index_x = index_x + 1;
				if (index_x == (C.right_most - C.left_most / C.wc()) )
					index_x = index_x - 1;

				if (index_y == 0)
					index_y = index_y + 1;
				if (index_y == (C.upper_most - C.below_most / C.wc()) )
					index_y = index_y - 1;

				new_left  = C.left_most + new_width * (index_x - 1);
				new_below = C.below_most + new_width * (index_y - 1);
				ans = R_CLOSEDBOX(R_COMPLEX(new_left+new_width, new_below+new_width),2*new_width);

				D = COMPONENT();
				D.add(ans);
				D = D.split();
				D.Nc = C.Nc * C.Nc;
				D.depth = C.depth+1;
			}

		}

	}
	return D;
}


std::vector< COMPONENT > 
Bisect(POLYNOMIAL P, COMPONENT Comp)
{
	COMPONENT C;
	std::vector< COMPONENT > U;
	R_CLOSEDBOX B;

	int flg;
	bool specialflg = true;

	C = Comp.split();  //C.depth ++

	for (int i =0; i< C.size(); i++)
	{
		B = C[i];
		if (softGTest(P,0, R_OPENDISC(B)) == false)
		{
			flg = 0;
			for(int j=0; j < (int)U.size(); j++)
			{	
				if(adj(U[j], B))
				{	
					U[j].add(B);
					flg = 1;
					break;
				}
			}
			if (flg == 0)	// If there does not exists a adjoint component in U
			{
				COMPONENT T = COMPONENT(B);
				T.depth = C.depth;
				U.push_back(T);
			}

		}
	}
	if ((int)U.size() ==1)
		specialflg = false;

	for (int i=0; i<(int)U.size(); i++)
	{
		if (specialflg)
			U[i].Nc = 4;
		else
			if (sqrt(C.Nc) > 4)
				U[i].Nc = sqrt(C.Nc);
			else
				U[i].Nc = 4;
	}
	return U;
}



std::vector<R_COMPLEX >
root_approximation_newton(int prec, POLYNOMIAL P)
{

	int n = P.degree;

	/*
		COMPUTING ROOT BOUND
	*/
	R_OPENDISC D = R_OPENDISC(R_COMPLEX(0,0),1);
	while(!softGTest(P, P.degree, D))
	{
		D = D.multiply(2);
	}

	// std::cout <<"dd\n";

	std::vector< COMPONENT > Q_out, Q_main;
	std::vector< COMPONENT > bisect_result;
    // std::priority_queue<COMPONENT , std::vector<COMPONENT >, component_wider> Q_main;


	COMPONENT C = COMPONENT(R_CLOSEDBOX(D.center, D.radius*2));
	C.Nc = 4;
	Q_main.push_back(C);
	COMPONENT Cprime;

	REAL epsilon = power(2,prec);


	int pi;

	RATIONAL tmp, qtmp;
	COMPONENT tc;

	cout << setRwidth(30);

	int fll;
	while(!Q_main.empty())
	{

// 		Swap Q1[end -1] <-> Q1[max wc]
		pi = 0;
		qtmp = 0;
		tmp = 0;
		for(int i=0; i< (int) Q_main.size(); i++)
		{
			qtmp = Q_main[i].Wc();
			if (qtmp > tmp)
			{
				pi = i;
				tmp = qtmp;
			}
		}
		tc = Q_main[pi];
		Q_main[pi] = Q_main[(int) Q_main.size() - 1];
		Q_main[(int) Q_main.size() -1] = tc;
// 		swap done


		C = Q_main.back();
		Q_main.pop_back();
		print(C);
		
		fll = 0;
		R_OPENDISC II = R_OPENDISC(C.Mc(),C.Rc());

		for(int i=0; i< (int) Q_main.size(); i++)
		{
			if(intersect(R_OPENDISC(Q_main[i].Mc(), Q_main[i].Rc()), II.multiply(4)))
			{
				fll = 1;
				break;
			}
		}

		if (fll == 0)
		{
			C.kc = softGStar(P, II, P.degree);
			
			if(C.kc >0)
			{
				if(choose(C.Wc() > epsilon, C.Wc() < epsilon*2) == 1)
				{
					Cprime= Newton(P,C,epsilon);

					if(Cprime.is_empty() == false)
					{
						cout <<" newton worked\n";
						Q_main.push_back(Cprime);
						continue;
					}
				}
				else if (C.Wc() <= 3 * C.wc())
				{
					Q_out.push_back(C);
					continue;
				}
			}
		}

		bisect_result = Bisect(P,C);
		for(int i=0; i<(int)bisect_result.size(); i++)
			Q_main.push_back(bisect_result[i]);


	}


	std::vector< R_COMPLEX > roots;
	for (int i =0; i<(int) Q_out.size(); i++)
		for (int j=0; j< Q_out[i].kc; j++)
			roots.push_back(Q_out[i].Mc());



	return roots;
}

