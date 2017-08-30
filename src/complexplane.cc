#include "../include/complexplane.h"
#include "../include/utilities.h"
#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <stdexcept>

#include <stack>
#include <queue>
#include <utility>
#include <vector>

using namespace iRRAM;

/*
 *
 * Member functions for Closed box in a complex plane
 *
 */
CLOSEDBOX::CLOSEDBOX()
{
	center = COMPLEX(0,0);
	width = 1;
	depth = 1;
}

CLOSEDBOX::CLOSEDBOX(COMPLEX c, REAL r)
{
	center = c;
	width = r;
	depth = 1;
}
CLOSEDBOX::~CLOSEDBOX()
{
}


CLOSEDBOX  CLOSEDBOX::subdivide(int n) const
{
	CLOSEDBOX S;
	if (n == 1)
		S =  CLOSEDBOX(center + COMPLEX(width/4,width/4), width / 2);
	else if (n == 2)
		S =  CLOSEDBOX(center + COMPLEX(width/4,0-width/4), width / 2);
	else if (n == 3)
		S =  CLOSEDBOX(center - COMPLEX(width/4,width/4), width / 2);
	else if (n == 4)
		S =  CLOSEDBOX(center + COMPLEX(0-width/4,width/4), width / 2);

	S.depth = depth + 1;
	return S;
}



/*
 *
 * Member functions for open disc
 *
 */
OPENDISC::OPENDISC()
{
	center = COMPLEX(0,0);
	radius = 1/2;
}

OPENDISC::OPENDISC(COMPLEX c, REAL r)
{
	center = c;
	radius = r;
}

OPENDISC::OPENDISC(CLOSEDBOX B)
{
	center = B.center;
	radius = B.width * 3 / 4;
}

OPENDISC::~OPENDISC()
{
}

OPENDISC  OPENDISC::multiply(REAL f) const
{
	return OPENDISC(center, radius*f); 
}


/*
 *
 * Member functions for closed box with rational center and width
 *
 */
R_CLOSEDBOX::R_CLOSEDBOX()
{
	center = R_COMPLEX(0,0);
	width = 1;
	depth = 1;
	id = "0";
}

R_CLOSEDBOX::R_CLOSEDBOX(R_COMPLEX c, RATIONAL r)
{
	center = c;
	width = r;
	depth = 1;
	id = "0";
}

R_CLOSEDBOX::~R_CLOSEDBOX()
{
}



R_CLOSEDBOX  R_CLOSEDBOX::subdivide(int n) const
{
	R_CLOSEDBOX S;
	if (n == 1)
		S =  R_CLOSEDBOX(center + R_COMPLEX(width/4,width/4), width / 2);
	else if (n == 2)
		S =  R_CLOSEDBOX(center + R_COMPLEX(width/4,0-width/4), width / 2);
	else if (n == 3)
		S =  R_CLOSEDBOX(center - R_COMPLEX(width/4,width/4), width / 2);
	else if (n == 4)
		S =  R_CLOSEDBOX(center + R_COMPLEX(0-width/4,width/4), width / 2);

	S.depth = depth + 1;
	S.id = id;
	S.id.append(std::to_string(n));
	return S;
}


/*
 *
 * Member functions for open disc with rational center and rational radius
 *
 */
R_OPENDISC::R_OPENDISC()
{
	center = R_COMPLEX();
	radius = RATIONAL(1, 2);
}

R_OPENDISC::R_OPENDISC(R_COMPLEX c, RATIONAL r)
{
	center = c;
	radius = r;
}

R_OPENDISC::R_OPENDISC(R_CLOSEDBOX B)
{
	center = B.center;
	radius = B.width * 3 / 4;
}

R_OPENDISC::~R_OPENDISC()
{
}

R_OPENDISC  R_OPENDISC::multiply(RATIONAL f) const
{
	return R_OPENDISC(center, radius*f); 
}

// COMPONENT MEMBER FUNCTIONS
COMPONENT::COMPONENT()
{
	depth = 0;
}

COMPONENT::COMPONENT(R_CLOSEDBOX B)
{
	upper_most = B.center.imag() + B.width / 2;
	right_most = B.center.real() + B.width / 2;	
	below_most = B.center.imag() - B.width / 2;
	left_most  = B.center.real() - B.width / 2;
	right_most_box = B;
	box_list.push_back(B);
}

COMPONENT::~COMPONENT(){}

void COMPONENT::add(R_CLOSEDBOX B)
{
	RATIONAL upper_most_ = B.center.imag() + B.width / 2;
	RATIONAL right_most_ = B.center.real() + B.width / 2;	
	RATIONAL below_most_ = B.center.imag() - B.width / 2;
	RATIONAL left_most_  = B.center.real() - B.width / 2;

	if (upper_most_ > upper_most)
		upper_most = upper_most_;
	if (right_most_ > right_most)
	{
		right_most = right_most_;
		right_most_box = B;
	}
	if (below_most_ < below_most)
		below_most = below_most_;
	if (left_most_ < left_most)
		left_most = left_most_;

	box_list.push_back(B);
}

int COMPONENT::size() const
{
	return box_list.size();
}

RATIONAL COMPONENT::Wc() const
{
	return maximum(right_most - left_most, upper_most - below_most);
}

RATIONAL COMPONENT::wc() const
{
	return box_list.back().width;
}

R_COMPLEX COMPONENT::Mc() const
{
	return R_COMPLEX((right_most+left_most)/2,(upper_most+below_most)/2);
}

RATIONAL COMPONENT::Rc() const
{
	return Wc() * 3 / 4;
}

bool COMPONENT::is_empty() const
{
	return box_list.empty();
}

COMPONENT COMPONENT::split() const
{
	COMPONENT N;

	N.Nc = Nc;
	N.kc = kc;
	N.depth = depth;

	for(int i=0; i < (int)box_list.size(); i++)
	{
		for(int p=1; p<5; p++)
			N.add(box_list[i].subdivide(p));
	}

	N.depth += 1;
	return N;

}

R_CLOSEDBOX COMPONENT::operator [](int i) const
{
	return box_list[i];
}

LAZY_BOOLEAN constains (OPENDISC D, R_CLOSEDBOX B)
{
	return imag(D.center) + D.radius < B.center.imag() + B.width/2 &&
			imag(D.center) - D.radius > B.center.imag() - B.width/2 &&
			real(D.center) + D.radius < B.center.real() + B.width/2 &&
			real(D.center) - D.radius > B.center.real() - B.width/2;

}

LAZY_BOOLEAN is_in (COMPONENT C, COMPLEX x)
{
	for (int i=0; i < C.size(); i++)
	{
		REAL dist = maximum(abs(C[i].center.real() -real(x)), abs(C[i].center.imag() - imag(x)));
		if( choose( dist < C[i].width, dist > C[i].width / 2 )==1)
			return true;
	}
	return false;


}



bool intersect(R_OPENDISC A, R_OPENDISC B)
{
	return (A.center.real() - B.center.real()) * (A.center.real() - B.center.real()) + 
	 (A.center.imag() - B.center.imag()) * (A.center.imag() - B.center.imag())
	 < (A.radius + B.radius) * (A.radius + B.radius);
}


bool intersect(R_CLOSEDBOX A, R_CLOSEDBOX B)
{
	return maximum(abs(A.center.real() - B.center.real()), abs(A.center.imag() - B.center.imag())) <= (A.width+B.width) / 2;
}


bool comp_disc_intersect(COMPONENT C, R_OPENDISC D)
{
	for (int i=0; i<C.size(); i++)
		if (intersect(C[i], D))
			return true;
	return false;
}

bool adj(R_CLOSEDBOX A, R_CLOSEDBOX B)
{
	return maximum(abs(A.center.real() - B.center.real()), abs(A.center.imag() - B.center.imag())) == (A.width+B.width) / 2;
}

bool adj(COMPONENT C, R_CLOSEDBOX I)
{
	for (int i=0; i<C.size(); i++)
		if (adj(C[i], I))
			return true;
	return false;
}

void print(COMPONENT C)
{
	for (int i=0; i<(int) C.box_list.size(); i++)
		cout<<REAL(C.box_list[i].center.real()) <<", "<<REAL(C.box_list[i].center.imag())<<", "<<REAL(C.box_list[i].width) <<"\n";
}
void printr(COMPONENT C)
{
	for (int i=0; i<(int) C.box_list.size(); i++)
		cout<<C.box_list[i].center.real() <<", "<<C.box_list[i].center.imag()<<", "<<C.box_list[i].width <<"\n";
}

