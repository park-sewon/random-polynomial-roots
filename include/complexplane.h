/*
the header file defines subsets of complex planes:
- CLOSEDBOX : closed square box represented by center : COMPLEX  and width : RATIONAL
- R_CLOSEDBOX : closed square box with represented by center : R_COMPLEX and width : R_COMPLEX
- OPENDISC : open disc represented by center : COMPLEX and radius : COMPLEX
- R_OPENDISC : open disc represented by center : R_COMPLEX and radius : R_COMPLEX
- COMPONENT : A set of equal-sized adjacent R_CLOSEDBOX represented by a vector of R_CLOSEDBOX
*/

#ifndef COMPLEXPLANE_H
#define COMPLEXPLANE_H

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "rcomplex.h"

#include <utility>
using namespace iRRAM;


/*
CONSTRUCTORS: 
	- CLOSEDBOX(): constructs a closed square box with its width 1 centered at the origin 
	- CLOSEDBOX(c, w): constructs a closed square box with its width w centered at c
MEMBERS:
	- center: represents center of the square box
	- width: represents width of the square box
	- depth: depth of the square box in a subdivision tree
	- id: represents path of the square box in a subdivision tree
FUNCTIONS:
	- subdivide( i= 1,2,3,4 ) : subdivide the square box w.r.t. i; clock-wise 
*/
class CLOSEDBOX
{
	public:
		COMPLEX center;
		REAL width;


		int depth;
		std::string id;

		CLOSEDBOX();
		CLOSEDBOX(COMPLEX, REAL);
		~CLOSEDBOX();

		CLOSEDBOX subdivide(int ) const;

};

// similar to CLOSEDBOX
class R_CLOSEDBOX
{
	public:
		R_COMPLEX center;
		RATIONAL width;


		int depth;
		std::string id;

		R_CLOSEDBOX();
		R_CLOSEDBOX(R_COMPLEX, RATIONAL);
		~R_CLOSEDBOX();

		R_CLOSEDBOX subdivide(int )const;
};

/*
CONSTRUCTORS: 
	- OPENDISC(): constructs a disc with its width 1/2 centered at the origin 
	- OPENDISC(c, r): constructs a disc with its width r centered at c
	- OPENDISC(C : CLOSEDBOX): constructs a disc which its inner radius properly contains C

MEMBERS:
	- center: represents center of the disc
	- radius: represents width of the disc

FUNCTIONS:
	- multiply( r ) : return a disc whose radius is multiplied by r 
*/
class OPENDISC
{
	public:
		COMPLEX center;
		REAL radius;


		OPENDISC();
		OPENDISC(COMPLEX, REAL);
		OPENDISC(CLOSEDBOX);
		~OPENDISC();

		OPENDISC multiply(REAL) const;
};


// similar to OPENDISC
class R_OPENDISC
{
	public:
		R_COMPLEX center;
		RATIONAL radius;

		R_OPENDISC();
		R_OPENDISC(R_COMPLEX, RATIONAL);
		R_OPENDISC(R_CLOSEDBOX);
		~R_OPENDISC();

		R_OPENDISC multiply(RATIONAL) const;

};




/*
CONSTRUCTORS: 
	- COMPONENT(): constructs an empty COMPONENT
	- COMPONENT(B : R_CLOSEDBOX): constructs a COMPONENT with a single element B

MEMBERS:
	- box_list : list of (rational) closed boxes in the component C
	- upper_most: max Im(z) : z \in \union C[i]
	- right_most: max Re(z) : z \in \union C[i]
	- left_most: min Re(z) : z \in \union C[i]
	- below_most: min Im(z) : z \in \union C[i]
	- right_most_box argmax_B Re(z) | z \in B \and B \in C

	- Nc: Newton step factor (it gets squared in every successful newton step. hence it needs to be INTEGER not int)
	- kc: number of roots that the component contains
	- depth: depth of the component in the compoenent tree 

FUNCTIONS:
	- add(B): add a (rational) closed box into the component
	- size(): number of boxes in the component
	- Wc(): width of the smallest square box that contains the component
	- wc(): width of a box that the component contains
	- Mc(): center of the smallest square box that contains the component
	- Rc(): radius of the disc whose inner radius properly contains the smallest square box that contains the component
	- is_empty(): return whether the component if empty
	- split(): subdivide every boxes in the component
	- [i]: return ith box in the component
*/

class COMPONENT
{
	public:
		std::vector<R_CLOSEDBOX> box_list;

		RATIONAL upper_most;
		RATIONAL right_most;
		RATIONAL left_most;
		RATIONAL below_most;
		R_CLOSEDBOX right_most_box;

		INTEGER Nc;
		int kc;
		int depth;

		COMPONENT();
		COMPONENT(R_CLOSEDBOX );
		~COMPONENT();

		void add(R_CLOSEDBOX );
		
		int size() const;
		RATIONAL Wc() const;		
		RATIONAL wc() const;		
		R_COMPLEX Mc() const;		
		RATIONAL Rc() const;		
		bool is_empty() const;		
		COMPONENT split() const;
		R_CLOSEDBOX operator [](int ) const;

};

LAZY_BOOLEAN constains (OPENDISC, R_CLOSEDBOX);
LAZY_BOOLEAN is_in (COMPONENT, COMPLEX);
bool intersect(R_OPENDISC , R_OPENDISC );
bool comp_disc_intersect(COMPONENT, R_OPENDISC );
bool adj(COMPONENT , R_CLOSEDBOX );
void print(COMPONENT);
void printr(COMPONENT);





#endif
