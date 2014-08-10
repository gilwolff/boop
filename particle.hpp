/*
 * particle class - stores index and position of particle
 * 
 */

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

class particle{
public:
    int index;
    boost::numeric::ublas::vector<double> r;
    particle(int i, double x, double y, double z) : r(3) {
	index=i;
	r(0)=x;
	r(1)=y;
	r(2)=z;
    }

    particle() : r(3) {}

    particle & operator= (const particle & other)
    {
        if (this != &other)
        {
	    index=other.index;
	    r(0)=other.r(0);
	    r(1)=other.r(1);
	    r(2)=other.r(2);
        }
        return *this;
    }
};

#endif
