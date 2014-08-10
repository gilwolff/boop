/*
 * periodic - 1. 3d array with periodic boundary conditions
 *            2. periodic distance ("displacement")
 */

#ifndef PERIODIC_HPP
#define PERIODIC_HPP

//===============
// extend boost's vector to include periodic boundary conditions
//===============
template<typename T>
class periodic_vector_3d{
public:
    int d1,d2,d3;
    std::vector<T> data;
public:
    periodic_vector_3d(int d1=0, int  d2=0, int d3=0) :
        d1(d1), d2(d2), d3(d3), data(d1*d2*d3)
	{}
    
    T & operator()(int i, int j, int k) {
        return data[((i%d1+d1)%d1)*d2*d3 + ((j%d2+d2)%d2)*d3 + ((k%d3+d3)%d3)];
    }
    
    T const & operator()(int i, int j, int k) const {
        return data[((i%d1+d1)%d1)*d2*d3 + ((j%d2+d2)%d2)*d3 + ((k%d3+d3)%d3)];
    }

    void resize(int a, int b, int c){
	d1=a;
	d2=b;
	d3=c;
	data.resize(d1*d2*d3);
    }

    class iterator : public std::iterator<std::random_access_iterator_tag, T> {
    public:
	int i,j,k;
    public:
	iterator(unsigned int n){
	    i=n/(d2*d3);
	    j=(n%(d2*d3))/d3;
	    k=n%d3;
	}
	~iterator() {}

    // 	// The assignment and relational operators are straightforward
    // 	Iterator& operator=(const Iterator& other)
    // 	    {
    // 		node_ = other.node_;
    // 		return(*this);
    // 	    }

    // 	bool operator==(const Iterator& other)
    // 	    {
    // 		return(node_ == other.node_);
    // 	    }

    // 	bool operator!=(const Iterator& other)
    // 	    {
    // 		return(node_ != other.node_);
    // 	    }

    // 	// Update my state such that I refer to the next element in the
    // 	// SQueue.
    // 	Iterator& operator++()
    // 	    {
    // 		if (node_ != NULL)
    // 		{
    // 		    node_ = node_->next_;
    // 		}
    // 		return(*this);
    // 	    }

    // 	Iterator& operator++(int)
    // 	    {
    // 		Iterator tmp(*this);
    // 		++(*this);
    // 		return(tmp);
    // 	    }

    // 	// Return a reference to the value in the node.  I do this instead
    // 	// of returning by value so a caller can update the value in the
    // 	// node directly.
    // 	T& operator*()
    // 	    {
    // 		return(node_->getVal());
    // 	    }

    // 	T* operator->()
    // 	    {
    // 		return(&*(SQueue<T>::Iterator)*this);
    // 	    }
    };

    // Iterator begin()
    // 	{
    // 	    return(Iterator(root_));
    // 	}

    // Iterator end()
    // 	{
    // 	    return(Iterator(NULL));
//	}
};

boost::numeric::ublas::vector<double> periodic_displacement(double L_x, double L_y, double L_z,
							    boost::numeric::ublas::vector<double> r1, 
							    boost::numeric::ublas::vector<double> r2){
    /*
     * @param r1,r2 positions of two particles. assumed r1.size()=r2.size()=3
     * @param L_x,L_y,L_z size of box in the three dimensions.
     * returns dr=r2-r1 or periodic displacement of r2-r1
     */
    boost::numeric::ublas::vector<double> dr(3);
    dr=r2-r1;
    dr(0)=std::min(dr(0),L_x-dr(0));
    dr(1)=std::min(dr(1),L_y-dr(1));
    dr(2)=std::min(dr(2),L_z-dr(2));
    return dr;
}

#endif
