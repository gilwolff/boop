/*
 * bond class - stores position (midpoint) of bond
 * and orientation (theta, phi)
 * 
 */

#ifndef BOND_HPP
#define BOND_HPP

class bond{
    int bond_index;
    boost::numeric::ublas::vector<double> r;
    double bond_theta; //angle relative to z-axis, in range [0,pi]
    double bond_phi;   //angle in xy plane, in range [0,2*pi)
public:
    bond(int i,
	 const boost::numeric::ublas::vector<double> & r2,
	 const boost::numeric::ublas::vector<double> r1=boost::numeric::ublas::zero_vector<double>(3)){
	bond_index=i;
	r.resize(3);
	r=0.5*(r1+r2);
	boost::numeric::ublas::vector<double> dr;
	dr.resize(3);
	dr=r2-r1;
	bond_phi=atan2(dr(1) , dr(0)); //atan2 return in range [-pi,pi) so...
	if(bond_phi<0) bond_phi += 2*PI; //now phi is in range [0,2*pi)
	double rho=sqrt( dr(0)*dr(0) + dr(1)*dr(1) );
	bond_theta=atan2(rho, dr(2)); //since rho>=0, atan2 returns result in [0,pi] range
    }
    bond(){
    }

    int index() const{
	return bond_index;
    }
    boost::numeric::ublas::vector<double> position() const{
	return r;
    }
    double theta() const{
	return bond_theta;
    }
    double phi() const{
	return bond_phi;
    }

};

#endif
