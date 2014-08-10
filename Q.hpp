/*
 * wrapper for boost's spherical harmonics.
 *
 * Q_{lm}(r) is defined as Y_{lm}(theta,phi)
 * where theta and phi are the orientational
 * coordinates of a bond located at position r
 *
 * 
 */

#ifndef Q_HPP
#define Q_HPP

std::complex<double> Qlm(int l, int m, const std::vector<bond> & b){
    /*
     * averages over Qlm for all bonds in the vector-of-bonds b
     */
    std::complex<double> result(0,0);
    for(int i=0; i<b.size(); i++){
	result+=boost::math::spherical_harmonic(l,m, b[i].theta(), b[i].phi());
    }
    return result*(double)(1.0/b.size());
}

double Ql(int l, const std::vector<bond> & b){
    /*
     * Ql is defined as (4*pi/(2l+1)) * sum( abs(Qlm)^2 )
     * where the sum is from m=-l to m=l
     * and Qlm, I remind you, is already averaged over all bonds
     */
    double result=0;
    for(int m=-l; m<=l; m++){
	result+= 4*PI/(2l+1) * std::norm(Qlm(l,m,b));
    }
    return result;
}
#endif
