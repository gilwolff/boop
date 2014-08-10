
#include <iostream>

#include "boop.hpp"

int main(){
  
    boop::box my_box(0,62.42,0,62.42,0,62.42);//63
    my_box.init_cells(3.41);//3.401);
    std::cout<<"init_cells exited successfully\n";
    //	my_box.set_max_particles(650);
    my_box.set_boundary_conditions('p','p','p');
    my_box.import("xyz_pre.xyz");
    //my_box.validate_import();
    my_box.compute_bonds();
    for(int i=2; i<=12; i+=2){
	std::cout<<"the Q"<<i<<" order parameter is "<<my_box.compute_Ql(i)<<"\n";
    }
    // boop::periodic_vector_3d<int> v(2,3,4);
    // for( int n=0; n<24; n++){
    // 	v.data[n]=n;
    // }
    // for(int n=0; n<24; n++){
    // 	boop::periodic_vector_3d<int>::iterator it(n);
    // 	std::cout<<"i,j,k="<<it.i<<","<<it.j<<","<<it.k<<"\n";
    // }

    return 0;
}
