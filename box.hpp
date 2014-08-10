/*
 * box - a container for the particles 
 * 
 */
#ifndef BOX_HPP
#define BOX_HPP

class box{
    double min_x,min_y,min_z,max_x,max_y,max_z;
    double L_x, L_y, L_z;
    int num_cells_x,num_cells_y,num_cells_z,num_cells;
    int num_particles, max_particles;
    double bond_cutoff;
    char boundary_condition_x, boundary_condition_y,boundary_condition_z;
    periodic_vector_3d<std::vector<particle> > cell;
    std::vector<bond> bonds;
public:
    // ==============================
    // constructor for box
    // ==============================
    box(double x0, double x1, double y0, double y1, double z0, double z1){
	min_x=x0;
	max_x=x1;
	min_y=y0;
	max_y=y1;
	min_z=z0;
	max_z=z1;
	L_x=max_x-min_x;
	L_y=max_y-min_y;
	L_z=max_z-min_z;
	num_particles=0;
	max_particles=-1; //i.e. default is no limit on num. of particles
	boundary_condition_x=boundary_condition_y=boundary_condition_z='o'; //default is open boundary conditions in all directions
    }

    // ==============================
    // set_max_particles
    // ==============================
    void set_max_particles(int max_N){
	/*
	 * @param max_N maximum number of particles
	 * effect: sets private variable max_particles to max_N
	 */
	max_particles=max_N;
    }
    // ==============================
    // set_boundary_conditions
    // ==============================
    void set_boundary_conditions(char c_x, char c_y, char c_z){
	/*
	 * @param c_x,c_y,c_z boundary conditions in x,y,z directions. 'o' for open, 'p' for periodic
	 * effect: 
	 *   changes boundary_condition_x,boundary_condition_y,boundary_condition_z.
	 *   if one of c_x,c_y,c_z is not 'o' or 'p' the program exits
	 */
	if((c_x !='o' && c_x!='p') || (c_y!='o' && c_y !='p') || (c_z!='o' && c_z!='p')){
	    std::cout<<"illegal boundary condition. should be \'o\' or \'p\'\n";
	    exit(1);
	}
	boundary_condition_x=c_x;
	boundary_condition_y=c_y;
	boundary_condition_z=c_z;
    }

    // ==============================
    // import
    // ==============================
    void import(const char* file_name){
	/*
	 * @param file_name name of file, assumed in xyz format, from which to read particle locations.
	 *   first two lines are ignored and the first letter of all consecutive lines ignored
	 * effect: 
	 *   stores particles in vector<particle> cell(i,j,k) 
	 *   number of particles stored in num_particles
	 *   number of particles in each cell is equivlaent to cell(i,j,k).size() 
	 */
	std::fstream file_xyz;
	file_xyz.open(file_name);
	char c;
	double x,y,z;
	file_xyz.ignore(256,'\n');
	file_xyz.ignore(256,'\n');
	while(file_xyz>>c>>x>>y>>z){
	    boost::numeric::ublas::vector<int> i(3);
	    particle p(num_particles,x,y,z);
	    i=identify_cell(p);
	    cell(i(0),i(1),i(2)).push_back(p);
	    num_particles++;
	    if(max_particles==num_particles)
		break;
	}
	file_xyz.close();
	
	std::cout<<"import counted "<<num_particles<<" and exited successfully\n";
    }
	
    // ==============================
    // init_cells
    // ==============================
    void init_cells(double bond_cut){
	/*
	 * input:
	 *   @param bond_cut the bond cutoff length
	 * effect:
	 *   initializes 3D array 'cell'to dimensions num_cells_x * num_cells_y * num_cells_z
	 *   changes num_cells_x,num_cells_y,num_cells_z
	 * postcondition:  
	 *   the length of each cell, in x,y,z directions must be no less than bond_cut
	 */
	bond_cutoff=bond_cut;
	if (bond_cutoff+1<=1) { std::cout<<"error: bond cutoff non-positive. exiting\n"; exit(1); }
	std::cout<<"starting init_cells\n";
	//the bond_cutoff serves as the (minimal) length of the cell
	num_cells_x = (int) (L_x/bond_cutoff);
	num_cells_y = (int) (L_y/bond_cutoff);
	num_cells_z = (int) (L_z/bond_cutoff);

	cell.resize(num_cells_x, num_cells_y, num_cells_z);
    }
	
    // ==============================
    // identify_cell
    // ==============================
    boost::numeric::ublas::vector<int> identify_cell(const particle & p){
	/*
	 * input:
	 *   @param p particle whose location we query to find out which cell it belongs to
	 * output:
	 *   pointer to cell(i,j,k)
	 */
	const boost::numeric::ublas::vector<double>& r=p.r;
	boost::numeric::ublas::vector<int> i(3);
	i(0)=(int) (num_cells_x*(r(0)-min_x)/L_x);
	i(1)=(int) (num_cells_y*(r(1)-min_y)/L_y);
	i(2)=(int) (num_cells_z*(r(2)-min_z)/L_z);
	if(i(0)>=num_cells_x || i(1)>= num_cells_y || i(2)>=num_cells_z){
	    std::cout<<"error: out of bounds in identify cell\n";
	    exit(1);
	}
		
	return i;
    }
	
	
    // ==============================
    // compute_bonds
    // ==============================
    void compute_bonds(){
	/*
	 * effect:
	 *   bonds[i] stores the i-th bond, and bonds.size() is the number of bonds.
	 *   if boundary_condition_x is set to 'o', then bulk-bulk and bulk-boundary bonds are computed in the x direction but not boundary-boundary
	 *   if boundary_condition_x is set to 'p', then bulk-bulk, bulk-boundary and boundary-boundary bonds are computed in the x direction
	 *   here a boundary cell in the x direction is a cell(i,j,k) with i==0 or i==num_cells_x-1
	 *
	 * method:
	 *   The plan is to cycle through cells. For each cell cycle through particles in it.
	 *   For each of those particles, cycle through all particles in all neighboring cells.
	 *   Mark a particle that checked all adjacent particles as "used" so bonds won't be counted twice.
	 *   For open boundary conditions cycle through bulk cells (neighbor cells may belong to boundary).
	 *   For periodic boundary conditions cycle through all cells
	 */
		
		
	std::vector<int> used(num_particles); //sets all values to zero
	int i, i_i, i_f, ii; 
	int j, j_i, j_f, jj;
	int k, k_i, k_f, kk;
	int count=0;
	particle p1,p2;
	boundary_condition_x=='o' ? (i_i=1,i_f=num_cells_x-2) : (i_i=0,i_f=num_cells_x-1);
	boundary_condition_y=='o' ? (j_i=1,j_f=num_cells_y-2) : (j_i=0,j_f=num_cells_y-1);
	boundary_condition_z=='o' ? (k_i=1,k_f=num_cells_z-2) : (k_i=0,k_f=num_cells_z-1);

	for(i=i_i; i<=i_f; i++) for(j=j_i; j<=j_f; j++) for(k=k_i; k<=k_f; k++){
		    //boost::numeric::ublas::vector<particle>::iterator it(cell(i,j,k).begin()), itend(cell(i,j,k).end());
		    for(int num1=0; num1<cell(i,j,k).size(); num1++){
			p1=cell(i,j,k)[num1];
			used[p1.index]=1;
			for (ii=i-1; ii<=i+1; ii++) for (jj=j-1; jj<=j+1; jj++) for (kk=k-1; kk<=k+1; kk++){
				    // boost::numeric::ublas::vector<particle>::iterator nit(cell(i+n[0],j+n[1],k+n[2]).begin()), nitend(cell(i+n[0],j+n[1],k+n[2]).end());
				    for(int num2=0; num2<cell(ii,jj,kk).size(); num2++){
					p2=cell(ii,jj,kk)[num2];
					boost::numeric::ublas::vector<double> dr(3);
					dr=periodic_displacement(L_x,L_y,L_z,p1.r,p2.r);
					if(!used[p2.index] && boost::numeric::ublas::norm_2(dr)<=bond_cutoff){
					    bonds.push_back(bond(count,dr));
					    ++count;
					}
				    }//end num2
				}//end ii,jj,kk
		    }//end num1
		}
	std::cout<<"counted "<<count<<" bonds\n";
	std::cout<<"compute_bonds exited successfully\n";
    }
    // ==============================
    // Compute Steinhardt's Ql order parameter
    // ==============================
    double compute_Ql(int l){
	/*
	 * input:
	 *   @param l the index of stienhardt's function $Q_{l}$
	 * output:
	 *   value of stienhardt's function $Q_{l}$
	 */
	return Ql(l, bonds);
    }


    // ==============================
    // validation stuff
    // ==============================
    void validate_init_cells(){
	/*
	 * 
	 */
	for(int i=0; i<num_cells_x; i++) for(int j=0; j<num_cells_y; j++) for(int k=0; k<num_cells_z; k++){
		    std::cout<<"hello from cell "<<i<<","<<j<<","<<k<<"\n";
		    std::cout<<"my boundaries are: x=["<<i*L_x/num_cells_x<<","<<(i+1)*L_x/num_cells_x<<") ";
		    std::cout<<"y=["<<j*L_y/num_cells_y<<","<<(j+1)*L_y/num_cells_y<<") ";
		    std::cout<<"z=["<<k*L_z/num_cells_z<<","<<(k+1)*L_z/num_cells_z<<")\n";
					
		    std::vector<particle>::iterator it(cell(i,j,k).begin()), itend(cell(i,j,k).end());
		    for(; it!=itend; it++){
			std::cout<<(*it).r<<" ";
		    }
		    std::cout<<"\n";
		}
    }

    void validate_bonds(){
	//boost::numeric::ublas::vector<bond>::iterator bit(bonds.begin()), bitend(bonds.end());
    }
    
    void validate_import(){
	for(int i=0; i<num_cells_x; i++) for(int j=0; j<num_cells_y; j++) for(int k=0; k<num_cells_z; k++){
		    std::vector<particle>::iterator it(cell(i,j,k).begin()), itend(cell(i,j,k).end());
		    std::cout<<"~~~~~~~~in cell "<<i<<","<<j<<","<<k<<"~~~~~~~~\n";
		    for(; it!=itend; it++){
			std::cout<<" particle "<<(*it).index; fflush(stdout);
			std::cout<<" has position "<<(*it).r<<"\n";
		    }
		}
    }
	
};

#endif
