/**
@file network.cpp
This file implements all functions defined in header network.h
*/
#include "network.h"
#include "gnuplot_i.hpp"
#include "helper_funcs.h"
#include "crack.h"
#include <fstream>
#include <sstream>
#include <iomanip>
// ----------------------------------------------------------------------- 
/// \brief Default constructor
// ----------------------------------------------------------------------- 
Network::Network() {

	initialized = false;

}

// ----------------------------------------------------------------------- 
/// \brief Copy constuctor for class Network
///
/// \param source --> copies from a Network object source
// ----------------------------------------------------------------------- 
Network::Network(Network const & source) {

	initialized = false;
	copy(source);
	initialized = true;

}

// ----------------------------------------------------------------------- 
/// \brief Constructor that instantiates a Network object according
/// to given \emph{fname} file. Please note that there are no checks made
/// on the GMSH .msh filepath and that responsibility is left to the user.
///
/// \param fname --> takes a filename and instantiates Network object
///
// ----------------------------------------------------------------------- 
Network::Network(string& fname) {

	initialized = false;
	load_network(fname);
	initialized = true;

}

// ----------------------------------------------------------------------- 
/// \brief Destructor for class network
///
// ----------------------------------------------------------------------- 
Network::~Network() {

	clear();

}


// ----------------------------------------------------------------------- 
/// \brief Overloads assignment operator (=) to implement copy constructor
///
// ----------------------------------------------------------------------- 
Network const & Network::operator=(Network const & other) {

	if (this != &other) {
		clear();
		copy(other);
	}

	return *this;
	
}

// ----------------------------------------------------------------------- 
/// \brief Called by destructor of Network objects to clear variables 
/// from memory
///
// ----------------------------------------------------------------------- 
void Network::clear() {

	free(R);
	R = NULL;
	free(edges);
	edges = NULL;
	free(forces);
	forces = NULL;
	free(damage);
	damage = NULL;
	free(L);
	L = NULL;
	free(PBC);
	PBC = NULL;
	free(lsideNodes);
	lsideNodes = NULL;
	free(rsideNodes);
	rsideNodes = NULL;
	free(tsideNodes);
	tsideNodes = NULL;
	free(bsideNodes);
	bsideNodes = NULL;
	free(moving_nodes);
	moving_nodes = NULL;

}

// ----------------------------------------------------------------------- 
/// \brief This function is for users to build custom networks. It has no 
/// body and implements nothing. 
///
// ----------------------------------------------------------------------- 
void Network::build_network() {
	//TODO: Add small-network code
	
}

// ----------------------------------------------------------------------- 
/// \brief This function allocates memory for all class data members. 
/// This is called by the Network constructor if the same is called using
/// a filepath string. User is discouraged from directly calling this function.
/// 
/// \param fname (string) --> string of the .msh filepath
/// 
// -----------------------------------------------------------------------
void Network::malloc_network(string& fname){
	
	read_n(n_nodes, n_elems, fname);

	int max_nodes_on_a_side = int(sqrt(n_nodes)*2.0);
	// add memory for side connections
	n_elems += 3*max_nodes_on_a_side;

	size_t sf = sizeof(float);
	size_t si = sizeof(int);
	size_t sb = sizeof(bool);

	/** modified by Dihan
	 *	alloacte memory for potential additional chains for arrays storing edges, forces, L, PBC, damage
	 */
	R = (float*)malloc(n_nodes*DIM*sf);
	edges = (int*)malloc((RANDOM_LONG+RANDOM_Y+n_elems)*2*si);
	forces = (float*)malloc(n_nodes*DIM*sf);
	damage = (float*)malloc((RANDOM_LONG+RANDOM_Y+n_elems)*sf);
	L = (float* )malloc((RANDOM_LONG+RANDOM_Y+n_elems)*sf);
	PBC = (bool* )malloc((RANDOM_LONG+RANDOM_Y+n_elems)*sb);

	// 	initialise n_xside for side nodes
	n_rside = 0;
	n_lside = 0;
	n_bside = 0;
	n_tside = 0;

	// adjust n_elems back to actual
	n_elems -= 3*max_nodes_on_a_side;
}

// ----------------------------------------------------------------------- 
/// \brief Gives the boundary nodes for a network.
///
// -----------------------------------------------------------------------
void Network::side_nodes(){

	// make the list of side nodes
	for(int i=0; i<n_nodes; i++){
		if(fabs(R[i*DIM]-0.0)<TOL){n_lside++;}
		if(fabs(R[i*DIM]-MAXBOUND_X)<TOL){n_rside++;}
		if(fabs(R[i*DIM + 1]-0.0)<TOL){n_bside++;}
		if(fabs(R[i*DIM + 1]-MAXBOUND_Y)<TOL){n_tside++;}
	}

	size_t si = sizeof(int);
	lsideNodes = (int* )malloc(n_lside*si);
	rsideNodes = (int* )malloc(n_rside*si);
	bsideNodes = (int* )malloc(n_bside*si);
	tsideNodes = (int* )malloc(n_tside*si);
	
	n_lside = 0; n_rside = 0; n_bside = 0; n_tside=0;

	for(int i=0; i<n_nodes; i++){
		if(fabs(R[i*DIM]-0.0)<TOL){lsideNodes[n_lside]=i; n_lside++;}
		if(fabs(R[i*DIM]-MAXBOUND_X)<TOL){rsideNodes[n_rside]=i; n_rside++;}
		if(fabs(R[i*DIM + 1]-0.0)<TOL){bsideNodes[n_bside]=i; n_bside++;}
		if(fabs(R[i*DIM + 1]-MAXBOUND_Y)<TOL){tsideNodes[n_tside]=i; n_tside++;}
	}
}

// ----------------------------------------------------------------------- 
/// \brief Removes duplicate edges. A much more efficient version can be made
/// of this program. But it is called only once in the simulation, so this 
/// should be fine.
// -----------------------------------------------------------------------
void Network::remove_duplicates(int& n_elems){
	int counter = 0;
	int node1, node2;
	int d_node1, d_node2, i,k;

	for(i=0;i<n_elems-1;i++){
		node1 = edges[2*i];
		node2 = edges[2*i+1];
		if(node1!=-2 && node2!=-2){
			for(int j=i+1; j<n_elems; j++){
				d_node1 = edges[2*j];
				d_node2 = edges[2*j+1];
				if((node1==d_node1 && node2==d_node2) || (node2==d_node1 && node1==d_node2)){
					edges[2*j] = -2; // -1 token reserved for broken edges. -2 used for readable code
					edges[2*j+1] = -2; // -1 token reserved for broken edges. -2 used for readable code
				}
			}
			counter++;
		}
	}
	// cout<<counter<<endl;
	int* buffer = new int[2*counter];
	for(i=0, k=0;i<n_elems;i++){
		node1 = edges[2*i];
		node2 = edges[2*i+1];
		if(node1!=-2 && node2!=-2){
			buffer[2*k] = node1;
			buffer[2*k+1] = node2;
			k++;
		}
	}

	for(i=0, k=0;i<n_elems;i++){
		if(k<counter){
			edges[2*i] = buffer[2*k];
			edges[2*i+1] = buffer[2*k+1];
			k++;
		}
		else{
			edges[2*i] = -1;
			edges[2*i+1] = -1;
		}
	}
	cout<<n_elems-counter<<" duplicates removed!\n";
	n_elems = counter;
	delete[] buffer;
	buffer = NULL;
}



// ----------------------------------------------------------------------- 
/// \brief Adds random long range edges to the network. 
///	\param n_add (int) --> specify how many links will be added
///	a new log file called add_link_tracker.txt will firstly be generated
/// add_link_tracker.txt will firstly record the average node-node distance and contour length of original network
/// add_link_tracker.txt will then record end node coordinate of each additional chain, its node-node distance and its contour length
/// the number of additional chains will be controlled under half of original regular chains
/// the node-node distance of additional chains will be 20 times longer than regular chains
/// the prestretch of additional chain will be set acoording to prestretch-rate defined by me (see param.h)
/// 
// -----------------------------------------------------------------------

void Network::add_long_range_egdes_random(int n_add, string folder_name){
	if (n_add == 0){
		return;
	}
	cout<<"Asked to add "<<n_add<<" more edges!\n";
	int node1, node2;
	float s = 0;

	if(n_add >= n_elems/2){
		cout<<"You're asking for too many long range bonds,"
				" can't do it! I will add none and continue..."<<endl;
	}
	string fname = std::string(FLDR_STRING)+folder_name + "/add_link_tracker.txt";
	ofstream logger;
	logger.open(fname, ios::trunc|ios_base::out);
	logger<< this->meanX <<"\t";
	logger<< (this->meanX/this->meanXL);
	logger<<"\n";
	srand(time(NULL));

	double pre_str_ct = 0;

	for(int i = 0; i < n_add; i++){
		while(s < 20*meanX){ 
			node1 = rand()%(n_nodes - 4) + 4;
			node2 = rand()%(n_nodes - 4) + 4;
			s = dist(&R[node1*DIM], &R[node2*DIM]);
			if (ismember(node1, lsideNodes, n_lside) || ismember(node1, rsideNodes, n_rside) || ismember(node1, tsideNodes, n_tside) || ismember(node1, bsideNodes, n_bside)){
				s = 0;
			}
			if (ismember(node2, lsideNodes, n_lside) || ismember(node2, rsideNodes, n_rside) || ismember(node2, tsideNodes, n_tside) || ismember(node2, bsideNodes, n_bside)){
				s = 0;
			}
		}

		edges[n_elems*2] = node1;
		edges[n_elems*2 + 1] = node2;
		L[n_elems] = ((1-PRESTRETCH)*s + PRESTRETCH*meanX)/meanXL;
		pre_str_ct += s/L[n_elems];

		logger<< node1 << ": " <<"\t";
		logger<< R[node1*DIM + 0] <<"\t";
		logger<< R[node1*DIM + 1] <<"\t";
		logger<< node2 << ": " <<"\t";
		logger<< R[node2*DIM + 0] <<"\t";
		logger<< R[node2*DIM + 1] <<"\t";
		logger<< s <<"\t";
		logger<< L[n_elems];
		logger<<"\n";

		// update PBC;
		PBC[n_elems] =  false;

		// update damage
		damage[n_elems] = 0;
		n_elems++;

		s = 0;
	}

    logger.close();
    pre_str_ct /= n_add;
    cout << "average additional link prestretch is "<<pre_str_ct<< endl;
	cout<<"Stored long link info in "<<fname<<"!\n";

}

// ----------------------------------------------------------------------- 
/// \brief Adds y-directional long range edges to the network. 
///	Basically same as the random direction function. 
/// Note that the orientation is not strictly vertical, since the mesh grids are not uniform
/// the horizontal coordination difference between two end is constrained by 3 times of node-node distance of regular chains
// -----------------------------------------------------------------------

void Network::add_long_range_egdes_y(int n_add, string folder_name){
	if (n_add == 0){
		return;
	}
	cout<<"Asked to add "<<n_add<<" more edges!\n";
	int node1, node2;
	float s = 0;
	float x_diff = MAXBOUND_X/2;
	if(n_add >= n_elems/2){
		cout<<"You're asking for too many long range bonds,"
				" can't do it! I will add none and continue..."<<endl;
	}
	string fname = std::string(FLDR_STRING)+ folder_name + "/add_link_tracker.txt";
	ofstream logger;
	logger.open(fname, ios::trunc|ios_base::out);
	logger<< this->meanX <<"\t";
	logger<< (this->meanX/this->meanXL);
	logger<<"\n";
	srand(time(NULL));

	double pre_str_ct = 0;

	for(int i = 0; i < n_add; i++){
		// while(s<L_MEAN/4 || s>L_MEAN){ 
		while(s<20*this->meanX || x_diff>3*meanX){ 
			node1 = rand()%(n_nodes - 4) + 4;
			node2 = rand()%(n_nodes - 4) + 4;
		
			s = dist(&R[node1*DIM], &R[node2*DIM]);
			x_diff = abs(R[node1*DIM]-R[node2*DIM]);
		}

		// add nodes to edges array
		edges[n_elems*2] = node1;
		edges[n_elems*2 + 1] = node2;


		float x1 = R[node1*DIM+0];
		float y1 = R[node1*DIM+1];
		float x2 = R[node2*DIM+0];
		float y2 = R[node2*DIM+1];

		L[n_elems] = ((1-PRESTRETCH)*s + PRESTRETCH*meanX)/meanXL;
		pre_str_ct += s/L[n_elems];

		// log file:
		logger<< R[node1*DIM + 0] <<"\t";
		logger<< R[node1*DIM + 1] <<"\t";
		logger<< R[node2*DIM + 0] <<"\t";
		logger<< R[node2*DIM + 1] <<"\t";
		logger<< s <<"\t";
		logger<< L[n_elems];
		logger<<"\n";
		PBC[n_elems] =  false;
		damage[n_elems] = 0;
		n_elems++;

		s = 0;
	}
    logger.close();
    pre_str_ct /= n_add;
    cout << "average additional link prestretch is "<<pre_str_ct<< endl;
	cout<<"Stored long link info in "<<fname<<"!\n";

}


// ----------------------------------------------------------------------- 
/// \brief Create patterned regions on this network
///	\param type (string) --> indicates the pattern type is "layer" or "spot" on this network
///	\param region_number (int) --> indicate how many patterned region is demanded for certain type on this network
///	\param rate (float) --> indicates how many times sparser the patterned region will be
// -----------------------------------------------------------------------
void Network::patterning(string type, int region_number, double rate){
	if (type == "none"){
		return;
	}

	if (type == "layer_h"){
		// int n_elems_per_layer = n_elems/(2*region_number+1);
		double unit_thick = MAXBOUND_Y/(region_number+(region_number+1)*rate);
		cout << "divide to "<<MAXBOUND_Y/unit_thick<<" parts"<<endl;

		for (int i=0; i<region_number; i++){// in one region
			for (int e=0; e<n_elems; e++){// check all elems, delete and record

				double cur_ini_y = R[2*edges[2*e+0]+1];
				double cur_end_y = R[2*edges[2*e+1]+1];
				double region_low = (i+1)*rate*unit_thick + i*unit_thick;
				double region_top = region_low+unit_thick;
				// cout <<"bound is "<<region_low<<" and "<<region_top<<endl;
				if ((cur_ini_y>region_low && cur_ini_y<region_top) && (cur_end_y>region_low && cur_end_y<region_top)){
					L[e] *= rate;
				}
			}
		}
	}


	if (type == "layer_v"){
		// int n_elems_per_layer = n_elems/(2*region_number+1);
		double unit_thick = MAXBOUND_X/(region_number+(region_number+1)*rate);
		cout << "divide to "<<MAXBOUND_X/unit_thick<<" parts"<<endl;

		for (int i=0; i<region_number; i++){// in one region
			for (int e=0; e<n_elems; e++){// check all elems, delete and record

				double cur_ini_x = R[2*edges[2*e+0]+0];
				double cur_end_x = R[2*edges[2*e+1]+0];
				double region_left = (i+1)*rate*unit_thick + i*unit_thick;
				double region_right = region_left+unit_thick;
				// cout <<"bound is "<<region_low<<" and "<<region_top<<endl;
				if ((cur_ini_x>region_left && cur_ini_x<region_right) && (cur_end_x>region_left && cur_end_x<region_right)){
					L[e] *= rate;
				}
			}
		}
	}
	if (type == "spot"){
		double gridX = MAXBOUND_X/(1+region_number);
		double gridY = MAXBOUND_Y/(1+region_number);
		double radius;
		// 0.382 is sqaure of golden ratio
		if (gridX > gridY){
			radius = 0.382*gridY;
		}else{
			radius = 0.382*gridX;
		}
		double center_x, center_y;
		int idx, x_idx, y_idx;
		// srand(time(NULL));
		for (int i=0; i<region_number*region_number; i++){
			// int idx = rand()%(region_number*region_number) + 1;
			// x_idx = floor((idx-1)/region_number + 1);
			// y_idx = idx%region_number + 1;
			center_x = (1+i%region_number)*gridX;
			center_y = (1+i/int(region_number))*gridY;

			for (int e=0; e<n_elems; e++){
				double cur_ini_x = R[2*edges[2*e+0]+0];
				double cur_ini_y = R[2*edges[2*e+0]+1];
				double cur_end_x = R[2*edges[2*e+1]+0];
				double cur_end_y = R[2*edges[2*e+1]+1];

				double ini_cen = sqrt((cur_ini_x-center_x)*(cur_ini_x-center_x) + (cur_ini_y-center_y)*(cur_ini_y-center_y));
				double end_cen = sqrt((cur_end_x-center_x)*(cur_end_x-center_x) + (cur_end_y-center_y)*(cur_end_y-center_y));
				if ((ini_cen < radius) && (end_cen < radius)){

					L[e] *= rate;

				}
			}

		}
	}

}

// void Network::long_crack(float x_start, float x_end, float y_start, float y_end){
// 	if (!MESH_CRACK){
// 		return;
// 	}
// 	int c_break = 0;
// 	for (int i=0; i<n_elems; i++){
// 		int node1 = edges[i*2+0];
// 		float x1 = R[node1*DIM+0];
// 		float y1 = R[node1*DIM+1];
// 		int node2 = edges[i*2+1];
// 		float x2 = R[node2*DIM+0];
// 		float y2 = R[node2*DIM+1];
// 		float x_mid = (x1+x2)/2;
// 		float y_mid = (y1+y2)/2;
// 		if (y_start<y_mid && y_end>y_mid && x_start<x_mid && x_end>=x_mid){
// 			edges[i*2+0] = -1;
// 			edges[i*2+1] = -1;
// 			c_break += 1;
// 		}
// 	}
// 	cout << "broke "<<c_break<<" links to create a long crack."<<endl;
// }

// ----------------------------------------------------------------------- 
/// \brief This function builds the network. It includes building PBC 
/// bonds. The function initializes all properties required for running
/// experiments on the instantiated Network code. This is called by constructor
/// of a Network object. User is discouraged from directly calling this function.
///
/// \param fname (string) --> string of the .msh filepath
// -----------------------------------------------------------------------
void Network::load_network(string& fname) {

	if (initialized) {
		clear();
	}
	//Check if file exists
	bool exists = does_file_exist(fname);
	if(!exists){
		cout<<"File does not exist!\n";
		return;
	}
	
	//Malloc all variables
	malloc_network(fname);

	cout<<"Malloc was successful!\n";

	cout<<"Reading the mesh...\n";
	take_input(R, edges, n_nodes, n_elems, fname);

	// remove duplicate edges
	remove_duplicates(n_elems);


	cout<<"Mesh read successfully!\n";
	cout<<"Number of nodes are: "<<n_nodes<<endl;
	cout<<"Number of elements are: "<<n_elems<<endl;


	// initialize boundary node arrays
	side_nodes();


	// moving nodes
	int c = 0;
	n_moving = n_nodes - n_tside - n_bside;
	moving_nodes = (int*)malloc(sizeof(int)*n_moving);
	for(int i =0; i<n_nodes; i++){
		if(!ismember(i, tsideNodes, n_tside) && !ismember(i, bsideNodes, n_bside)){
			moving_nodes[c] = i;
			c++;
		}
	
	}
	
	

	cout<<"We have "<<c<<" moving nodes\n";


	cout<<"Side nodes written successfully! \n";
	__init__(L, damage, PBC, n_elems);

	if(IMPLEMENT_PBC){
		this->make_edge_connections(15.0);
		cout<<"Number after new connections made: "<<n_elems<<endl;
	}
	cout<<__LINE__<<endl;
}

// ----------------------------------------------------------------------- 
/// \brief Implements the body of the copy constuctor. Allocates memory and
/// does the appropriate assignments
///
/// \param source (Network) --> an initialized object of type Network
// -----------------------------------------------------------------------
void Network::copy(Network const & source) {

	cracked = source.cracked;
	n_nodes = source.n_nodes;
	n_elems = source.n_elems;
	n_rside = source.n_rside;
	n_lside = source.n_lside;
	n_bside = source.n_bside;
	n_tside = source.n_tside;
	n_moving = source.n_moving;

	size_t sf = sizeof(float);
	size_t si = sizeof(int);
	size_t sb = sizeof(bool);

	R = (float*)malloc(n_nodes*DIM*sf);
	edges = (int*)malloc((RANDOM_LONG+RANDOM_Y+n_elems)*2*si);
	forces = (float*)malloc(n_nodes*DIM*sf);
	damage = (float*)malloc((RANDOM_LONG+RANDOM_Y+n_elems)*sf);
	L = (float* )malloc((RANDOM_LONG+RANDOM_Y+n_elems)*sf);
	PBC = (bool* )malloc((RANDOM_LONG+RANDOM_Y+n_elems)*sb);
	lsideNodes = (int* )malloc(n_lside*si);
	rsideNodes = (int* )malloc(n_rside*si);
	bsideNodes = (int* )malloc(n_bside*si);
	tsideNodes = (int* )malloc(n_tside*si);
	moving_nodes = (int *)malloc(n_moving*si);

	for (int i = 0; i < n_nodes * DIM; i++) {
		R[i] = source.R[i];
		forces[i] = source.forces[i];
	}
	
	for (int i = 0; i < (RANDOM_LONG+RANDOM_Y+n_elems) * 2; i++) {
		edges[i] = source.edges[i];
	}
	
	for (int i = 0; i < n_lside; i++) {
		lsideNodes[i] = source.lsideNodes[i];
	}
	for (int i = 0; i < n_rside; i++) {
		rsideNodes[i] = source.rsideNodes[i];
	}
	for (int i = 0; i < n_tside; i++) {
		tsideNodes[i] = source.tsideNodes[i];
	}
	for (int i = 0; i < n_bside; i++) {
		bsideNodes[i] = source.bsideNodes[i];
	}
	for (int i = 0; i < n_moving; i++) {
		moving_nodes[i] = source.moving_nodes[i];
	}

	for (int i = 0; i < (RANDOM_LONG+RANDOM_Y+n_elems); i++) {
		damage[i] = source.damage[i];
		L[i] = source.L[i];
		PBC[i] = source.PBC[i];
	}

}



// ----------------------------------------------------------------------- 
/// \brief Adds PBC bonds to the network according to dely_allowed.
///
/// \param dely_allowed(float) --> Connects two nodes i and j on the two 
/// sides of the network respectively iff $abs(y_i - y_j) < $dely_allowed
// -----------------------------------------------------------------------
void Network::make_edge_connections(float dely_allowed) {

	std::default_random_engine seed;
	std::uniform_real_distribution<float> generator(L_MEAN - L_STD, L_MEAN + L_STD);
	int nl, nr, lnode, rnode;
	for(nl= 0; nl < n_lside; nl++){
		lnode = lsideNodes[nl];
		for(nr= 0; nr < n_rside; nr++){
			rnode = rsideNodes[nr];
			if (fabs(R[lnode*DIM + 1] - R[rnode*DIM + 1]) < dely_allowed){
				edges[n_elems*2] = rnode;
				edges[n_elems*2 + 1] = lnode;
				//cout<<"Connected node "<<lnode<<" and "<<rnode<<"\n";
				L[n_elems] = generator(seed);
				damage[n_elems] = 0.0;
				PBC[n_elems] = true;
				n_elems += 1;
			}
		}
	}

}


// ----------------------------------------------------------------------- 
/// \brief Adds random cracks in the Network object. This is implemented 
/// by removing all connections and nodes that lie \emph{entirely} within
/// the ellipses specified in alist. Please note this assumption may be
/// relaxed and a probabilistic measure may be introduced to remove some
/// connections so as to simulate small world networks (cluster networks)
/// This is managed by the PROB_REMOVAL parameter in params.h 
///
/// \param alist (Cracklist) --> A list of ellipses
// -----------------------------------------------------------------------
void Network::apply_crack(Cracklist & alist) {

	cracked = true;
	std::default_random_engine seed;
	std::uniform_real_distribution<float> generator(0, 1);
	float equation = 0;
	//make clusters
	//vector<int> within_circle;
	vector<int> nodes_to_remove;
	int edges_removed = 0; 
	int node1, node2;
	Crack crack;
	float x;
	float y;
	float si,co;
	for(int ci = 0; ci < alist.n_cracks; ci++){
		crack.setter(alist[ci]);
		for(int i=0; i<n_nodes; i++){
			x = R[i*DIM];
			y = R[i*DIM+1];
			x -= crack.c[0];
			y -= crack.c[1];
			si = crack.trig[0];
			co = crack.trig[1];
			equation = pow((x*co + y*si)/crack.a[0], 2)+ \
						pow((x*si - y*co)/crack.a[1], 2);

			equation -= 1.0;

			if(equation<=0.0){
				nodes_to_remove.push_back(i);
			}
		}
	}

	// for(int i= 0; i < n_nodes; i++){
	// 	if(contains(within_circle, i)){
	// 		nodes_to_remove.push_back(i);
	// 	}
	// }

	for(int i = 0; i<n_elems; i++){
		node1 = edges[2*i];
		node2 = edges[2*i+ 1];
		if(contains(nodes_to_remove, node1) && contains(nodes_to_remove, node2)){
			if(generator(seed) < PROB_REMOVAL){
				edges[2*i] = -1;
				edges[2*i + 1] = -1;
				edges_removed += 1;
			}
		}
	}
	cout<<"Edges removed : "<<edges_removed<<endl;

}



// ----------------------------------------------------------------------- 
/// \brief Helper function that calls gnuplot to output a png file
///	\param png_size (string) --> the string in format "x_resolution,y_resolution" indicating the output figure size
/// \param png_path (string) --> the location+file name
/// \param xrange (string) --> [xmin:xmax]
/// \param fname (string) --> png_data.txt
/// \param pic_name (string) --> file name of png
// -----------------------------------------------------------------------
void pngPlotHelper(string png_size, string png_path, string xrange, string fname, string pic_name){
	Gnuplot gnu;
	gnu.cmd("set xrange " + xrange);
	gnu.cmd("load 'viridis.pal'");
	gnu.cmd("set key off");
	gnu.cmd("set colorbox default vertical");
	gnu.cmd("set cbrange [0:1]");
	gnu.cmd("set cbtics ('0.0' 0.0,'1.0' 1.0) offset 0,0.5 scale 0");
	gnu.cmd("set size ratio -1"); // to have the same scale on x and y axes
	gnu.cmd("plot '" + fname + "' every ::0::1 u 1:2:3 w l lc palette lw 2 title '"+\
		pic_name+"'");
	gnu.cmd("set term png size " + png_size);
	gnu.cmd("set output '"+ png_path + "'");

	gnu.cmd("replot");
	// gnu.cmd("set term x11");
	gnu.reset_plot();
}



// ----------------------------------------------------------------------- 
/// \brief Helper function that calls gnuplot to output a eps file
///	\param png_size (string) --> the string in format "x_resolution,y_resolution" indicating the output figure size
/// \param png_path (string) --> the location+file name
/// \param xrange (string) --> [xmin:xmax]
/// \param fname (string) --> eps_data.txt
/// \param pic_name (string) --> file name of eps
// -----------------------------------------------------------------------
void epsPlotHelper(string eps_size, string eps_path, string xrange, string fname, string pic_name){
	Gnuplot gnu;
	gnu.cmd("set xrange " + xrange);
	gnu.cmd("load 'viridis.pal'");
	gnu.cmd("set key off");
	gnu.cmd("set colorbox default vertical");
	gnu.cmd("set cbrange [0:1]");
	gnu.cmd("set cbtics ('0.0' 0.0,'1.0' 1.0) offset 0,0.5 scale 0");
	gnu.cmd("set size ratio -1"); // to have the same scale on x and y axes
	gnu.cmd("plot '" + fname + "' every ::0::1 u 1:2:3 w l lc palette lw 2 title '"+\
		pic_name+"'");
	// gnu.cmd("set term png size " + png_size);
	// gnu.cmd("set terminal postscript eps enhanced size "+png_size);
	gnu.cmd("set terminal postscript eps size "+eps_size);
	gnu.cmd("set output '"+ eps_path + "'");
	gnu.cmd("replot");
	// gnu.cmd("set term x11");
	gnu.reset_plot();
}

// ----------------------------------------------------------------------- 
/// \brief Plots the graph of the network. Note the edges in the graph just
/// signify connection between two nodes and do not represent the actual 
/// polymer geometry in any sense. Note that this plot will be saved as .eps 
/// file in the folder specified by FLDR_STRING in params.h
///
/// \param iter_step (int) --> The iteration number of the experiment. This
/// number will also be the filename for the .eps file
///
/// \param first_time (bool) --> This just effects what quantity is mapped
/// to the color of the graph. If True, colors represent magnitude of forces
/// in the polymer else, they represent the damage in a polymer connection
// -----------------------------------------------------------------------
void Network::plotNetwork(int iter_step, bool first_time, string folder_name){
	if (EPS == 0){
		return;
	}
	ofstream f;
	Gnuplot gnu;
	string fname = std::string(FLDR_STRING)+folder_name + "/" + "eps_data.txt";
	float c;
	f.open(fname,std::ofstream::out | std::ofstream::trunc);
	// std::default_random_engine seed;
	// std::uniform_real_distribution<float> arbitcolor(0.0,1.0);
	int node1, node2;
	if(first_time){
		float r1[DIM];
		float r2[DIM];
		float s;
		int j, k;
		for (j = 0; j < n_elems; j++){
		// read the two points that form the edge // 2 because 2 points make an edge! Duh.
			node1 = edges[j * 2]; 
			node2 = edges[j * 2 + 1];
			
			// // check if pair exists
			// if(node1 == -1 || node2 == -1) {
			// 	continue;
			// }

			// read the positions
			#pragma unroll
			for(k = 0; k<DIM; k++){
				r1[k] = R[node1*DIM + k]; 
				r2[k] = R[node2*DIM + k];
			}
			
			// check PBC_STATUS
			if (PBC[j]) {
				// add PBC_vector to get new node position
				#pragma unroll
				for (k = 0; k < DIM; k++){
					r2[k] += PBC_vector[k];
				}
				// get force on node1 due to node2
				s = dist(r1, r2);
			}
			else{
				s = dist(r1, r2);
			}
		
			c = s/L[j];
			for(int d = 0; d<DIM; d++){
					f<<R[node1*DIM+d]<<"\t";
				}
				f<<c<<endl;
			if(!PBC[j]){
				for(int d = 0; d<DIM; d++){
					f<<R[node2*DIM+d]<<"\t";
				}
			}
			else{
				for(int d = 0; d<DIM; d++){
					f<<(R[node1*DIM+d]+10)<<"\t";
				}				
			}
			f<<c<<endl<<endl;
		}
		f.close();
	}
	else{
	for(int i = 0; i<n_elems; i++){
		c = damage[i];
		node1 = edges[2*i];
		node2 = edges[2*i+1];

		if(node1!=-1 && node2!=-1){
			for(int d = 0; d<DIM; d++){
					f<<R[node1*DIM+d]<<"\t";
				}
				f<<c<<endl;
			if(!PBC[i]){
				for(int d = 0; d<DIM; d++){
					f<<R[node2*DIM+d]<<"\t";
				}
			}
			else{
				for(int d = 0; d<DIM; d++){
					f<<(R[node1*DIM+d]+10)<<"\t";
				}				
			}
		f<<c<<endl<<endl;
		}
	}
	f.close();
	}
	//Plot to h
	float aspect_ratio = (MAXBOUND_X)/(R[tsideNodes[0]*DIM+1]); 

	stringstream convert;

	convert<<"[-5:"<<MAXBOUND_X+5<<"]";
	string xrange = convert.str();
	convert.str(std::string());
	
	int x_res = 1280;
	convert<<x_res<<", "<<int((x_res-410)/aspect_ratio+410);
	string png_size = convert.str();
	convert.str(std::string());

	int eps_x_res = 40;
	int eps_y_res = eps_x_res/aspect_ratio;
	convert<<eps_x_res<<", "<<eps_y_res;
	string eps_size = convert.str();
	convert.str(std::string());

	std::stringstream ss;
	ss << std::setw(5) << std::setfill('0') << iter_step;
	std::string pic_name = ss.str();
	string eps_path = std::string(FLDR_STRING)+folder_name + "/graphs/" + pic_name + ".eps";
	epsPlotHelper(eps_size, eps_path, xrange, fname, pic_name);
	
}




// ----------------------------------------------------------------------- 
/// \brief Plots the graph of the network. Note the edges in the graph just
/// signify connection between two nodes and do not represent the actual 
/// polymer geometry in any sense. Note that this plot will be saved as .png 
/// file in the folder specified by FLDR_STRING in params.h
///
/// \param iter_step (int) --> The iteration number of the experiment. This
/// number will also be the filename for the .png file
///
/// \param first_time (bool) --> This just effects what quantity is mapped
/// to the color of the graph. If True, colors represent magnitude of forces
/// in the polymer else, they represent the damage in a polymer connection
// -----------------------------------------------------------------------
void Network::plotFrames(int iter_step, bool first_time, string folder_name){
	if (PNG == 0){
		return;
	}
	ofstream f;
	Gnuplot gnu;
	string fname = std::string(FLDR_STRING)+folder_name + "/" + "frame_data.txt";
	float c;
	f.open(fname,std::ofstream::out | std::ofstream::trunc);

	int node1, node2;
	if(first_time){
		float r1[DIM];
		float r2[DIM];
		float s;
		int j, k;
		for (j = 0; j < n_elems; j++){
		// read the two points that form the edge // 2 because 2 points make an edge! Duh.
			node1 = edges[j * 2]; 
			node2 = edges[j * 2 + 1];
			
			// read the positions
			#pragma unroll
			for(k = 0; k<DIM; k++){
				r1[k] = R[node1*DIM + k]; 
				r2[k] = R[node2*DIM + k];
			}
			
			// check PBC_STATUS
			if (PBC[j]) {
				// add PBC_vector to get new node position
				#pragma unroll
				for (k = 0; k < DIM; k++){
					r2[k] += PBC_vector[k];
				}
				// get force on node1 due to node2
				s = dist(r1, r2);
			}
			else{
				s = dist(r1, r2);
			}
		
			c = s/L[j];
			for(int d = 0; d<DIM; d++){
					f<<R[node1*DIM+d]<<"\t";
				}
				f<<c<<endl;
			if(!PBC[j]){
				for(int d = 0; d<DIM; d++){
					f<<R[node2*DIM+d]<<"\t";
				}
			}
			else{
				for(int d = 0; d<DIM; d++){
					f<<(R[node1*DIM+d]+10)<<"\t";
				}				
			}
			f<<c<<endl<<endl;
		}
		f.close();
	}
	else{
	for(int i = 0; i<n_elems; i++){
		c = damage[i];
		node1 = edges[2*i];
		node2 = edges[2*i+1];

		if(node1!=-1 && node2!=-1){
			for(int d = 0; d<DIM; d++){
					f<<R[node1*DIM+d]<<"\t";
				}
				f<<c<<endl;
			if(!PBC[i]){
				for(int d = 0; d<DIM; d++){
					f<<R[node2*DIM+d]<<"\t";
				}
			}
			else{
				for(int d = 0; d<DIM; d++){
					f<<(R[node1*DIM+d]+10)<<"\t";
				}				
			}
		f<<c<<endl<<endl;
		}
	}
	f.close();
	}
	//Plot to h
	float aspect_ratio = (MAXBOUND_X)/(R[tsideNodes[0]*DIM+1]); 

	stringstream convert;

	convert<<"[-5:"<<MAXBOUND_X+5<<"]";
	string xrange = convert.str();
	convert.str(std::string());
	
	int x_res = 1280;
	convert<<x_res<<", "<<int((x_res-410)/aspect_ratio+410);
	string png_size = convert.str();
	convert.str(std::string());


	int eps_x_res = 40;
	int eps_y_res = eps_x_res/aspect_ratio;
	convert<<eps_x_res<<", "<<eps_y_res;
	string eps_size = convert.str();
	convert.str(std::string());

	std::stringstream ss;
	ss << std::setw(5) << std::setfill('0') << int(iter_step/PNG);
	std::string pic_name = ss.str();
		string png_path = std::string(FLDR_STRING)+folder_name + "/frames/" + pic_name + ".png";
		pngPlotHelper(png_size, png_path, xrange, fname, pic_name);
	
}
// ----------------------------------------------------------------------- 
/// \brief Gets the active (load-bearing) edges in the Network object. 
/// Also prints out that info to STDOUT or piped OUT
///
/// 
// -----------------------------------------------------------------------
int Network::get_current_edges(){
	int n_elems_current = 0;
	for(int i=0; i<n_elems; i++){
		if(edges[i*2] != -1 && edges[i*2 + 1]!=-1){
			n_elems_current += 1;
		}
	}
	return n_elems_current;
}


// ----------------------------------------------------------------------- 
/// \brief record the number of remaining chains at this iteration in remain_chains array
///	\param remain_chains (int*) --> An array storing number of remaining chains at each iteration
/// \param iter (int) --> the iteration index
/// \param curr_n_edges (int) --> the number of remaining chains at this iteration
// -----------------------------------------------------------------------
void Network::get_edge_number(int* remain_chains, int iter, int curr_n_edges){
	remain_chains[iter] = curr_n_edges;
}

// ----------------------------------------------------------------------- 
/// \brief Prints out certain statistics of the network object.
///
/// 
// -----------------------------------------------------------------------
bool Network::get_stats(){
	/** add by Dihan
	 *	add 1 more features to this function
	 *	calculate the overall orientation parameter
	 *	average cos(theta)^2 for all elem (might need a helper function?)
	 *	OP = 0.5*(3*(average of square)-1)
	 */
	int node1, node2;
	int j, k, id; // loop variables
	float r1[DIM]; float r2[DIM] ;
	float mean_x = 0.0;
	float mean_t = 0.0;
	float var_x = 0.0;
	float var_t = 0.0;
	float max_x = 0.0;
	float max_t = 0.0;
	float s;
	int short_L = 0;
	float shortage = 0.0;
	int c = 0;

	if(n_elems <= 0){
		cout<<"Oops! Something's wrong! n_elems is not positive!"<<endl;
	}
	for (j = 0; j < n_elems; j++){
		// read the two points that form the edge // 2 because 2 points make an edge! Duh.
		node1 = edges[j * 2]; 
		node2 = edges[j * 2 + 1];
		
		// check if pair exists
		if(node1 == -1 || node2 == -1) {
			continue;
		}
		
		c++;


		// read the positions
		#pragma unroll
		for(k = 0; k<DIM; k++){
			r1[k] = R[node1*DIM + k]; 
			r2[k] = R[node2*DIM + k];
		}
		
		// check PBC_STATUS
		if (PBC[j]) {
			// add PBC_vector to get new node position
			#pragma unroll
			for (k = 0; k < DIM; k++){
				r2[k] += PBC_vector[k];
			}
			// check if L too short
			s = dist(r1, r2);
			mean_x += s;
			mean_t += s/L[j];
			if(s>max_x){max_x =s;}
			if(s/L[j]>max_t){max_t=s/L[j];}
			var_x += s*s;
			var_t += s*s/L[j]/L[j];
			if(s >= 0.99*L[j]){
				short_L++;
				shortage = s - L[j];
			}
		}
		else{
			s = dist(r1, r2);
			mean_x += s;
			mean_t += s/L[j];
			if(s>max_x){max_x =s;}
			if(s/L[j]>max_t){max_t=s/L[j];}
			var_x += s*s;
			var_t += s*s/L[j]/L[j];
			if(s >= 0.99*L[j]){
				short_L++;
				shortage = s - L[j];
			}
		}
	}
	if(c<=0){
		cout<<"No load bearing edges remain! Stopping simulation!\n";
		return true;
	}

	mean_x /= c;
	var_x /= c;
	var_x -= mean_x*mean_x;

	mean_t /= c;
	var_t /= c;
	var_t -= mean_t*mean_t;
	// cout<<"Number of edges: "<<c<<endl;
	// cout<<"Number of edges with x>0.99L: "<<short_L<<endl;
	// cout<<"node-node distance mean: "<<mean_x<<endl;
	this->meanX = mean_x;
	// cout<<"node-node distance std_dev: "<<sqrt(var_x)<<endl;
	// cout<<"node-node distance max: "<<max_x<<endl;
	cout<<"\n";
	// cout<<"node-node x/L mean: "<<mean_t<<endl;
	this->meanXL = mean_t;
	// cout<<"node-node x/L std_dev: "<<sqrt(var_t)<<endl;	
	// cout<<"node-node x/L max: "<<max_t<<endl;


	if((float)c/(float)n_elems < 0.02){
		cout<<"Too few edges remain in given mesh! Exiting...\n";
		return true;
	}
	return false;
}


// ----------------------------------------------------------------------- 
/// \brief Gets the total length of all polymer chains active in the network
/// 
// -----------------------------------------------------------------------
float Network::get_weight(){
	float sum_length = 0.0;
	for(int i = 0; i < n_elems; i++){
		if(edges[i*2] != -1 && edges[i*2 + 1]!=-1){
			sum_length += L[i];
		}
	}
	cout<<"\nNetwork weight is "<<sum_length<<"\n";
	return sum_length;
}

// ----------------------------------------------------------------------- 
/// \brief Sets the total length of all polymer chains according to supplied
/// weight parameter. Also prints out the new average polymer contour lengths
/// 
/// \param weight (float) --> Weight goal required.
// -----------------------------------------------------------------------
float Network::set_weight(float weight){
	float curr_weight = get_weight();
	float alpha = weight/curr_weight;
	for(int i = 0; i < n_elems; i++){
		if(edges[i*2] != -1 && edges[i*2 + 1]!=-1){
			L[i] *= alpha;
		}
	}
	cout<<"\n Network weight reset to "<<weight<<"\n";
	cout<<"Average L now is : "<<weight/n_elems<<"\n";
	return alpha;
}


// ----------------------------------------------------------------------- 
/// \brief Moves the top plate according to velocity supplied in vel.h
/// 
///
// -----------------------------------------------------------------------
void Network::move_top_plate(){
	int node;

	for(int i = 0; i<n_tside; i++){
		node = tsideNodes[i];
		#pragma unroll
		for(int d=0; d<DIM; d++){
			R[node*DIM + d] += TIME_STEP*vel[d]; 
		}
	}
}


// ----------------------------------------------------------------------- 
/// \brief Calculates forces along edges of the polymer network graph (i.e
/// forces in each polymer). The forces are then assigned to nodes which 
/// prepares the network for the optimization step. If update_damage is 
/// true, the appropriate damage metric is also updated. If damage exceeds
/// predefined thresholds, the edge is broken. In case of a sacNetwork object
/// hidden lengths are also opened according to damage in these hidden bonds
/// apart from end damage. The function does not return anything as it updates
/// the objects forces directly.
///
/// \param update_damage (bool) --> flag to update damage
// -----------------------------------------------------------------------
void Network::get_forces(bool update_damage) {
	int node1, node2;
	int j, k, id; // loop variables
	float r1[DIM]; float r2[DIM] ;
	float edge_force[DIM];
	float rhat[DIM];
	float s;
	float force;


	memset(forces, 0.0, n_nodes*DIM*sizeof(*forces));

	for (j = 0; j < n_elems; j++){
		// read the two points that form the edge // 2 because 2 points make an edge! Duh.
		node1 = edges[j * 2]; 
		node2 = edges[j * 2 + 1];
		
		// check if pair exists
		if(node1 == -1 || node2 == -1) {
			continue;
		}

		// read the positions
		#pragma unroll
		for(k = 0; k<DIM; k++){
			r1[k] = R[node1*DIM + k]; 
			r2[k] = R[node2*DIM + k];
		}
		
		// check PBC_STATUS
		if (PBC[j]) {
			// add PBC_vector to get new node position
			#pragma unroll
			for (k = 0; k < DIM; k++){
				r2[k] += PBC_vector[k];
			}
			// get force on node1 due to node2
			s = dist(r1, r2);
			unitvector(rhat, r1, r2);
			force = force_wlc(s, L[j]);
			if(force == 999999){edges[j*2] = -1; edges[j*2 +1] = -1; force =0.0;}
			convert_to_vector(edge_force, force, rhat);
			// subtract back the PBC_vector to get original node position
			// #pragma unroll
			// for (k = 0; k < DIM; k++){
			// 	r2[k] -= PBC_vector[k];
			// }
		}
		else{
			s = dist(r1, r2);

			/** add by Dihan
			 *	Critical edge case:
			 *	when two end node of one element overlap, update the force with a zero vector
			 *	aka not processing the update part.
			 */
			if (s == 0){
				continue;
			}

			unitvector(rhat, r1, r2);
			force = force_wlc(s, L[j]);
			// if(force == 999999){edges[j*2] = -1; edges[j*2 +1] = -1; force =0.0;}
			convert_to_vector(edge_force, force, rhat);
		}

		#pragma unroll
		for (k = 0; k < DIM; k++){
			forces[node1*DIM + k] -= edge_force[k];
			forces[node2*DIM + k] += edge_force[k];
		}

		
		//update damage if needed
		if (update_damage){
			if(RATE_DAMAGE){
				damage[j] += kfe(force)*TIME_STEP;
				if(damage[j] > 1.0){
					cout<<"RATE_DAMAGE: Breaking bond between "
					<<edges[j*2]<<" and "<<edges[2*j +1]<<" F, s/L = "<<force \
					<<", "<<s/L[j]<<"s: "<<s<<", L: "<<L[j]<<endl;
					edges[j*2] = -1; edges[j*2+1] = -1;
				}
			}
			else{
				damage[j] = s/L[j];
				if(damage[j] > 0.9){
					cout<<"Breaking bond between "
					<<edges[j*2]<<" and "<<edges[2*j +1]<<" F, s/L = "<<force \
					<<", "<<s/L[j]<<endl;
					edges[j*2] = -1; edges[j*2+1] = -1;
				}
			}
		}
	}

}


// ----------------------------------------------------------------------- 
/// \brief Equilibriates the Network object using RMSProp algorithm due to 
/// G. Hinton et al (a gradient based force minimization over each node). 
/// Stops when maximum force on node is less than TOL specified in params.h
/// or max_iter is reached
/// 
/// \param eta (float) --> learning rate
/// \param alpha (float) --> history weighting parameter
/// \param max_iter (int) --> maximum iterations to be allowed in descent
// -----------------------------------------------------------------------
void Network::optimize(float eta, float alpha, int max_iter){
	float* rms_history = new float[n_moving*DIM](); // () allows 0.0 initialization
	float g, delR;
	char p;
	int id, d, node;
	for(int step = 0; step < max_iter; step++){
		get_forces(false);
		// plotNetwork(step, true);
		// cin>>p;

		if(getabsmax(forces,n_nodes*DIM)>TOL){	
			for(id = 0; id < n_moving; id++){
				node = moving_nodes[id];
				#pragma unroll
				for(d = 0; d<DIM; d++){

					g = forces[DIM*node+d];
					rms_history[id*DIM + d] = alpha*rms_history[id*DIM + d] + (1-alpha)*g*g;
					delR = sqrt(1.0/(rms_history[id*DIM + d] + TOL))*eta*g;
					R[node*DIM + d] += delR;
				}

				/** add by Dihan
				 *	these lines are used to add roller constrain (lateral displacement) on the network
				 * 	for each interation, the x-coordinate of the node on left and right side bound will be adjust to zero or width of the network
				 */
				if (ROLLER && (ismember(node, lsideNodes, n_lside))) {
					R[node*2] = 0;	
				}
				if (ROLLER && (ismember(node, rsideNodes, n_rside))) {
					R[node*2] = MAXBOUND_X;	
				}	

				/** add by Dihan
				 *	at last step, loop through all moving nodes (we just neglect bottom ones)
				 *	calculate CC for each one (might need a helper function?)
				 *	sum up and take average at last
				 */

			}
		}
		else{
			break;
		}
	}
	//
	get_forces(true);
	delete[] rms_history;
}

// ----------------------------------------------------------------------- 
/// \brief Writes the forces on the top plate the supplied plate_forces array
/// 
/// \param plate_forces (float*) --> storage array for plate_forces
/// \param iter (int) --> force_{i} stored in plate_forces[iter*DIM + i]
/// i = 0 ... DIM - 1
/// \param max_iter (int) --> maximum iterations to be allowed in descent
// -----------------------------------------------------------------------
void Network::get_plate_forces(float* plate_forces, int iter){
	int node;
	for(int i=0; i<n_tside; i++){
		node = tsideNodes[i];
		plate_forces[iter*DIM+0] += forces[node*DIM + 0];
		plate_forces[iter*DIM+1] += forces[node*DIM + 1];
		if(DIM>2){
			plate_forces[iter*DIM+2] = forces[node*DIM + 2];
		}
	}
}


// ----------------------------------------------------------------------- 
/// \brief add by Dihan
///	the function will update arrays that store the info about long chains, at each itration
///	including internal forces, end node coordinate and orientation
///	\param long_link_forces (float*) --> the array storing internal forces in each chain at this iteration
///	\param long_link_node_pos (float*) --> the array storing end-node coordinate of each chain at this iteration
///	\param long_link_orient (float*) --> the array storing orientation of each chain at this iteration
// -----------------------------------------------------------------------
void Network::get_long_link_status(float* long_link_forces, float* long_link_node_pos, float* long_link_orient, int iter){
	int node1;
	int node2;
	int num = RANDOM_LONG+RANDOM_Y;
	float s;
	float Len;
	float force;
	float orient;
	//might need to change when 3D:
	// float* pos1 = new float[DIM];
	// float* pos2 = new float[DIM];
	for (int i=num; i>0; i--){
		node1 = edges[(n_elems-i)*2+0];
		node2 = edges[(n_elems-i)*2+1];

		Len = L[n_elems-i];

		s = dist(&R[node1*DIM], &R[node2*DIM]);
		orient = orientation(&R[node1*DIM], &R[node2*DIM]);
		force = force_wlc(s, Len);

		long_link_node_pos[iter*num + (num-i)] = s/Len;
		long_link_forces[iter*num + (num-i)] = force;
		long_link_orient[iter*num + (num-i)] = orient;
	}
}



// ----------------------------------------------------------------------- 
/// \brief Dumps the network information for an iteration in a file. 
/// Filename would be <FLDR_STRING>.txt with <FLDR_STRING taken from params.h
/// If <FLDR_STRING>.txt exists, <FLDR_STRING>_i.txt would be created where i
/// is an integer decided based on copies present in current folder. 
/// 
/// \param iter (int) --> The iteration number for the simulation
// -----------------------------------------------------------------------
void Network::dump(int iter, bool first_time){
	ofstream logger;
	std::time_t result = std::time(nullptr);
	string fname = FLDR_STRING;

	if(first_time){
		int c = filename(fname);
		if(c!=-1){
			fname = fname + "_" + std::to_string(c);
		}
	}
	fname = std::string(FLDR_STRING) + "/" + fname+"_dump.txt";
	logger.open(fname, ios_base::app);

	if(first_time){
		// If opened for the first time write the params
		logger<<"Network snapshots"<<iter<<"\n";
		logger<<"File created at "<<std::asctime(std::localtime(&result));
		logger<<"PARAMS:\n";

		logger<<"Sim dimension : "<<DIM<<"\n";
		logger<<"Simulation time : "<<SIM_TIME<<"\n";
		logger<<"Simulation time-step : "<<TIME_STEP<<"\n";
		logger<<"Velocity : "<<vel_x<<"\t"<<vel_y<<"\n";
		logger<<"MAXBOUND : "<<MAXBOUND_X<<"\t"<<MAXBOUND_Y<<"\n";

		logger<<"Disorder characteristics : "<<"\n";
		logger<<" -- L_MEAN : "<<L_MEAN<<"\n";
		logger<<" -- L_STD : "<<L_STD<<"\n";
		logger<<"Others : SACBONDS = "<<SACBONDS<<"; IMPLEMENT_PBC = "<<IMPLEMENT_PBC<<"\n";
		logger<<"Rate damage : "<<RATE_DAMAGE<<"\n";		
		logger<<"Cracked? : "<<CRACKED<<"; "<<PROB_REMOVAL<<"\n";
		logger<<"-------------------------------------------------\n\n";
		
		logger<<"CONSTANTS:\n";
		logger<<"Persistence length: "<<b_poly<<"\n";
		logger<<"Temperature: "<<T<<"\n";
		logger<<"End bond a: "<<ae<<"\n";
		logger<<"End bond delx: "<<delxe<<"\n";
		logger<<"Hidden bond a: "<<af<<"\n";
		logger<<"Hidden bond delx: "<<delxf<<"\n";
		logger<<"-------------------------------------------------\n\n";

	}
	if(!first_time){   
		logger<<"INFO FOR ITER\n"<<iter<<endl;
	    // Write R
	    logger<<"START NODE POSITIONS\n";
		for(int i = 0; i < n_nodes; i++){
			logger<<i<<"\t";
			for(int j = 0; j< DIM; j++){
				logger<<R[i*DIM+j]<<"\t";
			}
			logger<<"\n";
		}
		logger<<"END NODE POSITIONS\n";

		//write edges | damage | PBC_status 
		logger<<"START ACTIVE EDGES\n";
		for(int i = 0; i < n_elems; i++){
			if(edges[2*i]!=-1 && edges[2*i+1]!=-1){
				//Note i is then the original edge index
				logger<<i<<"\t"<<edges[2*i]<<"\t"<<edges[2*i+1]<<"\t"<<damage[i]<<"\t"<<PBC[i]<<"\n";
			}
		}
		logger<<"END ACTIVE EDGES\n";

		//write forces on nodes
		logger<<"START NODE FORCES\n";
		for(int i = 0; i < n_nodes; i++){
			logger<<i<<"\t";
			for(int j = 0; j< DIM; j++){
				logger<<forces[i*DIM+j]<<"\t";
			}
			logger<<"\n";
		}
		logger<<"END NODE FORCES\n\n";
	}

    logger.close();
	cout<<"Dumped everything in "<<fname<<"!\n";
}
