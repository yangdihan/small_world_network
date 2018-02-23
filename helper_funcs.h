/**
@file helper_functions.h
\brief Defines all heler functions needed for experiments. Inline helper functions
defined here.
*/

#ifndef __helper__
#define __helper__
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <random>
#include <time.h>
#include <vector>
#include <unistd.h>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <stddef.h>
// #include "vel.h"
#include "params.h"

using namespace std;

inline int get_num_vertices(int elem_type){
	switch(elem_type){
		case 1: return 2;
		case 2: return 3;
		case 3: return 4;
		//case 4: return 4;
		//case 5: return 8;
		//case 7: return 5;
		default: return -1;
	}
}

inline bool contains(vector<int>& vec, int elem){
	return (std::find(vec.begin(), vec.end(), elem) != vec.end());
}

inline float getnorm(const float* vec, const int dim = DIM){
	float s = 0,t;
	#pragma unroll
	for (int j = 0; j<DIM; j++){
		t = vec[j];
		s += t*t;
	}
	return sqrt(s);
}

void mapping(int& edge_counter, int elem_type, int n_vertices, stringstream& input, int* edges);

void mapping(int& edge_counter, int elem_type);

void read_n(int& n_nodes, int& n_elems, string& fname);

void take_input(float* R, int* edges, int n_nodes, int n_elems, string& fname);

void __init__(float* L, float* damage, bool* PBC, int n_elems);

int filename(const string&);

inline void normalize_vector(float* result, const float* vec){
	float norm = getnorm(vec);
	#pragma unroll
	for (int i = 0; i<DIM; i++){
		result[i] = vec[i] / norm;
	}
}

inline bool does_file_exist(string& fname){
	ifstream infile(fname);
	return infile.good();
}

inline void normalize_vector(float* vec){
	float norm = getnorm(vec);
	#pragma unroll
	for (int i = 0; i<DIM; i++){
		vec[i] = vec[i] / norm;
	}
}

inline void unitvector(float* result, float* r1, float* r2){
	#pragma unroll
	for (int j = 0; j<DIM; j++){
		result[j] = r1[j] - r2[j];
	}
	normalize_vector(result);
}

inline float force_wlc(float x, float L){
	float t = x / L;
	if (t < 0.99){ return kB*T / b_poly * (t + 1.0 / 4.0 / pow((1 - t), 2) - 1.0 / 4.0); }
	else { return 1000; }
}

inline void convert_to_vector(float* result, const float mag, const float* direction){
	#pragma unroll
	for (int i = 0; i<DIM; i++){
		result[i] = mag*direction[i];
	}
}

inline float dist(const float* r1, const float* r2){
	float s = 0.0, t;
	#pragma unroll
	for (int j = 0; j<DIM; j++){
		t = r1[j] - r2[j];
		s += t * t;
	}
	return sqrt(s);
}

inline float orientation(const float* r1, const float* r2){//return the angle
	float dx = r1[0] - r2[0];
	float dy = r1[1] - r2[1];
	float tan = dy/dx;
	float result = atan (tan) * 180 / PI;
	return result;
}

void forcevector(float* result, float* r1, float* r2, float L);


inline float kfe(float force_mag){
	return ae * exp(force_mag * delxe / kB / T);
}

inline bool ismember(int item, int* array,size_t size){
	bool is_inside = false; 
	for(int k=0; k<size; k++){
		if(array[k]==item){
			is_inside = true; 
			break;
		}
	}
	return is_inside;
}

float getabsmax(float* arr, size_t sizeofarr);

template <typename t>
void write_to_file(string& fname, t* arr, int rows, int cols){

	ofstream logger;
	std::time_t result = std::time(nullptr);
	logger.open(fname, ios::trunc|ios_base::out);


	logger<<"1D pulling of a 2D gel"<<"\n";
	logger<<"File created at "<<std::asctime(std::localtime(&result));
	logger<<"Rate damage : "<<RATE_DAMAGE<<"\n";

	logger<<"Sim dimension : "<<DIM<<"\n";
	logger<<"Simulation time : "<<SIM_TIME<<"\n";
	logger<<"Simulation time-step : "<<TIME_STEP<<"\n";
	logger<<"Velocity : "<<vel_x<<"\t"<<vel_y<<"\n";
	logger<<"MAXBOUND : "<<MAXBOUND_X<<"\t"<<MAXBOUND_Y<<"\n";

	logger<<"Disorder characteristics : "<<"\n";
	logger<<" -- L_MEAN : "<<L_MEAN<<"\n";
	logger<<" -- L_STD : "<<L_STD<<"\n";
	logger<<"Others : SACBONDS = "<<SACBONDS<<"; IMPLEMENT_PBC = "<<IMPLEMENT_PBC<<"\n";
	
	logger<<"Cracked? : "<<CRACKED<<"; "<<PROB_REMOVAL<<"\n";

    // logger.open(fname, ios::trunc|ios_base::out);
	for(int i =0; i < rows; i++){
		for(int j = 0; j< cols; j++){
			logger<<arr[i*cols + j]<<"\t";
		}
		logger<<"\n";
	}
    logger.close();
	cout<<"Stored everything in "<<fname<<"!\n";
}

template <typename t>
void write_edge_number(string& fname2, t* remain_chains, int row){
	ofstream logger;

	logger.open(fname2, ios::trunc|ios_base::out);
	for(int j = 0; j< row; j++){
		// cout << "#####" << remain_chains[j] << endl;
		logger<<remain_chains[j]<<"\n";
		// logger<<"\n";
	}
    logger.close();
	cout<<"finish store edge numbers "<<fname2<<"!\n";
}

template <typename t>
void write_long_link(string& fname3, t* long_link_forces, t* long_link_node_pos, t* long_link_orient, int row){
	ofstream logger;
	logger.open(fname3, ios::trunc|ios_base::out);
	int num = RANDOM_LONG+RANDOM_Y;
	for(int j = 0; j< row; j++){
		for (int i=0; i<num; i++){
			logger<< long_link_forces[j*num + i] <<"\t";
			logger<< long_link_node_pos[j*num + i] << "\t";
			logger<< long_link_orient[j*num + i] << "\t";
		}
		logger << "\n";
	}
	logger.close();
	cout<<"finish store long link info "<<fname3<<"!\n";
}

// void animation(){
// 	Gnuplot gnu;
// 	gnu.cmd(ENCODER = system('which mencoder'));
// 	if (strlen(ENCODER)==0) print '=== mencoder not found ===';
// 	gnu.cmd(std::string(FLDR_STRING)+'/*.png -mf fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o animation.avi');
// 	// system(CMD)
// }

inline float kf(float force){
	return af*expf(force*delxf/kB/T);
}

void __init__(float* L, int* m, float* damage, float* sacdamage, bool* PBC, int n_elems);

#endif
