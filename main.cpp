/** @mainpage
@authors {Konik Kothari kkothar3@illinois.edu, Sahil Gupta sjgupta2@illinois.edu}
@version Revision 0.0
@brief This codebase simulates discrete polymer network as
graphs with and without cross-linker bond dynamics. Please 
feel free to add or modify the code as per your usage with the
only condition that this codebase is referenced.
The authors have tried to provide modules for users to build 
their own custom networks and design experiments with them.
The authors do not claim any responsibility for correctness or
accuracy of this code and welcome the community to critique and
improve the code as they wish. This code was written under advisors 
Yuhang Hu and Ahmed Elbanna at University of Illinois at Urbana-Champaign.
@date Thursday April 13, 2017
*/

/**
@file main.cpp
\brief This file includes the code to check whether an MPI implementation 
has been asked by the user and run the experiment accordingly.
*/

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
#include <mpi.h>
#include <assert.h>

#include "mpi_network.h"
#include "network.h"
using namespace std;


int main(int argc, char* argv[]) {
	/**< 
	\brief  Takes in 2 command line arguments, <br>
	[1] filename (string) of the .msh file
	generated by GMSH, this file creates the network <br>
	[2] MPI_flag (int) if want to run the code on MPI. 
	Note the syntax for MPI runs.

	Note: mpirun -np <numprocs> <executable> <filename.msh> 1 
	*/
	string path = argv[1];

/** add by Dihan
 *	create main folder to store everything inside
 *	naming: meshSize_crack_longLinkNumber_prestretchValue
 *	create two folders: 
 * 	"frames" to store the .png format pictures as frames for final animation;
 *	"graphs" to store the .eps files for publication use.
 */
	string folder_name = "/"+to_string(int(MAXBOUND_X))+"_"+to_string(int(MAXBOUND_Y));
	if (MESH_CRACK){
		folder_name += "_crack";
	}
	if (RANDOM_LONG > 0){
		folder_name += "_"+to_string(int(RANDOM_LONG))+"_"+to_string(PRESTRETCH);
	}
	if (RANDOM_Y > 0){
		folder_name += "_"+to_string(int(RANDOM_Y))+"_v"+to_string(PRESTRETCH);
	}
	if (std::string(PATTERN_TYPE) != "none"){
		folder_name += "_"+std::string(PATTERN_TYPE)+"_num"+to_string(PATTERN_NUM)+"_rate"+to_string(PATTERN_RATE);
	}
	folder_name += "_vel_"+to_string(int(vel_y));
	string arg0 = "mkdir ./"+std::string(FLDR_STRING)+folder_name;
	system(arg0.c_str());

	if (EPS != 0){
		string arg1 = "mkdir ./"+std::string(FLDR_STRING)+folder_name+"/graphs";
		system(arg1.c_str());
	}
	if (PNG != 0){
		string arg2 = "mkdir ./"+std::string(FLDR_STRING)+folder_name+"/frames";
		system(arg2.c_str());
	}
	

	#if SACBONDS
	#define DECL_NET sacNetwork test_network(path)
	#else
	#define DECL_NET Network test_network(path)
	#endif
	DECL_NET;

/** add by Dihan
 *	call the patterning function to add pattern according to pre-defined patterning condition from param.h
 */
	test_network.patterning(PATTERN_TYPE, PATTERN_NUM, PATTERN_RATE);
	if(CRACKED){
		// Specific crack
		Crack a;
		a.setter(MAXBOUND_X, MAXBOUND_Y/2.0, MAXBOUND_X/10.0, MAXBOUND_Y/10.0, 0.0, 1.0);
		Cracklist definite_cracks(a);
		test_network.apply_crack(definite_cracks);
		cout<<__LINE__<<endl;
	}



	if (*argv[2]!='1') {
		cout<<"Running serial version of the code: \n";
	
		float weight_multiplier;
		float weight = test_network.get_weight();
		if (weight<WEIGHT_GOAL){
			weight_multiplier = test_network.set_weight(WEIGHT_GOAL);
		}
		bool should_stop = test_network.get_stats();
		
	/** add by Dihan 
	 *	call functions to add random long links to the network when RANDOM_LONG or RANDOM_Y is not zero
	 *	store info for long links including internal forces, node position and orientation in three array
	 *	arrays are initialized here
	 */
		bool long_links = (RANDOM_LONG+RANDOM_Y > 0);
		// float* long_link_forces;
		// float* long_link_node_pos;
		// float* long_link_orient;
		int curr_n_edges = test_network.get_current_edges();

		if (long_links){
			test_network.add_long_range_egdes_random(RANDOM_LONG, folder_name);
			test_network.add_long_range_egdes_y(RANDOM_Y, folder_name);
		// 	long_link_forces = (float*)malloc(sizeof(float)*(RANDOM_LONG+RANDOM_Y)*STEPS);
		// 	memset(long_link_forces, 0.0, sizeof(float)*(RANDOM_LONG+RANDOM_Y)*STEPS);
		// 	long_link_node_pos = (float*)malloc(sizeof(float)*(RANDOM_LONG+RANDOM_Y)*STEPS);
		// 	memset(long_link_node_pos, 0.0, sizeof(float)*(RANDOM_LONG+RANDOM_Y)*STEPS);
		// 	long_link_orient = (float*)malloc(sizeof(float)*(RANDOM_LONG+RANDOM_Y)*STEPS);
		// 	memset(long_link_orient, 0.0, sizeof(float)*(RANDOM_LONG+RANDOM_Y)*STEPS);
		}

		// int old_n_edges = test_network.get_current_edges();
		// int curr_n_edges = old_n_edges;
		// int total_n_edges = curr_n_edges;

		if(should_stop){cout<<"Simulation needs to stop!\n";return 0;}

		float* plate_forces;
		plate_forces = (float*)malloc(sizeof(float)*DIM*STEPS);
		memset(plate_forces, 0.0, STEPS*DIM*sizeof(float));

	/** add by Dihan
	 *	store number of remaining chains in a array of int
	 *	initialize array here
	 */
		// int* remain_chains;
		// remain_chains = (int*)malloc(sizeof(int)*STEPS);
		// memset(remain_chains, 0, STEPS*sizeof(int));

	/** add by Dihan
	 *	determine whether need to output .png and .eps files at certain frequency
	 */
		test_network.plotNetwork(0, true, folder_name);
		test_network.plotFrames(0, true, folder_name);

		// test_network.dump(0, true);

		// For time is endless but your time...not so much!  
		clock_t t = clock();
		float old_time_per_iter;
		float new_time_per_iter = 100.0f;
		cout<<"\n Will run for "<<STEPS<<":\n";
		

	/** add by Dihan
	 *	initialize float about top plate force, and bool to control whether to move the top plate at given rate
	 *	The bool is set to false at first: we will let the network to equilibrium itself before start moving it
	 *	The strategy is if at next iteration, the top plate force is no longer changing, then we start moving top plate since next iteration
	 */
		bool should_move = false;
		int first_few = STEPS/100;
		double plate_force_new;
		double plate_force_min;

		for(int i = 1; i<STEPS; i++ ){
			// The "optimize-move_top_plate-get_plate_forces" do all the work
			test_network.optimize();
			test_network.get_plate_forces(plate_forces, i);

		/** add by Dihan
		 *	if force start to increase, we should move top plate
		 */
			plate_force_new = -plate_forces[DIM*i + 1];
			if (i == 1 || plate_force_new < plate_force_min){
				plate_force_min = plate_force_new;
			}else{
				should_move = true;
			}
			if (should_move){
				test_network.move_top_plate();
			}

		/** add by Dihan
		 *	update info about long_link_forces, long_link_node_pos, long_link_orient at this iteration in corresponding array
		 */
			// if (long_links){
			// 	test_network.get_long_link_status(long_link_forces, long_link_node_pos, long_link_orient, i);
			// }

			should_stop = test_network.get_stats();

		/** add by Dihan
		 *	update number of remaining chains at this iteration in corresponding array
		 */
			curr_n_edges = test_network.get_current_edges();

			// test_network.get_edge_number(remain_chains, i, curr_n_edges);

		/** add by Dihan
		 *	output png or eps file if should at this iteration
		 */
			if (PNG != 0){
				if (i%int(PNG) == 0){
					test_network.plotFrames(i, false, folder_name);
				}
			}
			if (EPS != 0){
				if (i%int(EPS) == 0){
					test_network.plotNetwork(i, false, folder_name);
				}
			}
			
			// if(curr_n_edges<=old_n_edges){
			// 	test_network.plotNetwork(i, false);
			// 	test_network.dump(i);
			// }
			// test_network.dump(i);
			
			if(should_stop){
				break;
			}
			new_time_per_iter = float(clock()-t)/CLOCKS_PER_SEC;
			if(i==0){
				old_time_per_iter = new_time_per_iter;
			}
			// if(new_time_per_iter < 0.1*old_time_per_iter){
			// 	cout<<"Seems like very few edges remain! \n";
			// 	break;
			// }
			cout<<"Step "<<(i+1)<<" took "<<float(clock()-t)/CLOCKS_PER_SEC<<" s\n";
			t = clock();

		}

	/** add by Dihan
	 *	define name of output files and generate corresponding txt files
	 *	forces.txt stores the vertical and horizontal forces on top plate
	 *	remain_chains.txt stores the number of remaining chains at each iteration
	 *	add_long_link_info.txt stores the info about additional long chains at each iteration
	 */
		string fname = std::string(FLDR_STRING)+folder_name + "/" + "forces.txt";
		// string fname2 = std::string(FLDR_STRING) + "/" + "remain_chains.txt";
		// string fname3 = std::string(FLDR_STRING) + "/" + "add_long_link_info.txt";
		write_to_file<float>(fname, plate_forces, STEPS, DIM);
		// write_edge_number<int>(fname2, remain_chains, STEPS);
		// if (long_links){
		// 	// write_long_link<float>(fname3, long_link_forces, long_link_node_pos, long_link_orient, STEPS);
		// 	free(long_link_forces);
		// 	free(long_link_node_pos);
		// 	free(long_link_orient);
		// }

	/** add by Dihan
	 *	generate .mp4 animation from .png files in frames folder
	 *	copy plot.py file to same folder which is used for ploting force-stretch diagram
	 */
		if (PNG != 0){
			cout << "#### rendering animation ####" << endl;
			string arg3 = "ffmpeg -f image2 -r "+std::to_string(int(20/PNG))+" -i ./"+std::string(FLDR_STRING)+folder_name+"/frames/%05d.png -vcodec mpeg4 -y ./"+std::string(FLDR_STRING)+folder_name+"/movie.mp4";
			system(arg3.c_str());
			string arg4 = "rm ./"+std::string(FLDR_STRING)+folder_name+"/frames/*.png";
			system(arg4.c_str());
			cout << "#### animation saved as movie.mp4 ####" << endl;
		}
		string arg5 = "cp plot.py ./"+std::string(FLDR_STRING)+folder_name+"/";
		system(arg5.c_str());


		free(plate_forces);
		// free(remain_chains);


	}
	else {

		// Settle in! This code is murky af
		cout<<"Running MPI version of the code: "<<endl;
		MPI_Init(NULL, NULL);
		
		// Get the number of processes
	  	int world_size;
	  	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	  	// Get the rank of the process
	  	int world_rank;
	  	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	  	MPI_Network * main_network = new MPI_Network(test_network);
		
	  	float* plate_forces = NULL;

  
		
		bool should_stop = main_network->get_stats();
		if(main_network->get_weight()==0){
			cout<<"Problem with MPI constructor! Reading 0 weight! Exiting"<<endl;
			MPI_Finalize();
			return 0;
		}

		if(should_stop || world_size % 2 == 1) {
			//Always take even number of processors
			cout<<"Got stop signal! stop flag: \t"<<should_stop<<\
				"world_size: "<<world_size<<endl;
			MPI_Finalize();
			return 0;
		}
	    
	    // world_rank 0 handles forces and R sync
		if (world_rank == 0) {
			plate_forces = (float*)malloc(sizeof(float)*DIM*STEPS/NSYNC);
			memset(plate_forces, 0, STEPS/NSYNC*DIM*sizeof(float));
		}
		
		// cout<<"world rank:  "<<world_rank<<main_network->n_elems<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
		main_network->init_MPI(world_rank, world_size);

		cout<<__LINE__<<endl;
		MPI_Bcast(main_network->L, main_network->n_elems, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(main_network->PBC, main_network->n_elems, MPI_C_BOOL, 0, MPI_COMM_WORLD); 

		size_t r_size = main_network->n_nodes * DIM * world_size;
		float * R_buffer; 
		R_buffer = (float*)malloc(r_size*sizeof(float));//buffer to gather the R from all nodes
		
		float * forces_buffer; 
		forces_buffer = (float*)malloc(r_size*sizeof(float));//buffer to gather the forces from all nodes
		// Force buffer needed to calculate plate_forces
		
		int * chunk_nodes_buffer = new int[main_network->chunk_nodes_len*world_size];
		cout<<"world rank: "<<world_rank<< " chunk len = "<<main_network->chunk_nodes_len<<endl;
				
		MPI_Gather(main_network->chunk_nodes, main_network->chunk_nodes_len, MPI_INT, chunk_nodes_buffer, main_network->chunk_nodes_len, MPI_INT, 0, MPI_COMM_WORLD);

		// Uniqueness of partition check
		if (world_rank == 0) {
			int chunk_sum = 0;
			for (int i = 0; i < main_network->chunk_nodes_len * world_size; i++) {
				if (chunk_nodes_buffer[i] != -1) {
					chunk_sum += chunk_nodes_buffer[i];
				}
			}
			
			int nn = main_network->n_nodes;
			int id;
			for(int d = 0 ; d<main_network->chunk_edges_len; d++){
				id = main_network->chunk_edges[d];
				if(id!=-1){
					if(main_network->edges[2*id] >= nn || main_network->edges[2*id + 1] >= nn){
								cout<<"Node is "<<main_network->edges[2*id]<<" for index "<<id<<endl;
					}	
				}
			}
			if (chunk_sum != (nn*nn - nn)/2 ) {
				cout << chunk_sum << " | "<< (nn*nn - nn)/2<<endl; 
				cout << "Uneven chunk partitioning" << endl;

			}
		}
		
		cout << "World rank proc "<<world_rank << " starting the loop:" << endl;

		int iter = 0; // needed to write forces later

		clock_t t = clock(); 
		for(iter = 0; iter<STEPS; iter++){
			if((iter+1)%100 == 0){
				cout<<"That took "<<(clock()-t)/CLOCKS_PER_SEC<<" s\n";
				t = clock();  // reset clock
				if(world_rank==0){
					cout<<iter+1<<endl; 
					main_network->plotNetwork(iter, false, folder_name);
					main_network->get_stats();
				}
			}
			main_network->optimize();
			MPI_Barrier(MPI_COMM_WORLD);


			MPI_Gather(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, R_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);

			if((iter+1)%NSYNC == 0){
				MPI_Gather(main_network->forces, main_network->n_nodes * DIM, MPI_FLOAT, forces_buffer, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
			}

			// syncing R and forces
			if (world_rank == 0) {
				int node_to_sync  = 0;
				for (int i = 0; i < world_size; i += 1) {
					for (int j = i*main_network->chunk_nodes_len; j < (i+1)*main_network->chunk_nodes_len; j++) {
						node_to_sync = chunk_nodes_buffer[j];
						if (node_to_sync == -1) {
							break;
						}
						else{
							main_network->R[DIM * node_to_sync] = R_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync];
							main_network->R[DIM * node_to_sync + 1] = R_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync + 1];
							if((iter+1)%NSYNC == 0){
								main_network->forces[DIM * node_to_sync] = forces_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync];
								main_network->forces[DIM * node_to_sync + 1] = forces_buffer[main_network->n_nodes * DIM * i + DIM * node_to_sync + 1];
							}
						}
					}
					
				}
				if((iter+1)%NSYNC == 0){
					cout << "Synced forces" << endl;
					main_network->get_plate_forces(plate_forces, iter/NSYNC);
				}
				main_network->move_top_plate();
			}

			MPI_Bcast(main_network->R, main_network->n_nodes * DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);
			
			//Barrier required as Bcast is not synchronous and does not block other processes from continuing
			MPI_Barrier(MPI_COMM_WORLD);

		} // the simulation loop ends here


		if (world_rank == 0) {
			string sb = SACBONDS ? "true" : "false" ; 
			string fname = FLDR_STRING + std::to_string(L_STD/L_MEAN) + "_" + sb + ".txt";
			write_to_file<float>(fname, plate_forces, STEPS, DIM);
			free(plate_forces);
			plate_forces = NULL;
		}
		
		//Needed to not have double free, corruption error
		delete[] R_buffer;
		R_buffer = NULL;
		delete[] forces_buffer;
		forces_buffer = NULL;
		delete[] chunk_nodes_buffer;
		chunk_nodes_buffer = NULL;
		delete main_network;
		main_network = NULL;
		cout << "Made it to the end! Exiting now." << endl;
		MPI_Finalize();

	}

	return 0;
}
