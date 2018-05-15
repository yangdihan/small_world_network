#ifndef __params__
#define __params__
// ssh dyang42@cc-login.campuscluster.illinois.edu
// PATH=/home/dyang42/scratch/MPICH/mpich-install/bin/:$PATH
// export PATH
// module load intel/17.0 
// module load gcc/6.2.0 
// module load mvapich2/2.2-intel-17.0
// qsub -l walltime=03:00:00,nodes=2:ppn=2 serial.pbs 
// qsub -l walltime=00:30:00,nodes=2:ppn=2 serial.pbs 
// qstat -u dyang42
// FREQUENTLY USED
// Please refer to original Documentation file
// 1. PBC for long chains, give up mirror case
// 2. figure out nonlinearity for both overall network and single chain
// 3. sigma = E*eps
// 4. RATE_DAMAGE = true
#define MAXBOUND_X 500.0f
#define MAXBOUND_Y 800.0f
#define MESH_CRACK false
#define TIME_STEP 1e-2
#define SIM_TIME 40
#define L_MEAN 50.0f
#define L_STD 5.0f
#define FLDR_STRING "playground"

// add by Dihan:
// I moved the defination of parameters vel_x & vel_y from the separate vel.h file here.
#define vel_x 0.0f
#define vel_y 50.0f

// If the total weight is less than this value, the weight of each chain will be multiplied by certain factor to reach this threshold. Because when the total weight is too small, many flexure properties will not be stable and conclusive 
#define WEIGHT_GOAL 1.0e6

// bool that makes left and right boundary of the network constrained in lateral displacement
#define ROLLER false

// bool that makes nodes too close to top/bottom the ghost nodes.
#define GHOST true

// int that indicates the number of random long chains imposed
#define RANDOM_LONG 50

// int that indicates the number of random long chains imposed in vertical directions only
#define RANDOM_Y 0

// float that indicates the "prestretch rate". I defined "prestretch rate" as a parameter to describe the extent of prestretch in addtional chains that 0 means the prestetch is same as the original short chains; 1 means the contour length is same as the original short chains.
#define PRESTRETCH 0.5

// int that indicates the frequency of output a EPS format snapshot of the network configuration (every how many iterations)
#define EPS 0

// int that indicates the frequency of output a PNG format snapshot of the network configuration (every how many iterations)
#define PNG 2

// string that indicates the pattern style is "layer" or "spot". "layer" means the network is straitified by sparse and dense regions; "spot" means the there are round-shape sparse regions that distribute over the network.
#define PATTERN_TYPE "none"

// int that indicates the number of patterned regions: how many stripes for "layer" or how many round-regions for "spot".
#define PATTERN_NUM 3

// float that indicates how many times sparser the patterned region will be
#define PATTERN_RATE 2

// PBC
#define IMPLEMENT_PBC true

// Breaking Criteria
#define RATE_DAMAGE true

// other constant definations:
#define PI 3.141592653
#define DIM 2
#define TOL 1e-6
#define STEPS int(SIM_TIME/TIME_STEP)
#define PAD_X MAXBOUND_X*1.03
#define PAD_Y MAXBOUND_Y*1.03
#define SACBONDS false
#define CRACKED false

#if CRACKED
#define PROB_REMOVAL 1.0
#else
#define PROB_REMOVAL 0.0
#endif

#define __constants__

#define kB 1.38064852e-5
#define b_poly 0.1
#define T 300.0
#define ae 0.1
#define delxe 0.15
#define af 0.1
#define delxf 0.25

#endif
