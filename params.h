#ifndef __params__
#define __params__

// FREQUENTLY USED
#define MAXBOUND_X 50.0f
#define MAXBOUND_Y 50.0f
#define TIME_STEP 1e-2
#define SIM_TIME 20
#define L_MEAN 150.0f
#define L_STD 25.0f
#define FLDR_STRING "test_fix_side"
//add by Dihan
#define SIDE_BC false
#define ROLLER true
#define RANDOM_LONG 3
#define RANDOM_Y 0

#define DIM 2
#define TOL 1e-6
#define STEPS int(SIM_TIME/TIME_STEP)
#define PAD_X MAXBOUND_X*1.03
#define PAD_Y MAXBOUND_Y*1.03
#define SACBONDS false
#define IMPLEMENT_PBC false
#define CRACKED false
#define RATE_DAMAGE true


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
