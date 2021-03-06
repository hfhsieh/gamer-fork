#ifndef __CUFLU_H__
#define __CUFLU_H__



// *********************************************************************
// ** This header will be included by all CPU/GPU fluid/ELBDM solvers **
// *********************************************************************


// include "Macro.h" and "Typedef.h" here since the header "GAMER.h" is NOT included in GPU solvers
#ifdef __CUDACC__
# include "Macro.h"
# include "Typedef.h"
#else
# include "GAMER.h"
#endif


// allow GPU to output messages in the debug mode
#ifdef GAMER_DEBUG
#  include "stdio.h"
#endif


// faster integer multiplication in Fermi
#if ( defined __CUDACC__  &&  __CUDA_ARCH__ >= 200 )
#  define __umul24( a, b )   ( (a)*(b) )
#  define  __mul24( a, b )   ( (a)*(b) )
#endif


// bitwise reproducibility in flux and electric field fix-up operations
#ifdef BITWISE_REPRODUCIBILITY
#  define BIT_REP_FLUX
#endif

// enable BIT_REP_ELECTRIC by default even when BITWISE_REPRODUCIBILITY is off
// --> it ensures that the B field on the common interface between two nearby patches are fully
//     consistent with each other (even the round-off errors are the same)
//     --> reducing the div(B) errors significantly
#ifdef MHD
//#ifdef BITWISE_REPRODUCIBILITY
#  define BIT_REP_ELECTRIC
//#endif
#endif



// #################################
// ## macros for different models ##
// #################################

// 1. hydro macro
//=========================================================================================
#if   ( MODEL == HYDRO )

// size of different arrays
// ** to reduce the GPU memory consumption, large arrays in the fluid solvers are reused as much as possible
// ** --> the strides of arrays can change when accessed by different routines for different purposes

// N_SLOPE_PPM          : size of Slope_PPM[]
// N_FC_VAR             : size of FC_Var[]
// N_FC_FLUX            : size of FC_Flux[]
// N_FL_FLUX/N_HF_FLUX  : for accessing FC_Flux[]
//                        --> may be different from N_FC_FLUX
//                            --> for example, in MHM_RP FC_Flux[] is also linked to Half_Flux[] used by
//                                Hydro_RiemannPredict_Flux() and Hydro_RiemannPredict()
//                            --> for the latter two routines, Half_Flux[] is accessed with N_HF_FLUX
//                                that is smaller than N_FC_FLUX
// N_HF_VAR             : for accessing PriVar_Half[], which is linked to PriVar[] with the size FLU_NXT^3
//                        --> also for accessing FC_B_Half[] in MHD

// NWAVE                : number of characteristic waves
// NCOMP_TOTAL_PLUS_MAG : total number of fluid variables plus magnetic field
// MAG_OFFSET           : array offset of magnetic field for arrays with the size NCOMP_TOTAL_PLUS_MAG
//
#define NCOMP_TOTAL_PLUS_MAG     ( NCOMP_TOTAL + NCOMP_MAG )

#ifdef MHD
#  define NWAVE                  ( NCOMP_FLUID + 2 )
#  define MAG_OFFSET             ( NCOMP_TOTAL )
#else
#  define NWAVE                  ( NCOMP_FLUID )
#  define MAG_OFFSET             ( NULL_INT )
#endif

#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

#  if   ( FLU_SCHEME == MHM )

#     define N_FC_VAR            ( PS2 + 2 )
#     define N_FL_FLUX           ( PS2 + 1 )
#     define N_FC_FLUX           ( N_FL_FLUX )

#  elif ( FLU_SCHEME == MHM_RP )

#    ifdef MHD
#     define N_FC_VAR            ( PS2 + 2 )
#     define N_FL_FLUX           ( PS2 + 2 )
#     define N_HF_FLUX           ( FLU_NXT )
#    else
#     define N_FC_VAR            ( PS2 + 2 )
#     define N_FL_FLUX           ( PS2 + 1 )
#     define N_HF_FLUX           ( FLU_NXT - 1 )
#    endif
#     define N_FC_FLUX           ( N_HF_FLUX )
#     define N_HF_VAR            ( FLU_NXT - 2 )

#  elif ( FLU_SCHEME == CTU )

#    ifdef MHD
#     define N_FC_VAR            ( PS2 + 4 )
#     define N_FL_FLUX           ( N_FC_VAR - 2 )
#     define N_HF_VAR            ( PS2 + 2 )
#    else
#     define N_FC_VAR            ( PS2 + 2 )
#     define N_FL_FLUX           ( N_FC_VAR )
#    endif
#     define N_HF_FLUX           ( N_FC_VAR )
#     define N_FC_FLUX           ( N_HF_FLUX )

#  endif // FLU_SCHEME

#  define N_SLOPE_PPM            ( N_FC_VAR + 2 )

#  ifdef MHD
#   define N_HF_ELE              ( N_FC_FLUX - 1 )
#   define N_FL_ELE              ( N_FL_FLUX - 1 )
#   define N_EC_ELE              ( N_FC_FLUX - 1 )
#  else
#   define N_EC_ELE              0
#  endif

#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )


// check non-physical negative values (e.g., negative density) for the fluid solver
#if ( defined GAMER_DEBUG  &&  MODEL == HYDRO )
#  define CHECK_NEGATIVE_IN_FLUID
#endif

#ifdef CHECK_NEGATIVE_IN_FLUID
#  include "stdio.h"
   bool Hydro_CheckNegative( const real Input );
#endif


// perform spatial data reconstruction in characteristic variables (default: primitive variables)
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
#  ifndef GRAVITY
#     define CHAR_RECONSTRUCTION
#  endif
#endif


// verify that the density and pressure in the intermediate states of Roe's Riemann solver are positive.
// --> if either is negative, we switch to other Riemann solvers (EXACT/HLLE/HLLC/HLLD)
#if (  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  &&  RSOLVER == ROE  )
#  ifdef MHD
//#     define CHECK_INTERMEDIATE    HLLD
#     define CHECK_INTERMEDIATE    HLLE
#  else
//#     define CHECK_INTERMEDIATE    HLLC
#     define CHECK_INTERMEDIATE    HLLE
#  endif
#endif


// use Eulerian with Y factor for Roe Solver in MHD
#if (  defined MHD  &&  ( RSOLVER == ROE || RSOLVER == HLLE )  )
#  define EULERY
#endif


// do not use the reference states for HLL solvers during the data reconstruction, as suggested in ATHENA
#if (  defined RSOLVER  &&  ( RSOLVER == HLLE || RSOLVER == HLLC || RSOLVER == HLLD )  )

#  define HLL_NO_REF_STATE

// include waves both from left and right directions during the data reconstruction, as suggested in ATHENA
#  ifdef HLL_NO_REF_STATE
#     define HLL_INCLUDE_ALL_WAVES
#  endif

#endif


// 2. ELBDM macro
//=========================================================================================
#elif ( MODEL == ELBDM )


#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL



// ###############################################################################################
// number of threads in x/y directions for different solvers:
// Better to be a multiple of myWarpAllocationGranularity (WAG) * limitThreadsPerWarp(TPW)
// --> Please refer to CUDA_Occupancy_Calculator
//
// ComputeCompability  WAG   TPW   WAG*TPW
//       2.x             2    32        64
//       3.x             4    32       128
//       5.x             4    32       128
// ###############################################################################################

// 1. hydro solver
//=========================================================================================
#if ( MODEL == HYDRO )
#if   ( FLU_SCHEME == RTVD )

#     define FLU_BLOCK_SIZE_X       FLU_NXT

#  ifdef FLOAT8
#     define FLU_BLOCK_SIZE_Y       4
#  else
#     define FLU_BLOCK_SIZE_Y       8
#  endif

#elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )

#  if   ( GPU_ARCH == FERMI )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512
#     endif
#  elif ( GPU_ARCH == KEPLER )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512
#     endif
#  elif ( GPU_ARCH == MAXWELL )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == PASCAL )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == VOLTA )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == TURING )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  else
#     define FLU_BLOCK_SIZE_X       NULL_INT
#     ifdef GPU
#     error : UNKNOWN GPU_ARCH !!
#     endif
#  endif

#     define FLU_BLOCK_SIZE_Y       1

#elif ( FLU_SCHEME == CTU )

#  if   ( GPU_ARCH == FERMI )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512
#     endif
#  elif ( GPU_ARCH == KEPLER )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512
#     endif
#  elif ( GPU_ARCH == MAXWELL )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == PASCAL )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == VOLTA )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == TURING )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  else
#     define FLU_BLOCK_SIZE_X       NULL_INT
#     ifdef GPU
#     error : UNKNOWN GPU_ARCH !!
#     endif
#  endif

#     define FLU_BLOCK_SIZE_Y       1

#else
#  error : ERROR : unsupported hydro scheme in the makefile !!
#endif


// 2. ELBDM kinematic solver
//=========================================================================================
#elif ( MODEL == ELBDM )
#     define FLU_BLOCK_SIZE_X       PS2

#  if   ( GPU_ARCH == FERMI )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    4
#     else
#        define FLU_BLOCK_SIZE_Y    8
#     endif

#  elif ( GPU_ARCH == KEPLER )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    16    // not optimized yet
#     else
#        define FLU_BLOCK_SIZE_Y    32    // not optimized yet
#     endif

#  elif ( GPU_ARCH == MAXWELL )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    16    // not optimized yet
#     else
#        define FLU_BLOCK_SIZE_Y    32    // not optimized yet
#     endif

#  elif ( GPU_ARCH == PASCAL )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    16    // not optimized yet
#     else
#        define FLU_BLOCK_SIZE_Y    32    // not optimized yet
#     endif

#  elif ( GPU_ARCH == VOLTA )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    16    // not optimized yet
#     else
#        define FLU_BLOCK_SIZE_Y    32    // not optimized yet
#     endif

#  elif ( GPU_ARCH == TURING )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    16    // not optimized yet
#     else
#        define FLU_BLOCK_SIZE_Y    32    // not optimized yet
#     endif

#  else
#        define FLU_BLOCK_SIZE_Y    NULL_INT
#        ifdef GPU
#        error : UNKNOWN GPU_ARCH !!
#        endif
#  endif

#endif // MODEL


// 3. dt solver for fluid
//=========================================================================================
#     define DT_FLU_BLOCK_SIZE      512

// use shuffle reduction in the KEPLER and later GPUs
#  if ( GPU_ARCH == KEPLER  ||  GPU_ARCH == MAXWELL  ||  GPU_ARCH == PASCAL  ||  GPU_ARCH == VOLTA  ||  GPU_ARCH == TURING )
#     define DT_FLU_USE_SHUFFLE
#  endif



// warp size (which must be the same as the CUDA predefined constant "warpSize")
// --> please refer to https://en.wikipedia.org/wiki/CUDA#Version_features_and_specifications
//     for information on warp size
#ifdef __CUDACC__
#if ( GPU_ARCH == FERMI  ||  GPU_ARCH == KEPLER  ||  GPU_ARCH == MAXWELL  ||  GPU_ARCH == PASCAL  ||  GPU_ARCH == VOLTA  ||  GPU_ARCH == TURING )
// CUPOT.h will define WARP_SIZE as well
#  ifndef WARP_SIZE
#  define WARP_SIZE 32
#  endif
#elif defined GPU
#  error : UNKNOWN GPU_ARCH !!
#endif
#endif // #ifdef __CUDACC__



// #########################
// ## CPU/GPU integration ##
// #########################

// GPU device function specifier
#ifdef __CUDACC__
# define GPU_DEVICE __forceinline__ __device__
#else
# define GPU_DEVICE
#endif

// unified CPU/GPU loop
#ifdef __CUDACC__
# define CGPU_LOOP( var, niter )    for (int (var)=threadIdx.x; (var)<(niter); (var)+=blockDim.x)
#else
# define CGPU_LOOP( var, niter )    for (int (var)=0;           (var)<(niter); (var)++          )
#endif



#endif // #ifndef __CUFLU_H__
