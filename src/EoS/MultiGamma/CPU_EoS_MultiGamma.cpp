#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#include "CUFLU_Shared_FluUtility.cu"
#endif

#if ( MODEL == HYDRO )



#ifdef __CUDACC__

__device__ static real EoS_DensEint2Pres_MultiGamma( const real Dens, const real Eint, const real Passive[],
                                          const double AuxArray_Flt[], const int AuxArray_Int[],
                                          const real *const Table[EOS_NTABLE_MAX] );
__device__ static real EoS_DensEint2Temp_MultiGamma( const real Dens, const real Eint, const real Passive[],
                                          const double AuxArray_Flt[], const int AuxArray_Int[],
                                          const real *const Table[EOS_NTABLE_MAX] );
__device__ static void EoS_General_MultiGamma( const int Mode, real Out[], const real In_Flt[], const int In_Int[],
                                               const double AuxArray_Flt[], const int AuxArray_Int[],
                                               const real *const Table[EOS_NTABLE_MAX] );

#else

static real EoS_DensEint2Pres_MultiGamma( const real Dens, const real Eint, const real Passive[],
                                          const double AuxArray_Flt[], const int AuxArray_Int[],
                                          const real *const Table[EOS_NTABLE_MAX] );
static real EoS_DensEint2Temp_MultiGamma( const real Dens, const real Eint, const real Passive[],
                                          const double AuxArray_Flt[], const int AuxArray_Int[],
                                          const real *const Table[EOS_NTABLE_MAX] );
static void EoS_General_MultiGamma( const int Mode, real Out[], const real In_Flt[], const int In_Int[],
                                    const double AuxArray_Flt[], const int AuxArray_Int[],
                                    const real *const Table[EOS_NTABLE_MAX] );

#endif



/********************************************************
1. Ideal gas EoS with a constant adiabatic index (EOS_MULTIGAMMA)

2. This file is shared by both CPU and GPU

   GPU_EoS_MultiGamma.cu -> CPU_EoS_MultiGamma.cpp

3. Three steps are required to implement an EoS

   I.   Set EoS auxiliary arrays
   II.  Implement EoS conversion functions
   III. Set EoS initialization functions
********************************************************/



// =============================================
// I. Set EoS auxiliary arrays
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_SetAuxArray_MultiGamma
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by EoS_Init_MultiGamma()
//                2. AuxArray_Flt/Int[] have the size of EOS_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//                4. Physical constants such as Const_amu/Const_kB should be set to unity when disabling OPT__UNIT
//                5. Do not change the order of AuxArray_Flt[]
//                6. Ref: Ref: Wurster 2016, PASA, 34 (arXiv: 1608.00983), eq. 44
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void EoS_SetAuxArray_MultiGamma( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[0] = 2.0 * MULTIGAMMA_GAMMA1;
   AuxArray_Flt[1] =       MULTIGAMMA_GAMMA2;
   AuxArray_Flt[2] =       MULTIGAMMA_GAMMA3;
   AuxArray_Flt[3] = MULTIGAMMA_TURNOVER1 * MOLECULAR_WEIGHT * MU_NORM / UNIT_D;
   AuxArray_Flt[4] = MULTIGAMMA_TURNOVER2 * MOLECULAR_WEIGHT * MU_NORM / UNIT_D;
   AuxArray_Flt[5] = MULTIGAMMA_TURNOVER3 * MOLECULAR_WEIGHT * MU_NORM / UNIT_D;
   AuxArray_Flt[6] = MULTIGAMMA_TEMPBASE;
   AuxArray_Flt[7] = ( OPT__UNIT ) ? Const_kB * UNIT_M / ( MOLECULAR_WEIGHT * MU_NORM * UNIT_E )
                                   : MULTIGAMMA_TEMPBASE / MOLECULAR_WEIGHT;

} // FUNCTION : EoS_SetAuxArray_MultiGamma
#endif // #ifndef __CUDACC__



// =============================================
// II. Implement EoS conversion functions
//     (1) EoS_DensEint2Pres_*
//     (2) EoS_DensPres2Eint_*
//     (3) EoS_DensPres2CSqr_*
//     (4) EoS_DensEint2Temp_* [OPTIONAL]
//     (5) EoS_DensTemp2Pres_* [OPTIONAL]
//     (6) EoS_DensEint2Entr_* [OPTIONAL]
//     (7) EoS_General_*       [OPTIONAL]
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Pres_MultiGamma
// Description :  Convert gas mass density and internal energy density to gas pressure
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_MultiGamma() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Pres_MultiGamma( const real Dens, const real Eint, const real Passive[],
                                          const double AuxArray_Flt[], const int AuxArray_Int[],
                                          const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density", TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   const real _m_kB = (real)AuxArray_Flt[7];
   const real Temp  = EoS_DensEint2Temp_MultiGamma( Dens, Eint, Passive, AuxArray_Flt, AuxArray_Int, Table );
         real Pres  = Temp * Dens * _m_kB;

   return Pres;

} // FUNCTION : EoS_DensEint2Pres_MultiGamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2Eint_MultiGamma
// Description :  Convert gas mass density and pressure to gas internal energy density
//
// Note        :  1. See EoS_DensEint2Pres_MultiGamma()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas internal energy density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2Eint_MultiGamma( const real Dens, const real Pres, const real Passive[],
                                          const double AuxArray_Flt[], const int AuxArray_Int[],
                                          const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density",  TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( Pres, "input pressure", (real)0.0,   HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   real In_Flt[1] = { Dens };
   real Out[1];

   EoS_General_MultiGamma( NULL_INT, Out, In_Flt, NULL, AuxArray_Flt, AuxArray_Int, Table );

   const real _Gamma_m1 = (real)1.0 / ( Out[0] - (real)1.0 );
         real Eint      = Pres * _Gamma_m1;

   return Eint;

} // FUNCTION : EoS_DensPres2Eint_MultiGamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2CSqr_MultiGamma
// Description :  Convert gas mass density and pressure to sound speed squared
//
// Note        :  1. See EoS_DensEint2Pres_MultiGamma()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Sound speed squared
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2CSqr_MultiGamma( const real Dens, const real Pres, const real Passive[],
                                          const double AuxArray_Flt[], const int AuxArray_Int[],
                                          const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density",  TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( Pres, "input pressure", (real)0.0,   HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   real In_Flt[1] = { Dens };
   real Out[1];

   EoS_General_MultiGamma( NULL_INT, Out, In_Flt, NULL, AuxArray_Flt, AuxArray_Int, Table );

   const real Gamma = Out[0];
         real Cs2   = Gamma * Pres / Dens;

   return Cs2;

} // FUNCTION : EoS_DensPres2CSqr_MultiGamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Temp_MultiGamma
// Description :  Convert gas mass density and internal energy density to gas temperature
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_MultiGamma() for the values stored in AuxArray_Flt/Int[]
//                3. Temperature is in kelvin
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas temperature in kelvin
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Temp_MultiGamma( const real Dens, const real Eint, const real Passive[],
                                          const double AuxArray_Flt[], const int AuxArray_Int[],
                                          const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density", TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   const real PowIndex_1     = (real)AuxArray_Flt[0];
   const real PowIndex_2     = (real)AuxArray_Flt[1];
   const real PowIndex_3     = (real)AuxArray_Flt[2];
   const real Turnover_Dens1 = (real)AuxArray_Flt[3];
   const real Turnover_Dens2 = (real)AuxArray_Flt[4];
   const real Turnover_Dens3 = (real)AuxArray_Flt[5];
   const real Temp_Base      = (real)AuxArray_Flt[6];
         real Temp;

   Temp = Temp_Base
        * SQRT( (real)1.0 + POW( Dens / Turnover_Dens1, PowIndex_1 )  )
        *  POW( (real)1.0 +      Dens / Turnover_Dens2, PowIndex_2    )
        *  POW( (real)1.0 +      Dens / Turnover_Dens3, PowIndex_3    );

   return Temp;

} // FUNCTION : EoS_DensEint2Temp_MultiGamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensTemp2Pres_MultiGamma
// Description :  Convert gas mass density and temperature to gas pressure
//
// Note        :  1. See EoS_SetAuxArray_MultiGamma() for the values stored in AuxArray_Flt/Int[]
//                2. Temperature is in kelvin
//
// Parameter   :  Dens       : Gas mass density
//                Temp       : Gas temperature in kelvin
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensTemp2Pres_MultiGamma( const real Dens, const real Temp, const real Passive[],
                                          const double AuxArray_Flt[], const int AuxArray_Int[],
                                          const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density",     TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( Temp, "input temperature", (real)0.0,   HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   const real _m_kB = (real)AuxArray_Flt[7];
         real Pres  = Dens * Temp * _m_kB;

   return Pres;

} // FUNCTION : EoS_DensTemp2Pres_MultiGamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Entr_MultiGamma
// Description :  Convert gas mass density and internal energy density to gas entropy
//                --> Here entropy is defined as "pressure / density^(Gamma-1)" (i.e., entropy per volume)
//
// Note        :  1. See EoS_SetAuxArray_MultiGamma() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas entropy
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Entr_MultiGamma( const real Dens, const real Eint, const real Passive[],
                                          const double AuxArray_Flt[], const int AuxArray_Int[],
                                          const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_IsUnphysical_Single( Dens, "input density", TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   real In_Flt[1] = { Dens };
   real Out[1];

   EoS_General_MultiGamma( NULL_INT, Out, In_Flt, NULL, AuxArray_Flt, AuxArray_Int, Table );

   const real Gamma_m1 = Out[0] - (real)1.0;
         real Pres     = EoS_DensEint2Pres_MultiGamma( Dens, Eint, Passive, AuxArray_Flt, AuxArray_Int, Table );
         real Entr     = Pres * POW( Dens, -Gamma_m1 );

   return Entr;

} // FUNCTION : EoS_DensEint2Entr_MultiGamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_General_MultiGamma
// Description :  General EoS converter: In_*[] -> Out[]
//
// Note        :  1. See EoS_DensEint2Pres_MultiGamma()
//                2. In_*[] and Out[] must NOT overlap
//                3. Compute the adiabatic index for given density
//                   --> In_Flt[0] = mass density in code unit
//                       Out   [0] = adiabatic index
//
// Parameter   :  Mode       : To support multiple modes in this general converter
//                Out        : Output array
//                In_*       : Input array
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Out[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void EoS_General_MultiGamma( const int Mode, real Out[], const real In_Flt[], const int In_Int[],
                                    const double AuxArray_Flt[], const int AuxArray_Int[],
                                    const real *const Table[EOS_NTABLE_MAX] )
{

   real Dens = In_Flt[0];

   const real PowIndex_1     = (real)AuxArray_Flt[0];
   const real PowIndex_2     = (real)AuxArray_Flt[1];
   const real PowIndex_3     = (real)AuxArray_Flt[2];
   const real Turnover_Dens1 = (real)AuxArray_Flt[3];
   const real Turnover_Dens2 = (real)AuxArray_Flt[4];
   const real Turnover_Dens3 = (real)AuxArray_Flt[5];
   const real Term_1         = POW( Dens / Turnover_Dens1, PowIndex_1 );
   const real Term_2         = Dens / Turnover_Dens2;
   const real Term_3         = Dens / Turnover_Dens3;
         real Adia_Index;


   // compute the adiabatic index using gamma = 1.0 + (rho / T) * (dT / drho)
   Adia_Index = (real)1.0
              + PowIndex_1 * (  Term_1 / ( (real)1.0 + Term_1 )  )
              + PowIndex_2 * (  Term_2 / ( (real)1.0 + Term_2 )  )
              + PowIndex_3 * (  Term_3 / ( (real)1.0 + Term_3 )  );

   Out[0] = Adia_Index;

} // FUNCTION : EoS_General_MultiGamma



// =============================================
// III. Set EoS initialization functions
// =============================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE EoS_DE2P_t EoS_DensEint2Pres_Ptr = EoS_DensEint2Pres_MultiGamma;
FUNC_SPACE EoS_DP2E_t EoS_DensPres2Eint_Ptr = EoS_DensPres2Eint_MultiGamma;
FUNC_SPACE EoS_DP2C_t EoS_DensPres2CSqr_Ptr = EoS_DensPres2CSqr_MultiGamma;
FUNC_SPACE EoS_DE2T_t EoS_DensEint2Temp_Ptr = EoS_DensEint2Temp_MultiGamma;
FUNC_SPACE EoS_DT2P_t EoS_DensTemp2Pres_Ptr = EoS_DensTemp2Pres_MultiGamma;
FUNC_SPACE EoS_DE2S_t EoS_DensEint2Entr_Ptr = EoS_DensEint2Entr_MultiGamma;
FUNC_SPACE EoS_GENE_t EoS_General_Ptr       = EoS_General_MultiGamma;

//-----------------------------------------------------------------------------------------
// Function    :  EoS_SetCPU/GPUFunc_MultiGamma
// Description :  Return the function pointers of the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init_MultiGamma()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      EoS_SetFunc_MultiGamma( CPU_FuncPtr, GPU_FuncPtr );
//
//                3. Call-by-reference
//
// Parameter   :  EoS_DensEint2Pres_CPU/GPUPtr : CPU/GPU function pointers to be set
//                EoS_DensPres2Eint_CPU/GPUPtr : ...
//                EoS_DensPres2CSqr_CPU/GPUPtr : ...
//                EoS_DensEint2Temp_CPU/GPUPtr : ...
//                EoS_DensTemp2Pres_CPU/GPUPtr : ...
//                EoS_DensEint2Entr_CPU/GPUPtr : ...
//                EoS_General_CPU/GPUPtr       : ...
//
// Return      :  EoS_DensEint2Pres_CPU/GPUPtr, EoS_DensPres2Eint_CPU/GPUPtr,
//                EoS_DensPres2CSqr_CPU/GPUPtr, EoS_DensEint2Temp_CPU/GPUPtr,
//                EoS_DensTemp2Pres_CPU/GPUPtr, EoS_DensEint2Entr_CPU/GPUPtr,
//                EoS_General_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void EoS_SetGPUFunc_MultiGamma( EoS_DE2P_t &EoS_DensEint2Pres_GPUPtr,
                                EoS_DP2E_t &EoS_DensPres2Eint_GPUPtr,
                                EoS_DP2C_t &EoS_DensPres2CSqr_GPUPtr,
                                EoS_DE2T_t &EoS_DensEint2Temp_GPUPtr,
                                EoS_DT2P_t &EoS_DensTemp2Pres_GPUPtr,
                                EoS_DE2S_t &EoS_DensEint2Entr_GPUPtr,
                                EoS_GENE_t &EoS_General_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Pres_GPUPtr, EoS_DensEint2Pres_Ptr, sizeof(EoS_DE2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2Eint_GPUPtr, EoS_DensPres2Eint_Ptr, sizeof(EoS_DP2E_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2CSqr_GPUPtr, EoS_DensPres2CSqr_Ptr, sizeof(EoS_DP2C_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Temp_GPUPtr, EoS_DensEint2Temp_Ptr, sizeof(EoS_DE2T_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensTemp2Pres_GPUPtr, EoS_DensTemp2Pres_Ptr, sizeof(EoS_DT2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Entr_GPUPtr, EoS_DensEint2Entr_Ptr, sizeof(EoS_DE2S_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_General_GPUPtr,       EoS_General_Ptr,       sizeof(EoS_GENE_t) )  );
}

#else // #ifdef __CUDACC__

void EoS_SetCPUFunc_MultiGamma( EoS_DE2P_t &EoS_DensEint2Pres_CPUPtr,
                               EoS_DP2E_t &EoS_DensPres2Eint_CPUPtr,
                               EoS_DP2C_t &EoS_DensPres2CSqr_CPUPtr,
                               EoS_DE2T_t &EoS_DensEint2Temp_CPUPtr,
                               EoS_DT2P_t &EoS_DensTemp2Pres_CPUPtr,
                               EoS_DE2S_t &EoS_DensEint2Entr_CPUPtr,
                               EoS_GENE_t &EoS_General_CPUPtr )
{
   EoS_DensEint2Pres_CPUPtr = EoS_DensEint2Pres_Ptr;
   EoS_DensPres2Eint_CPUPtr = EoS_DensPres2Eint_Ptr;
   EoS_DensPres2CSqr_CPUPtr = EoS_DensPres2CSqr_Ptr;
   EoS_DensEint2Temp_CPUPtr = EoS_DensEint2Temp_Ptr;
   EoS_DensTemp2Pres_CPUPtr = EoS_DensTemp2Pres_Ptr;
   EoS_DensEint2Entr_CPUPtr = EoS_DensEint2Entr_Ptr;
   EoS_General_CPUPtr       = EoS_General_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void EoS_SetAuxArray_MultiGamma( double [], int [] );
void EoS_SetCPUFunc_MultiGamma( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t &, EoS_DE2T_t &, EoS_DT2P_t &, EoS_DE2S_t &, EoS_GENE_t & );
#ifdef GPU
void EoS_SetGPUFunc_MultiGamma( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t &, EoS_DE2T_t &, EoS_DT2P_t &, EoS_DE2S_t &, EoS_GENE_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  EoS_Init_MultiGamma
// Description :  Initialize EoS
//
// Note        :  1. Set auxiliary arrays by invoking EoS_SetAuxArray_*()
//                   --> It will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU EoS routines by invoking EoS_SetCPU/GPUFunc_*()
//                3. Invoked by EoS_Init()
//                   --> Enable it by linking to the function pointer "EoS_Init_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void EoS_Init_MultiGamma()
{

// check
#  ifndef BAROTROPIC_EOS
#  error : ERROR : must enable BAROTROPIC_EOS in the Makefile for the Multi_Gamma EoS !!
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );


   EoS_SetAuxArray_MultiGamma( EoS_AuxArray_Flt, EoS_AuxArray_Int );
   EoS_SetCPUFunc_MultiGamma( EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                              EoS_DensPres2CSqr_CPUPtr, EoS_DensEint2Temp_CPUPtr,
                              EoS_DensTemp2Pres_CPUPtr, EoS_DensEint2Entr_CPUPtr,
                              EoS_General_CPUPtr );
#  ifdef GPU
   EoS_SetGPUFunc_MultiGamma( EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr,
                              EoS_DensPres2CSqr_GPUPtr, EoS_DensEint2Temp_GPUPtr,
                              EoS_DensTemp2Pres_GPUPtr, EoS_DensEint2Entr_GPUPtr,
                              EoS_General_GPUPtr );
#  endif

} // FUNCTION : EoS_Init_MultiGamma

#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO )
