#include "GAMER.h"
#include "TestProb.h"
#include "NuclearEoS.h"



// problem-specific global variables
// =======================================================================================
typedef int CCSN_t;
const CCSN_t
   Migration_Test = 0
  ,Post_Bounce    = 1
  ;

typedef int CCSN_Mag_t;
const CCSN_Mag_t
   Liu2008  = 0
  ,Suwa2007 = 1
  ;

static CCSN_t     CCSN_Prob;                       // target CCSN problem
static char       CCSN_Name[100];                  // name of the target CCSN problem
static char       CCSN_Prof_File[MAX_STRING];      // filename of input profile
static double    *CCSN_Prof = NULL;                // radial profile of initial condition
static int        CCSN_Prof_NBin;                  // number of radial bins in the input profile
static int        CCSN_NCol;                       // number of columns read from the input profile
static int       *CCSN_TargetCols = new int [7];   // index of columns read from the input profile
static int        CCSN_ColIdx_Dens;                // column index of density in the input profile
static int        CCSN_ColIdx_Pres;                // column index of pressure in the input profile
static int        CCSN_ColIdx_Velr;                // column index of radial velocity in the input profile

#ifdef MHD
static CCSN_Mag_t CCSN_Mag;                        // target magnetic field profile
static double     CCSN_Mag_B0;                     // magnetic field strength
static double     CCSN_Mag_np;                     // dependence of magnetic field on density
static double     CCSN_Mag_R0;                     // characteristic radius of magnetic field
#endif

static bool       CCSN_GW_OUTPUT;                  // output GW signals (0=off, 1=on)
static double     CCSN_GW_DT;                      // output GW signals every CCSN_GW_DT interval

static int        CCSN_Eint_Mode;                  // Mode of obtaining internal energy in SetGridIC()
                                                   // ( 1=Temp Mode: Eint(dens, temp, [Ye])
                                                   //   2=Pres Mode: Eint(dens, pres, [Ye]) )

static double     CCSN_MaxRefine_RadFac;           // factor that determines the maximum refinement level based on distance from the box center
       double     CCSN_LB_TimeFac;                 // factor that scales the dt constrained by lightbulb scheme
// =======================================================================================


// problem-specific function prototypes
void Record_CCSN_CentralDens();
void Record_CCSN_GWSignal();
void Mis_GetTimeStep_User_Lightbulb( const int lv, const double dTime_dt );



//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "CCSN_Prob",             &CCSN_Prob,             -1,            0,                1                 );
   ReadPara->Add( "CCSN_Prof_File",         CCSN_Prof_File,        Useless_str,   Useless_str,      Useless_str       );
#  ifdef MHD
   ReadPara->Add( "CCSN_Mag",              &CCSN_Mag,              1,             0,                1                 );
   ReadPara->Add( "CCSN_Mag_B0",           &CCSN_Mag_B0,           1.0e14,        0.0,              NoMax_double      );
   ReadPara->Add( "CCSN_Mag_np",           &CCSN_Mag_np,           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "CCSN_Mag_R0",           &CCSN_Mag_R0,           1.0e8,         Eps_double,       NoMax_double      );
#  endif
   ReadPara->Add( "CCSN_GW_OUTPUT",        &CCSN_GW_OUTPUT,        false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "CCSN_GW_DT",            &CCSN_GW_DT,            1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "CCSN_Eint_Mode",        &CCSN_Eint_Mode,        2,             1,                2                 );
   ReadPara->Add( "CCSN_MaxRefine_RadFac", &CCSN_MaxRefine_RadFac, 0.15,          0.0,              NoMax_double      );
   ReadPara->Add( "CCSN_LB_TimeFac",       &CCSN_LB_TimeFac,       0.1,           Eps_double,       1.0               );


   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   switch ( CCSN_Prob )
   {
      case Migration_Test : CCSN_NCol = 4;
                            CCSN_TargetCols[0] =  0;  CCSN_TargetCols[1] =  1;  CCSN_TargetCols[2] =  2;  CCSN_TargetCols[3] =  3;
                            CCSN_TargetCols[4] = -1;  CCSN_TargetCols[5] = -1;  CCSN_TargetCols[6] = -1;
                            CCSN_ColIdx_Dens   =  2;  CCSN_ColIdx_Pres   =  3;  CCSN_ColIdx_Velr   =  1;
                            sprintf( CCSN_Name, "GREP migration test" );
                            break;

      case Post_Bounce    : CCSN_NCol = 7;
                            CCSN_TargetCols[0] =  0;  CCSN_TargetCols[1] =  2;  CCSN_TargetCols[2] =  3;  CCSN_TargetCols[3] =  4;
                            CCSN_TargetCols[4] =  5;  CCSN_TargetCols[5] =  6;  CCSN_TargetCols[6] =  7;
                            CCSN_ColIdx_Dens   =  1;  CCSN_ColIdx_Pres   =  5;  CCSN_ColIdx_Velr   =  3;
                            sprintf( CCSN_Name, "Post bounce test" );
                            break;

      default             : Aux_Error( ERROR_INFO, "unsupported CCSN problem (%d) !!\n", CCSN_Prob );
   } // switch ( CCSN_Prob )

// (1-3) check the runtime parameters
   if (  ( CCSN_Prob == Migration_Test )  &&  ( CCSN_Eint_Mode == 1 )  )
      Aux_Error( ERROR_INFO, "Temperature mode for obtaining internal energy is not supported in Migration Test yet!!\n" );

#  if ( EOS != EOS_NUCLEAR )
   if ( CCSN_Eint_Mode == 0 )
      Aux_Error( ERROR_INFO, "Temperature mode only works with EOS_NUCLEAR !!\n" );
#  endif

   if (  OPT__DT_USER  &&  ( SrcTerms.Lightbulb == 0 )  )
      Aux_Error( ERROR_INFO, "OPT__DT_USER only supports SRC_LIGHTBULB == 1 !!\n" );


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                         = %d\n",      TESTPROB_ID           );
      Aux_Message( stdout, "  target CCSN problem                     = %s\n",      CCSN_Name             );
      Aux_Message( stdout, "  initial profile                         = %s\n",      CCSN_Prof_File        );
#     ifdef MHD
      Aux_Message( stdout, "  magnetic field profile                  = %d\n",      CCSN_Mag              );
      Aux_Message( stdout, "  magnetic field strength                 = %13.7e\n",  CCSN_Mag_B0           );
      Aux_Message( stdout, "  dependence of magnetic field on density = %13.7e\n",  CCSN_Mag_np           );
      Aux_Message( stdout, "  characteristic radius of magnetic field = %13.7e\n",  CCSN_Mag_R0           );
#     endif
      Aux_Message( stdout, "  output GW signals                       = %d\n",      CCSN_GW_OUTPUT        );
      Aux_Message( stdout, "  sampling interval of GW signals         = %13.7e\n",  CCSN_GW_DT            );
      Aux_Message( stdout, "  mode for obtaining internal energy      = %d\n",      CCSN_Eint_Mode        );
      Aux_Message( stdout, "  radial factor for maximum refine level  = %13.7e\n",  CCSN_MaxRefine_RadFac );
      Aux_Message( stdout, "  scaling factor for lightbulb dt         = %13.7e\n",  CCSN_LB_TimeFac       );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const double  BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   const double *Table_R      = CCSN_Prof +                0*CCSN_Prof_NBin;
   const double *Table_Velr   = CCSN_Prof + CCSN_ColIdx_Velr*CCSN_Prof_NBin;
   const double *Table_Dens   = CCSN_Prof + CCSN_ColIdx_Dens*CCSN_Prof_NBin;
   const double *Table_Pres   = CCSN_Prof + CCSN_ColIdx_Pres*CCSN_Prof_NBin;

   const double x0 = x - BoxCenter[0];
   const double y0 = y - BoxCenter[1];
   const double z0 = z - BoxCenter[2];
   const double r  = sqrt(  SQR( x0 ) + SQR( y0 ) + SQR( z0 )  );

   double Dens, Velr, Pres, Momx, Momy, Momz, Eint, Etot, Ye, Temp, Entr;

   Dens = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, Table_Dens, r );
   Velr = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, Table_Velr, r );
   Pres = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, Table_Pres, r );

   if ( Dens == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed for density at radius %13.7e !!\n", r );
   if ( Velr == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed for radial velocity at radius %13.7e !!\n", r );
   if ( Pres == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed for pressure at radius %13.7e !!\n", r );

   if ( CCSN_Prob == Post_Bounce )
   {
      Ye   = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, CCSN_Prof+2*CCSN_Prof_NBin, r );
      Temp = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, CCSN_Prof+4*CCSN_Prof_NBin, r );  // in Kelvin
      Entr = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, CCSN_Prof+6*CCSN_Prof_NBin, r );

      if ( Ye   == NULL_REAL )
         Aux_Error( ERROR_INFO, "interpolation failed for Ye at radius %13.7e !!\n", r );
      if ( Temp == NULL_REAL )
         Aux_Error( ERROR_INFO, "interpolation failed for temperature at radius %13.7e !!\n", r );
      if ( Entr == NULL_REAL )
         Aux_Error( ERROR_INFO, "interpolation failed for entropy at radius %13.7e !!\n", r );
   }


   Momx = Dens*Velr*x0/r;
   Momy = Dens*Velr*y0/r;
   Momz = Dens*Velr*z0/r;

// calculate the internal energy
#  if ( EOS == EOS_NUCLEAR )
   real *Passive = new real [NCOMP_PASSIVE];

   Passive[ YE - NCOMP_FLUID ] = Ye*Dens;
#  else
   real *Passive = NULL;
#  endif

   if ( CCSN_Eint_Mode == 1 )   // Temperature Mode
   {
#     if ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
      const int  NTarget = 1;
#     else
      const int  NTarget = 0;
#     endif
            int  TmpIn_Int[NTarget+1];
            real TmpIn_Flt[3], TmpOut[NTarget+1];

      TmpIn_Flt[0] = (real)Dens;
      TmpIn_Flt[1] = (real)Temp;
      TmpIn_Flt[2] = (real)Ye;

      TmpIn_Int[0] = NTarget;
#     if ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
      TmpIn_Int[1] = NUC_VAR_IDX_EORT;
#     endif

#     if ( EOS == EOS_NUCLEAR )
      EoS_General_CPUPtr( NUC_MODE_TEMP, TmpOut, TmpIn_Flt, TmpIn_Int, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#     endif
      Eint = TmpOut[0];
   }

   else                         // Pressure Mode
   {
      Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   }


   Etot = Hydro_ConEint2Etot( Dens, Momx, Momy, Momz, Eint, 0.0 );   // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = Momx;
   fluid[MOMY] = Momy;
   fluid[MOMZ] = Momz;
   fluid[ENGY] = Etot;
#  if ( EOS == EOS_NUCLEAR )
   for (int v=0; v<NCOMP_PASSIVE; v++)   fluid[ NCOMP_FLUID + v ] = Passive[v];
#  endif


   if ( Passive != NULL )   delete [] Passive;

} // FUNCTION : SetGridIC



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC_Liu2008
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                2. Generate the poloidal B field from a variant form of the vector potential in
//                   Liu+ 2008, Phys. Rev. D78, 024012:
//                       A_phi = B0 * varpi^2 * (1 - rho / rho_max)^np * (P / P_max)
//                       A_r = A_theta = 0
//                   where
//                       varpi^2 =  x^2 + y^2
//                       A_x      = -(y / varpi^2) * A_phi
//                       A_y      =  (x / varpi^2) * A_phi
//                       A_z      =  0
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void SetBFieldIC_Liu2008( real magnetic[], const double x, const double y, const double z, const double Time,
                          const int lv, double AuxArray[] )
{

   const double  BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   const double *Table_R      = CCSN_Prof +                0*CCSN_Prof_NBin;
   const double *Table_Dens   = CCSN_Prof + CCSN_ColIdx_Dens*CCSN_Prof_NBin;
   const double *Table_Pres   = CCSN_Prof + CCSN_ColIdx_Pres*CCSN_Prof_NBin;

   const double x0 = x - BoxCenter[0];
   const double y0 = y - BoxCenter[1];
   const double z0 = z - BoxCenter[2];

// approximate the central density and pressure by the data at the first row
   const double dens_c = Table_Dens[0];
   const double pres_c = Table_Pres[0];
   const double B0     = CCSN_Mag_B0 / UNIT_B;

// use finite difference to compute the B field
   double delta = amr->dh[MAX_LEVEL];
   double r,    dens,    pres;
   double r_xp, dens_xp, pres_xp;
   double r_yp, dens_yp, pres_yp;
   double r_zp, dens_zp, pres_zp;

   r    = sqrt(  SQR( x0         ) + SQR( y0         ) + SQR( z0         )  );
   r_xp = sqrt(  SQR( x0 + delta ) + SQR( y0         ) + SQR( z0         )  );
   r_yp = sqrt(  SQR( x0         ) + SQR( y0 + delta ) + SQR( z0         )  );
   r_zp = sqrt(  SQR( x0         ) + SQR( y0         ) + SQR( z0 + delta )  );

   dens    = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, Table_Dens, r    );
   dens_xp = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, Table_Dens, r_xp );
   dens_yp = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, Table_Dens, r_yp );
   dens_zp = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, Table_Dens, r_zp );

   if ( dens    == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed for dens    at radius %13.7e !!\n", r    );
   if ( dens_xp == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed for dens_xp at radius %13.7e !!\n", r_xp );
   if ( dens_yp == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed for dens_yp at radius %13.7e !!\n", r_yp );
   if ( dens_zp == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed for dens_zp at radius %13.7e !!\n", r_zp );

   pres    = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, Table_Pres, r    );
   pres_xp = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, Table_Pres, r_xp );
   pres_yp = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, Table_Pres, r_yp );
   pres_zp = Mis_InterpolateFromTable( CCSN_Prof_NBin, Table_R, Table_Pres, r_zp );

   if ( pres    == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed for pres    at radius %13.7e !!\n", r    );
   if ( pres_xp == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed for pres_xp at radius %13.7e !!\n", r_xp );
   if ( pres_yp == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed for pres_yp at radius %13.7e !!\n", r_yp );
   if ( pres_zp == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed for pres_zp at radius %13.7e !!\n", r_zp );


   double dAy_dx = (  ( x0 + delta )*POW( 1.0 - dens_xp/dens_c, CCSN_Mag_np )*( pres_xp / pres_c )   \
                   -  ( x0         )*POW( 1.0 - dens   /dens_c, CCSN_Mag_np )*( pres    / pres_c ) ) \
                 / delta;

   double dAx_dy = ( -( y0 + delta )*POW( 1.0 - dens_yp/dens_c, CCSN_Mag_np )*( pres_yp / pres_c )   \
                   - -( y0         )*POW( 1.0 - dens   /dens_c, CCSN_Mag_np )*( pres    / pres_c ) ) \
                 / delta;

   double dAphi_dz = ( POW( 1.0 - dens_zp/dens_c, CCSN_Mag_np )*( pres_zp / pres_c )   \
                     - POW( 1.0 - dens   /dens_c, CCSN_Mag_np )*( pres    / pres_c ) ) \
                   / delta;

   magnetic[MAGX] = -x0 * B0 * dAphi_dz;
   magnetic[MAGY] = -y0 * B0 * dAphi_dz;
   magnetic[MAGZ] =       B0 * ( dAy_dx - dAx_dy );

} // FUNCTION : SetBFieldIC_Liu2008



//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC_Suwa2007
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                2. Generate the poloidal B field from the vector potential in
//                   Suwa+ 2007, PASJ, 59, 771:
//                       A_phi = 0.5 * B0 * ( R0^3 / (r^3 + R0^3) ) * r * sin(theta)
//                       A_r = A_theta = 0
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void SetBFieldIC_Suwa2007( real magnetic[], const double x, const double y, const double z, const double Time,
                           const int lv, double AuxArray[] )
{

   const double  BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };

   const double x0 = x - BoxCenter[0];
   const double y0 = y - BoxCenter[1];
   const double z0 = z - BoxCenter[2];
   const double r  = sqrt(  SQR( x0 ) + SQR( y0 ) + SQR( z0 )  );

   const double B0     = CCSN_Mag_B0 / UNIT_B;
   const double R0     = CCSN_Mag_R0 / UNIT_L;
   const double R0_cub = CUBE( R0 );
   const double frac   = 1.0 / (  1.0 + CUBE( r / R0 )  );

   magnetic[MAGX] = 1.5 * B0 * SQR( frac ) * ( r * x0 * z0               ) / R0_cub;
   magnetic[MAGY] = 1.5 * B0 * SQR( frac ) * ( r * y0 * z0               ) / R0_cub;
   magnetic[MAGZ] =       B0 *      frac
                  - 1.5 * B0 * SQR( frac ) * ( r * x0 * x0 + r * y0 * y0 ) / R0_cub;

} // FUNCTION : SetBFieldIC_Suwa2007
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_IC_Prof_CCSN
// Description :  Load input table file for initial condition
//-------------------------------------------------------------------------------------------------------
void Load_IC_Prof_CCSN()
{

   const bool RowMajor_No  = false;           // load data into the column major
   const bool AllocMem_Yes = true;            // allocate memory for CCSN_Prof

   CCSN_Prof_NBin = Aux_LoadTable( CCSN_Prof, CCSN_Prof_File, CCSN_NCol, CCSN_TargetCols, RowMajor_No, AllocMem_Yes );

// convert radius, density, radial velocity, and pressure to code units
   double *Table_R    = CCSN_Prof +                0*CCSN_Prof_NBin;
   double *Table_Velr = CCSN_Prof + CCSN_ColIdx_Velr*CCSN_Prof_NBin;
   double *Table_Dens = CCSN_Prof + CCSN_ColIdx_Dens*CCSN_Prof_NBin;
   double *Table_Pres = CCSN_Prof + CCSN_ColIdx_Pres*CCSN_Prof_NBin;

   for (int b=0; b<CCSN_Prof_NBin; b++)
   {
      Table_R   [b] /= UNIT_L;
      Table_Dens[b] /= UNIT_D;
      Table_Velr[b] /= UNIT_V;
      Table_Pres[b] /= UNIT_P;
   }

} // FUNCTION : Load_IC_Prof_CCSN()



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_CCSN
// Description :  Interface for invoking multiple record functions
//-------------------------------------------------------------------------------------------------------
void Record_CCSN()
{

// 1. record the maximum density
   Record_CCSN_CentralDens();


#ifdef GRAVITY
// 2. record the GW signal
   if ( CCSN_GW_OUTPUT )
   {
      static double DumpTime_Last = -__DBL_MAX__;

             double PhysicalTime = Time[0];
             double DumpTime = int( PhysicalTime / CCSN_GW_DT ) * CCSN_GW_DT;

//    output data when (a) the elapsed time > GW_OUTPUT_DT
//                     (b) the physical time is almost a multiple of GW_OUTPUT_DT
      if ( fabs( PhysicalTime - DumpTime_Last ) > CCSN_GW_DT  ||
           fabs( PhysicalTime - DumpTime      ) < 1.0e-2 * CCSN_GW_DT )
      {
         DumpTime_Last = PhysicalTime;
         Record_CCSN_GWSignal();
      }
   }
#endif

} // FUNCTION : Record_CCSN()



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_User_CCSN
// Description :  Check if the element (i,j,k) of the input data satisfies the user-defined flag criteria
//
// Note        :  1. Invoked by "Flag_Check" using the function pointer "Flag_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this function will become useless
//                2. Enabled by the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//                Threshold   : User-provided threshold for the flag operation, which is loaded from the
//                              file "Input__Flag_User"
//                              In order of radius_min, radius_max, threshold_dens
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_User_CCSN( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   bool Flag = false;

   const double dh        = amr->dh[lv];
   const double Center[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   const double Pos   [3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   const double dx = Center[0] - Pos[0];
   const double dy = Center[1] - Pos[1];
   const double dz = Center[2] - Pos[2];
   const double r  = sqrt(  SQR( dx ) + SQR( dy ) + SQR( dz )  );

   const real (*Rho )[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS];


// TODO: fine-tune the criteria
// always allow the highest level can be reached in the region with r < 30 km
   if ( r * UNIT_L < 3e6 )
   {
      Flag = true;
   }

   else
   {
//    if density is larger than the threshold (equivalent to Input__Flag_Rho)
//    and the cell width at lv+1 is larger than the threshold `r * CCSN_MaxRefine_RadFac` (equivalent to Flag_Region)
      if (  ( Rho[k][j][i] > Threshold[0] )  &&  ( 0.5 * dh > r * CCSN_MaxRefine_RadFac )  )   Flag = true;
   }

   return Flag;

} // FUNCTION : Flag_User_CCSN



//-------------------------------------------------------------------------------------------------------
// Function    :  End_CCSN
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_CCSN()
{

   delete [] CCSN_Prof;         CCSN_Prof       = NULL;
   delete [] CCSN_TargetCols;   CCSN_TargetCols = NULL;

} // FUNCTION : End_CCSN



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_CCSN
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_CCSN()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();

// load initial condition from file
   if ( OPT__INIT != INIT_BY_RESTART )   Load_IC_Prof_CCSN();

// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr         = SetGridIC;
#  ifdef MHD
   switch ( CCSN_Mag )
   {
      case Liu2008  : Init_Function_BField_User_Ptr  = SetBFieldIC_Liu2008;    break;
      case Suwa2007 : Init_Function_BField_User_Ptr  = SetBFieldIC_Suwa2007;   break;
   }
#  endif
   Flag_User_Ptr                  = Flag_User_CCSN;
   Aux_Record_User_Ptr            = Record_CCSN;
   End_User_Ptr                   = End_CCSN;
// estimate the evolution time-step constrained by the lightbulb source term
#  if ( EOS == EOS_NUCLEAR )
   if ( SrcTerms.Lightbulb )
   {
      Mis_GetTimeStep_User_Ptr    = Mis_GetTimeStep_User_Lightbulb;  // option: OPT__DT_USER;
   }
#  endif // #if ( EOS == EOS_NUCLEAR )
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_CCSN
