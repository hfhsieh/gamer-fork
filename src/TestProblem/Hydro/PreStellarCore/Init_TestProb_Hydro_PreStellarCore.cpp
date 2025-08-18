#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double PSC_Mass;             // mass   of the spherical molecular cloud, in solar mass
static double PSC_Radius;           // radius of the spherical molecular cloud, in AU
static double PSC_Ratio_Dens_Edge;  // ratio of central to edge    density
static double PSC_Ratio_Dens_Env;   // ratio of edge to background density
static double PSC_Ratio_Energy;     // ratio of rotational to gravitational energy

       double PSC_Mass_Code;        // mass   of the spherical molecular cloud, in code unit
       double PSC_Radius_Code;      // radius of the spherical molecular cloud, in code unit
       double PSC_RhoBase_Code;     // density of central plateau, in code unit
       double PSC_RadBase_Code;     // radius  of central plateau, in code unit
       double PSC_RhoEnv_Code;      // background density, in code unit
       double PSC_OmegaBase;        // angular velocity of the spherical molecular cloud

       bool   PSC_Prof;             // output spherically averaged profile at each global step
       int    PSC_Prof_Center;      // center of spherically averaged profile
       bool   PSC_Prof_LogBin;      // log/linear bins
       double PSC_Prof_LogBinRatio; // ratio of adjacent log bins
       double PSC_Prof_MaxRadius;   // maximum radius in the radial profile, in code units
       double PSC_Prof_MinBinSize;  // minimum bin size, in code units
// =======================================================================================




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

#  if ( EOS != EOS_MULTIGAMMA )
   Aux_Error( ERROR_INFO, "EOS != EOS_MULTIGAMMA !!\n" );
#  endif

#  ifndef BAROTROPIC_EOS
   Aux_Error( ERROR_INFO, "BAROTROPIC_EOS must be enabled !!\n" );
#  endif


   if ( MPI_Rank == 0 )
   {
      if ( !OPT__FLAG_JEANS )
         Aux_Error( ERROR_INFO, "OPT__FLAG_JEANS must be enabled !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  LoadInputTestProb
// Description :  Read problem-specific runtime parameters from Input__TestProb and store them in HDF5 snapshots (Data_*)
//
// Note        :  1. Invoked by SetParameter() to read parameters
//                2. Invoked by Output_DumpData_Total_HDF5() using the function pointer Output_HDF5_InputTest_Ptr to store parameters
//                3. If there is no problem-specific runtime parameter to load, add at least one parameter
//                   to prevent an empty structure in HDF5_Output_t
//                   --> Example:
//                       LOAD_PARA( load_mode, "TestProb_ID", &TESTPROB_ID, TESTPROB_ID, TESTPROB_ID, TESTPROB_ID );
//
// Parameter   :  load_mode      : Mode for loading parameters
//                                 --> LOAD_READPARA    : Read parameters from Input__TestProb
//                                     LOAD_HDF5_OUTPUT : Store parameters in HDF5 snapshots
//                ReadPara       : Data structure for reading parameters (used with LOAD_READPARA)
//                HDF5_InputTest : Data structure for storing parameters in HDF5 snapshots (used with LOAD_HDF5_OUTPUT)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LoadInputTestProb( const LoadParaMode_t load_mode, ReadPara_t *ReadPara, HDF5_Output_t *HDF5_InputTest )
{

#  ifndef SUPPORT_HDF5
   if ( load_mode == LOAD_HDF5_OUTPUT )   Aux_Error( ERROR_INFO, "please turn on SUPPORT_HDF5 in the Makefile for load_mode == LOAD_HDF5_OUTPUT !!\n" );
#  endif

   if ( load_mode == LOAD_READPARA     &&  ReadPara       == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_READPARA and ReadPara == NULL !!\n" );
   if ( load_mode == LOAD_HDF5_OUTPUT  &&  HDF5_InputTest == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_HDF5_OUTPUT and HDF5_InputTest == NULL !!\n" );

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// --> LOAD_PARA() is defined in "include/TestProb.h"
// ********************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",        &VARIABLE,                  DEFAULT,      MIN,              MAX               );
// ********************************************************************************************************************************
   LOAD_PARA( load_mode, "PSC_Mass",               &PSC_Mass,                  1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "PSC_Radius",             &PSC_Radius,                5.0e3,        Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "PSC_Ratio_Dens_Edge",    &PSC_Ratio_Dens_Edge,       1.0e1,        2.0,              NoMax_double      );
   LOAD_PARA( load_mode, "PSC_Ratio_Dens_Env",     &PSC_Ratio_Dens_Env,        1.0e2,        1.0,              NoMax_double      );
   LOAD_PARA( load_mode, "PSC_Ratio_Energy",       &PSC_Ratio_Energy,          1.0e-2,       0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "PSC_Prof",               &PSC_Prof,                  false,        Useless_bool,     Useless_bool      );
   LOAD_PARA( load_mode, "PSC_Prof_Center",        &PSC_Prof_Center,           2,            1,                4                 );
   LOAD_PARA( load_mode, "PSC_Prof_LogBin",        &PSC_Prof_LogBin,           false,        Useless_bool,     Useless_bool      );
   LOAD_PARA( load_mode, "PSC_Prof_LogBinRatio",   &PSC_Prof_LogBinRatio,      1.01,         1.0,              NoMax_double      );
   LOAD_PARA( load_mode, "PSC_Prof_MaxRadius",     &PSC_Prof_MaxRadius,       -1.0,          NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "PSC_Prof_MinBinSize",    &PSC_Prof_MinBinSize,      -1.0,          NoMin_double,     NoMax_double      );

} // FUNCITON : LoadInputTestProb



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
//                3. Must call EoS_Init() before calling any other EoS routine
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

   LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
// unit conversion
   PSC_Mass_Code   = PSC_Mass * Const_Msun / UNIT_M;
   PSC_Radius_Code = PSC_Radius * Const_au / UNIT_L;

// compute the density and radius of central plateau
   const double PSC_Ratio_Rad_Edge = sqrt( PSC_Ratio_Dens_Edge - 1.0 );  // edge-to-central radius ratio

   PSC_RadBase_Code = PSC_Radius_Code / PSC_Ratio_Rad_Edge;
   PSC_RhoBase_Code = PSC_Mass_Code
                    / (  4.0 * M_PI * CUBE(PSC_RadBase_Code) * ( PSC_Ratio_Rad_Edge - atan(PSC_Ratio_Rad_Edge) )  );

// compute the background density
   PSC_RhoEnv_Code = PSC_RhoBase_Code / ( PSC_Ratio_Dens_Edge * PSC_Ratio_Dens_Env );

// compute the angular velocity
   PSC_OmegaBase = sqrt( 3.0 * PSC_Ratio_Energy * NEWTON_G * PSC_Mass_Code / CUBE(PSC_Radius_Code) );


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

   if ( PSC_Prof_MaxRadius < 0.0 ) {
      PSC_Prof_MaxRadius = sqrt(3.0) * amr->BoxSize[0];
      PRINT_RESET_PARA( END_T, PSC_Prof_MaxRadius, "" );
   }

   if ( PSC_Prof_MinBinSize < 0.0 ) {
      PSC_Prof_MinBinSize = amr->dh[MAX_LEVEL];
      PRINT_RESET_PARA( END_T, PSC_Prof_MinBinSize, "" );
   }

// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                                  = %d\n",      TESTPROB_ID               );
      Aux_Message( stdout, "  mass   of spherical molecular cloud (solar mass) = % 14.7e\n", PSC_Mass                  );
      Aux_Message( stdout, "  radius of spherical molecular cloud         (AU) = % 14.7e\n", PSC_Radius                );
      Aux_Message( stdout, "  ratio of central to edge    density              = % 14.7e\n", PSC_Ratio_Dens_Edge       );
      Aux_Message( stdout, "  ratio of edge to background density              = % 14.7e\n", PSC_Ratio_Dens_Env        );
      Aux_Message( stdout, "  ratio of edge to central    radius               = % 14.7e\n", PSC_Ratio_Rad_Edge        );
      Aux_Message( stdout, "  ratio of rotational to gravitational energy      = % 14.7e\n", PSC_Ratio_Energy          );
      Aux_Message( stdout, "  central plateau density                  (g/cm3) = % 14.7e\n", PSC_RhoBase_Code * UNIT_D );
      Aux_Message( stdout, "  background density                       (g/cm3) = % 14.7e\n", PSC_RhoEnv_Code  * UNIT_D );
      Aux_Message( stdout, "  angular velocity                         (rad/s) = % 14.7e\n", PSC_OmegaBase    / UNIT_T );
      Aux_Message( stdout, "  dump spherically average profile                 = % d\n",     PSC_Prof                  );
      Aux_Message( stdout, "  center of spherically average profile            = % d\n",     PSC_Prof_Center           );
      Aux_Message( stdout, "  log/linear bins in the profile                   = % d\n",     PSC_Prof_LogBin           );
      Aux_Message( stdout, "  ratio of adjacent log bins                       = % 14.7e\n", PSC_Prof_LogBinRatio      );
      Aux_Message( stdout, "  maximum radius in the profile                    = % 14.7e\n", PSC_Prof_MaxRadius        );
      Aux_Message( stdout, "  minimum bin size                                 = % 14.7e\n", PSC_Prof_MinBinSize       );
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
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
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

   const double x0 = x - amr->BoxCenter[0];
   const double y0 = y - amr->BoxCenter[1];
   const double z0 = z - amr->BoxCenter[2];
   const double r  = sqrt(  SQR( x0 ) + SQR( y0 ) + SQR( z0 )  );

   double Dens, MomX, MomY, MomZ, Temp, Pres, Eint, Etot;

// assume the molecular cloud undergoes rigid-body rotation along the z direction, if applicable
   if ( r > PSC_Radius_Code )
   {
      Dens = PSC_RhoEnv_Code;
      MomX = 0.0;
      MomY = 0.0;
      MomZ = 0.0;
   }

   else
   {
      Dens = PSC_RhoBase_Code / (  1.0 + SQR( r / PSC_RadBase_Code )  );
      MomX = -Dens * PSC_OmegaBase * y0;
      MomY =  Dens * PSC_OmegaBase * x0;
      MomZ = 0.0;
   }

// REVISE: support DensTemp2Eint in EoS_General
   Temp = EoS_DensEint2Temp_CPUPtr( Dens, NULL_REAL, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );
   Pres = EoS_DensTemp2Pres_CPUPtr( Dens, Temp,      NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres,      NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{

   /*
// example
   magnetic[MAGX] = 1.0;
   magnetic[MAGY] = 2.0;
   magnetic[MAGZ] = 3.0;
   */

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep_PSC
// Description :  Criteria to estimate the evolution time-step based on the local free-fall time.
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Invoked by Mis_GetTimeStep() using the function pointer "Mis_GetTimeStep_User_Ptr",
//                   which must be set by a test problem initializer
//                3. Enabled by the runtime option "OPT__DT_USER"
//
// Parameter   :  lv       : Target refinement level
//                dTime_dt : dTime/dt (== 1.0 if COMOVING is off)
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double Mis_GetTimeStep_PSC( const int lv, const double dTime_dt )
{

   const double Factor_FF      = 3.0 * M_PI / 32.0;
   const double Factor_Scaling = 0.1;

// allocate memory for per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;
#  else
   const int NT = 1;
#  endif

   double  dt_PSC     = HUGE_NUMBER;
   double *OMP_dt_PSC = new double [NT];


#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    initialize arrays
      OMP_dt_PSC[TID] = HUGE_NUMBER;

#     pragma omp for schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {

         for (int k=0; k<PS1; k++)  {
         for (int j=0; j<PS1; j++)  {
         for (int i=0; i<PS1; i++)  {

            const real   Dens        = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
            const double dt_ThisCell = SQRT(  Factor_FF / ( NEWTON_G * Dens )  );


            OMP_dt_PSC[TID] = FMIN( OMP_dt_PSC[TID], dt_ThisCell );

         }}} // i,j,k
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // OpenMP parallel region


// find the minimum over all OpenMP threads
   for (int TID=0; TID<NT; TID++)   dt_PSC = FMIN( dt_PSC, OMP_dt_PSC[TID] );

// free per-thread arrays
   delete [] OMP_dt_PSC;


// find the minimum over all MPI processes
#  ifndef SERIAL
   MPI_Allreduce( MPI_IN_PLACE, &dt_PSC, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
#  endif


   return Factor_Scaling * dt_PSC;

} // FUNCTION : Mis_GetTimeStep_PSC



//-------------------------------------------------------------------------------------------------------
// Function    :  Record_PSC
// Description :  Record the spherically averaged profiles
//-------------------------------------------------------------------------------------------------------
void Record_PSC()
{

   if ( PSC_Prof )
   {

//    (1-1) find the location of peak density
      double    Center[3];
      Extrema_t Extrema;

      switch ( PSC_Prof_Center )
      {
         case 1: // box center
         {
            for (int i=0; i<3; i++)   Center[i] = amr->BoxCenter[i];
         }
         break;

         case 2: // density maximum
         {
            Extrema.Field     = _DENS;
            Extrema.Radius    = __FLT_MAX__;
            Extrema.Center[0] = amr->BoxCenter[0];
            Extrema.Center[1] = amr->BoxCenter[1];
            Extrema.Center[2] = amr->BoxCenter[2];

            Aux_FindExtrema( &Extrema, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );

            for (int i=0; i<3; i++)   Center[i] = Extrema.Coord[i];

//          shift the center to the box center if it coincides with one of the innermost cells
            const double Extrema_dh = amr->dh[ Extrema.Level ];

            if (  fabs( Center[0] - amr->BoxCenter[0] ) < Extrema_dh  &&
                  fabs( Center[1] - amr->BoxCenter[1] ) < Extrema_dh  &&
                  fabs( Center[2] - amr->BoxCenter[2] ) < Extrema_dh    )
               for (int i=0; i<3; i++)   Center[i] = amr->BoxCenter[i];
         }
         break;

         case 3: // potential minimum
         {
            Extrema.Field     = _POTE;
            Extrema.Radius    = __FLT_MAX__;
            Extrema.Center[0] = amr->BoxCenter[0];
            Extrema.Center[1] = amr->BoxCenter[1];
            Extrema.Center[2] = amr->BoxCenter[2];

            Aux_FindExtrema( &Extrema, EXTREMA_MIN, 0, TOP_LEVEL, PATCH_LEAF );

            for (int i=0; i<3; i++)   Center[i] = Extrema.Coord[i];
         }
         break;

         case 4: // CoM
         {
            const double CoM_ref[3]  = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
            const double CoM_MaxR    = __FLT_MAX__;
            const double CoM_MinRho  = 0.0;
            const long   CoM_Field   = _DENS;
            const double CoM_TolErrR = __FLT_MAX__;
            const int    CoM_MaxIter = 1;

            double FinaldR;
            int    FinalNIter;

            Aux_FindWeightedAverageCenter( Center, CoM_ref, CoM_MaxR, CoM_MinRho, CoM_Field, CoM_TolErrR,
                                           CoM_MaxIter, &FinaldR, &FinalNIter );
         }
         break;

         default:
            Aux_Error( ERROR_INFO, "unsupported %s = %d !!\n", "PSC_Prof_Center", PSC_Prof_Center );
      }


//    (2) compute spherically averaged profile
      const bool        RemoveEmpty_Yes = true;
      const double      PrepTime_No     = -1.0;
      const int         NVar            = 2;
      const int         MinLv           = 0;
      const int         MaxLv           = MAX_LEVEL;
      const PatchType_t PatchType  = PATCH_LEAF;

      Profile_t  Dens, Vrad;
      Profile_t *Prof_List[] = { &Dens, &Vrad };
      long       TVar     [] = { _DENS, _VELR };

      Aux_ComputeProfile( Prof_List, Center, PSC_Prof_MaxRadius, PSC_Prof_MinBinSize, PSC_Prof_LogBin,
                          PSC_Prof_LogBinRatio, RemoveEmpty_Yes, TVar, NVar, MinLv, MaxLv, PatchType, PrepTime_No );

//    (3) dump data
      if ( MPI_Rank == 0 )
      {
         char FileName[MAX_STRING];

         sprintf( FileName, "%s/Profile_SphAve_%06ld", OUTPUT_DIR, Step );
         FILE *File = fopen( FileName, "w" );

//       metadata
         Aux_Message( File, "# Step             : %ld\n",                  Step                            );
         Aux_Message( File, "# Time             : %13.7e\n",               Time[0]                         );
         Aux_Message( File, "# Center Method    : %d\n",                   PSC_Prof_Center                 );
         Aux_Message( File, "# Center           : %13.7e %13.7e %13.7e\n", Center[0], Center[1], Center[2] );
         Aux_Message( File, "# Maximum Radius   : %13.7e\n",               Dens.MaxRadius                  );
         Aux_Message( File, "# Minimum Bin Size : %13.7e\n",               PSC_Prof_MinBinSize             );
         Aux_Message( File, "# LogBin           : %d\n",                   PSC_Prof_LogBin                 );
         Aux_Message( File, "# LogBinRatio      : %13.7e\n",               PSC_Prof_LogBinRatio            );
         Aux_Message( File, "# NBin             : %d\n",                   Dens.NBin                       );
         Aux_Message( File, "# ------------------------------------------------------------------------\n" );
         Aux_Message( File, "%5s %9s %22s %22s %22s\n",
                            "# [1]", "[2]", "[3]", "[4]", "[5]" );
         Aux_Message( File, "%5s %9s %22s %22s %22s\n",
                            "# Bin", "NCell", "Bin_Center", "Density", "Vrad" );

//       data
         for (int i=0; i<Dens.NBin; i++)
         fprintf( File, "%5d %9ld %22.14e %22.14e %22.14e\n",
                        i, Dens.NCell[i], Dens.Radius[i], Dens.Data[i] * UNIT_D, Vrad.Data[i] * UNIT_V );

         fclose( File );
      }

   } // if ( PSC_Prof )

} // FUNCTION : Record_PSC



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_PreStellarCore
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_PreStellarCore()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();

// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
   Init_BField_ByVecPot_User_Ptr = NULL; // option: OPT__INIT_BFIELD_BYVECPOT=2;  example: Model_Hydro/MHD_Init_BField_ByVecPot_Function.cpp
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
   Mis_GetTimeStep_User_Ptr      = Mis_GetTimeStep_PSC;
   Aux_Record_User_Ptr           = Record_PSC;
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_PreStellarCore
