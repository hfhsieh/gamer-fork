#include "GAMER.h"

#if ( defined GRAVITY  &&  defined GREP )


Profile_t *DensAve [NLEVEL+1][2];
Profile_t *EngyAve [NLEVEL+1][2];
Profile_t *VrAve   [NLEVEL+1][2];
Profile_t *PresAve [NLEVEL+1][2];
Profile_t *Phi_eff [NLEVEL  ][2];

       int    GREP_LvUpdate;
       int    GREPSg     [NLEVEL];
static double GREPSgTime [NLEVEL][2];

static double Center [3];
static double MaxRadius;
static double MinBinSize;

// temporary switch for test different approaches: do interpolation in ComputeProfile() or in stored profiles
static bool Do_TEMPINT_in_ComputeProfile = true;


static void Init_GREP_Profile();
static void Update_GREP_Profile( const int level, const int Sg, const double PrepTime );
static void Combine_GREP_Profile( Profile_t *Prof[][2], const int level, const int Sg, const double PrepTime,
                                  const bool RemoveEmpty );
extern void SetTempIntPara( const int lv, const int Sg_Current, const double PrepTime, const double Time0, const double Time1,
                            bool &IntTime, int &Sg, int &Sg_IntT, real &Weighting, real &Weighting_IntT );


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_GREffPot
// Description :  Compute the spherical-averaged profile, and GR effective potential.
//                Then set up the CPU/GPU arrays for the GR potential correction
//
// Note        :  1. Invoked by Init_GAMER() and EvolveLevel()
//                2. Enabled by the macros GRAVITY and GREP
//                3. The total averaged profile is stored at QUANT[NLEVEL]
//-------------------------------------------------------------------------------------------------------
void Init_GREffPot( const int level, const double TimeNew )
{

   if ( level == -1 )
   {
//    initialize the Center, MaxRadius, MinBinSize, and GREP profiles at the first call;
      Init_GREP_Profile();

      for (int lv=0; lv<NLEVEL; lv++)
      {
         if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Lv %2d ... ", lv );

         Init_GREffPot( lv, Time[lv] );

         if ( MPI_Rank == 0 )   Aux_Message( stdout, "done\n" );
      }
   }

   else
   {
//    compare the input TimeNew with the stored time to reduce unnecessary calculations
      bool FluUpdateTime;
      int  Sg;

      if      (  Mis_CompareRealValue( TimeNew, GREPSgTime[level][0], NULL, false )  )
      {
         FluUpdateTime = false;
         Sg            = 0;
      }

      else if (  Mis_CompareRealValue( TimeNew, GREPSgTime[level][1], NULL, false )  )
      {
         FluUpdateTime = false;
         Sg            = 1;
      }

      else
      {
         FluUpdateTime = true;
         Sg            = 1 - GREPSg[level];
      }

      FluUpdateTime = true;

      GREP_LvUpdate         = level;
      GREPSg    [level]     = Sg;
      GREPSgTime[level][Sg] = TimeNew;


      if ( FluUpdateTime )
      {
//       update the spherical-averaged profile
         if ( Do_TEMPINT_in_ComputeProfile )
            Update_GREP_Profile( level, Sg, TimeNew );  // do interpolation in ComputeProfile()
         else
            Update_GREP_Profile( level, Sg, -1.0    );  // do interpolation in the stored profiles

//       combine the profile at each level
         Combine_GREP_Profile( DensAve, level, Sg, TimeNew, true );
         Combine_GREP_Profile( EngyAve, level, Sg, TimeNew, true );
         Combine_GREP_Profile( VrAve  , level, Sg, TimeNew, true );
         Combine_GREP_Profile( PresAve, level, Sg, TimeNew, true );

//       compute the GR effective potential
         CPU_ComputeEffPot( DensAve[NLEVEL][Sg], EngyAve[NLEVEL][Sg], VrAve[NLEVEL][Sg], PresAve[NLEVEL][Sg],
                            Phi_eff[level ][Sg] );
      }

//    initialize the auxiliary GPU arrays
#     ifdef GPU
      CUAPI_Init_GREffPot();
#     endif
   } // if ( level == -1 )

} // FUNCTION : Init_GREffPot



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_GREP_Profile
// Description :  Initialize the Profiles_t objects
//-------------------------------------------------------------------------------------------------------
static void Init_GREP_Profile()
{

// initialize static variables
   for (int lv=0; lv<NLEVEL; lv++)
   {
      GREPSg[lv] = 0;
      for (int Sg=0; Sg<2; Sg++)   GREPSgTime[lv][Sg] = -__FLT_MAX__;
   }

// initialize the Center, MaxRadius, and MinBinSize
   switch ( GREP_CENTER_METHOD )
   {
      case 1:
         for (int i=0; i<3; i++)   Center[i] = amr->BoxCenter[i];
      break;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "GREP_CENTER_METHOD", GREP_CENTER_METHOD );
   }

// defaults to the distance between the center and the farthest box vertex
   MaxRadius  = ( GREP_MAXRADIUS > 0.0 )  ? GREP_MAXRADIUS
                                          : SQRT( SQR( MAX( amr->BoxSize[0] - Center[0], Center[0] ) )
                                          +       SQR( MAX( amr->BoxSize[1] - Center[1], Center[1] ) )
                                          +       SQR( MAX( amr->BoxSize[2] - Center[2], Center[2] ) ) );

   MinBinSize = ( GREP_MINBINSIZE > 0.0 ) ? GREP_MINBINSIZE : amr->dh[MAX_LEVEL];

// declare pointer array of GREP profiles
   for (int Sg=0; Sg<2; Sg++)
   for (int lv=0; lv<=NLEVEL; lv++)
   {
      DensAve [lv][Sg] = new Profile_t();
      EngyAve [lv][Sg] = new Profile_t();
      VrAve   [lv][Sg] = new Profile_t();
      PresAve [lv][Sg] = new Profile_t();

      if ( lv < NLEVEL )
      Phi_eff [lv][Sg] = new Profile_t();
   }

} // FUNCTION : Init_GREP_Profile



//-------------------------------------------------------------------------------------------------------
// Function    :  Update_GREP_Profile
// Description :  Update the spherical-averaged profiles
//
// Note        :  1. The contribution from non-leaf patches at the specified level is stored QUANT[NLEVEL]
//-------------------------------------------------------------------------------------------------------
static void Update_GREP_Profile( const int level, const int Sg, const double PrepTime )
{

   long       TVar         [] = {               _DENS,             _VELR,               _PRES,           _EINT_DER };
   Profile_t *Prof_Leaf    [] = { DensAve[ level][Sg], VrAve[ level][Sg], PresAve[ level][Sg], EngyAve[ level][Sg] };
   Profile_t *Prof_NonLeaf [] = { DensAve[NLEVEL][Sg], VrAve[NLEVEL][Sg], PresAve[NLEVEL][Sg], EngyAve[NLEVEL][Sg] };

   if ( Do_TEMPINT_in_ComputeProfile )
//    update the profile from leaf patches at level equal/below the specified level
      Aux_ComputeProfile( Prof_Leaf,    Center, MaxRadius, MinBinSize, GREP_LOGBIN, GREP_LOGBINRATIO, false, TVar,
                          4, -1,    level, PATCH_LEAF,    PrepTime );
   else
   {
//    update the profile from leaf patches at the current level
      Aux_ComputeProfile( Prof_Leaf,    Center, MaxRadius, MinBinSize, GREP_LOGBIN, GREP_LOGBINRATIO, false, TVar,
                          4, level, -1,    PATCH_LEAF,    PrepTime );

//    also update the USG profile to account for the correction from finer level
      if ( ( level < TOP_LEVEL )  &&  ( GREPSgTime[level][Sg_USG] >= 0 ) )
      {
         int               Sg_USG    = 1 - Sg;
         double          Time_USG    = GREPSgTime[level][Sg_USG];
         Profile_t *Prof_Leaf_USG [] = { DensAve[ level][Sg_USG], VrAve[ level][Sg_USG], PresAve[ level][Sg_USG], EngyAve[ level][Sg_USG] };

         Aux_ComputeProfile( Prof_Leaf_USG,    Center, MaxRadius, MinBinSize, GREP_LOGBIN, GREP_LOGBINRATIO, false, TVar,
                             4, level, -1,    PATCH_LEAF,    Time_USG );
      }
   }

// update the profile from the non-leaf patches at the current level
   Aux_ComputeProfile   ( Prof_NonLeaf, Center, MaxRadius, MinBinSize, GREP_LOGBIN, GREP_LOGBINRATIO, false, TVar,
                          4, level, -1,    PATCH_NONLEAF, PrepTime );

} // FUNCTION : Update_GREP_Profile



//-------------------------------------------------------------------------------------------------------
// Function    :  Combine_GREP_Profile
// Description :  Combine the separated spherical-averaged profiles
//
// Note        :  1. The total averaged profile is stored at QUANT[NLEVEL]
//
// Parameter   :  Prof        : Profiles_t object array to be combined
//                RemoveEmpty : true  --> remove empty bins from the data
//                              false --> these empty bins will still be in the profile arrays with
//                                        Data[empty_bin]=Weight[empty_bin]=NCell[empty_bin]=0
//-------------------------------------------------------------------------------------------------------
void Combine_GREP_Profile( Profile_t *Prof[][2], const int level, const int Sg, const double PrepTime,
                           const bool RemoveEmpty )
{

   Profile_t *Prof_NonLeaf = Prof[NLEVEL][Sg];

// combine contribution from leaf patches equal/below the specified level, and non-leaf patches at lv=level
   if ( Do_TEMPINT_in_ComputeProfile )
   {
      Profile_t *Prof_Leaf = Prof[level][Sg];

//    temporal interpolation is done in Aux_ComputeProfile()
      for (int b=0; b<Prof_Leaf->NBin; b++)
      {
         if ( Prof_Leaf->NCell[b] == 0L )  continue;

         Prof_NonLeaf->Data  [b]  = Prof_NonLeaf->Weight[b] * Prof_NonLeaf->Data[b]
                                  + Prof_Leaf   ->Weight[b] * Prof_Leaf   ->Data[b];
         Prof_NonLeaf->Weight[b] += Prof_Leaf   ->Weight[b];
         Prof_NonLeaf->NCell [b] += Prof_Leaf   ->NCell [b];

         Prof_NonLeaf->Data  [b] /= Prof_NonLeaf->Weight[b];
      }
   }

   else
   {
      for (int lv=0; lv<=level; lv++)
      {
//       temporal interpolation parameters
         bool FluIntTime;
         int  FluSg, FluSg_IntT;
         real FluWeighting, FluWeighting_IntT;

         SetTempIntPara( lv, GREPSg[lv], PrepTime, GREPSgTime[lv][0], GREPSgTime[lv][1],
                         FluIntTime, FluSg, FluSg_IntT, FluWeighting, FluWeighting_IntT );

         Profile_t *Prof_Leaf      = Prof[lv][FluSg];
         Profile_t *Prof_Leaf_IntT = ( FluIntTime ) ? Prof[lv][FluSg_IntT] : NULL;


//CHECK: if the radius of each bin is the same in Prof_Leaf and Prof_Leaf_IntT?
//       add contributions from leaf patches at level=lv
         for (int b=0; b<Prof_Leaf->NBin; b++)
         {
            if ( Prof_Leaf->NCell[b] == 0L )  continue;

            Prof_NonLeaf->Data  [b]  = ( FluIntTime )
                                     ?                       Prof_NonLeaf  ->Weight[b] * Prof_NonLeaf  ->Data[b]
                                       + FluWeighting      * Prof_Leaf     ->Weight[b] * Prof_Leaf     ->Data[b]
                                       + FluWeighting_IntT * Prof_Leaf_IntT->Weight[b] * Prof_Leaf_IntT->Data[b]
                                     :                       Prof_NonLeaf  ->Weight[b] * Prof_NonLeaf  ->Data[b]
                                       +                     Prof_Leaf     ->Weight[b] * Prof_Leaf     ->Data[b];

            Prof_NonLeaf->Weight[b] += ( FluIntTime )
                                     ?   FluWeighting      * Prof_Leaf     ->Weight[b]
                                       + FluWeighting_IntT * Prof_Leaf_IntT->Weight[b]
                                     :                       Prof_Leaf     ->Weight[b];

            Prof_NonLeaf->NCell [b] += Prof_Leaf   ->NCell [b];
            Prof_NonLeaf->Data  [b] /= Prof_NonLeaf->Weight[b];
         } // for (int b=0; b<Prof_Leaf->NBin; b++)
      } // for (int lv=0; lv<=level, lv++)
   } // if ( Do_TEMPINT_in_ComputeProfile )


// remove the empty bins in Prof_NonLeaf
   if ( RemoveEmpty )
   {
      for (int b=0; b<Prof_NonLeaf->NBin; b++)
      {
         if ( Prof_NonLeaf->NCell[b] != 0L )   continue;

   //    for cases of consecutive empty bins
         int b_up;
         for (b_up=b+1; b_up<Prof_NonLeaf->NBin; b_up++)
            if ( Prof_NonLeaf->NCell[b_up] != 0L )   break;

         const int stride = b_up - b;

         for (int b_up=b+stride; b_up<Prof_NonLeaf->NBin; b_up++)
         {
            const int b_up_ms = b_up - stride;

            Prof_NonLeaf->Radius[b_up_ms] = Prof_NonLeaf->Radius[b_up];
            Prof_NonLeaf->Data  [b_up_ms] = Prof_NonLeaf->Data  [b_up];
            Prof_NonLeaf->Weight[b_up_ms] = Prof_NonLeaf->Weight[b_up];
            Prof_NonLeaf->NCell [b_up_ms] = Prof_NonLeaf->NCell [b_up];
         }

   //    reset the total number of bins
         Prof_NonLeaf->NBin -= stride;
      } // for (int b=0; b<Prof_NonLeaf->NBin; b++)

//    update the maximum radius since the last bin may have not been removed
      const int LastBin = Prof_NonLeaf->NBin-1;

      Prof_NonLeaf->MaxRadius = ( Prof_NonLeaf->LogBin ) ? Prof_NonLeaf->Radius[LastBin] * sqrt( Prof_NonLeaf->LogBinRatio )
                                                         : Prof_NonLeaf->Radius[LastBin] + 0.5*MinBinSize;
   }

} // FUNCTION : Combine_GREP_Profile



#endif // #if ( defined GRAVITY  &&  defined GREP )
