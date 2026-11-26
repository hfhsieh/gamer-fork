#include "GAMER.h"

#if ( MODEL == HYDRO )


extern double CCSN_PeakDens_CGS;
extern double CCSN_Rsh_Ave;

extern bool   CCSN_Ref_CC_Flag1;
extern bool   CCSN_Ref_CC_Flag2;
extern int    CCSN_Ref_CC_MaxLv1;
extern int    CCSN_Ref_CC_MaxLv2;
extern double CCSN_Ref_CC_Dens1;
extern double CCSN_Ref_CC_Dens2;
extern double CCSN_Ref_MaxLv_Rad_Code;
extern double CCSN_Ref_BaseRad_Code;
extern double CCSN_Ref_AngRes_Min_Code;
extern double CCSN_Ref_AngRes_Max_Code;
extern double CCSN_Ref_Shock_AngResMax_Code;
extern double CCSN_Ref_Shock_Rsh_Code;
extern double CCSN_Ref_Shock_Rad_Code;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_CoreCollapse
// Description :  Check if the element (i,j,k) of the input data satisfies the user-defined flag criteria
//
// Note        :  1. Invoked by "Flag_Check" using the function pointer "Flag_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this function will become useless
//                2. Enabled by the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k     : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the target patch
//                PID       : ID of the target patch
//                Threshold : User-provided threshold for the flag operation, which is loaded from the
//                            file "Input__Flag_User"
//                            In order of radius_min, radius_max, threshold_dens
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_CoreCollapse( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   bool Flag     = false;
   bool Is_MaxLv = false;

   const double dh     = amr->dh[lv];
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
#  ifdef GREP
   const double dR [3] = { Pos[0]-GREP_Center[0],    Pos[1]-GREP_Center[1],    Pos[2]-GREP_Center[2]    };
#  else
   const double dR [3] = { Pos[0]-amr->BoxCenter[0], Pos[1]-amr->BoxCenter[1], Pos[2]-amr->BoxCenter[2] };
#  endif
   const double R      = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );

   const double CentralDens = CCSN_PeakDens_CGS;


// (1) check if the allowed maximum level is reached
   if      ( CCSN_Ref_CC_Flag1  &&  CentralDens < CCSN_Ref_CC_Dens1 )
      Is_MaxLv = lv >= CCSN_Ref_CC_MaxLv1;

   else if ( CCSN_Ref_CC_Flag2  &&  CentralDens < CCSN_Ref_CC_Dens2 )
      Is_MaxLv = lv >= CCSN_Ref_CC_MaxLv2;


   if ( !Is_MaxLv )
   {
//    (2) check if the minimum angular resolution is reached
      if ( CCSN_Ref_AngRes_Min_Code > 0.0  &&  R * CCSN_Ref_AngRes_Min_Code < dh )
         Flag = true;

//    (3) always refine the innermost cells
      if ( R < amr->dh[lv] )
         Flag = true;

//    (4) always refine to the allowed maximum level for cells with r <= CCSN_Ref_MaxLv_Rad
      if ( R <= CCSN_Ref_MaxLv_Rad_Code )
         Flag = true;
   }


   return Flag;

} // FUNCTION : Flag_CoreCollapse



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_PostBounce
// Description :  Check if the element (i,j,k) of the input data satisfies the user-defined flag criteria
//
// Note        :  1. Invoked by "Flag_Check" using the function pointer "Flag_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this function will become useless
//                2. Enabled by the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k     : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the target patch
//                PID       : ID of the target patch
//                Threshold : User-provided threshold for the flag operation, which is loaded from the
//                            file "Input__Flag_User"
//                            In order of radius_min, radius_max, threshold_dens
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_PostBounce( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   bool Flag = false;

   const double dh     = amr->dh[lv];
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
#  ifdef GREP
   const double dR [3] = { Pos[0]-GREP_Center[0],    Pos[1]-GREP_Center[1],    Pos[2]-GREP_Center[2]    };
#  else
   const double dR [3] = { Pos[0]-amr->BoxCenter[0], Pos[1]-amr->BoxCenter[1], Pos[2]-amr->BoxCenter[2] };
#  endif
   const double R = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );


// (1) always refine the innermost cells
   if ( R < amr->dh[lv] )
      Flag = true;

// (2) always refine to the maximum level for cells with r <= CCSN_Ref_MaxLv_Rad
   if ( R <= CCSN_Ref_MaxLv_Rad_Code )
      Flag = true;

// (3) check if the minimum angular resolution is reached
   if ( CCSN_Ref_AngRes_Min_Code > 0.0  &&  R * CCSN_Ref_AngRes_Min_Code < dh )
      Flag = true;


   return Flag;

} // FUNCTION : Flag_PostBounce



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Region_CCSN
// Description :  Check if the element (i,j,k) of the input patch is within the regions allowed to be refined
//                for the CCSN test problem
//
// Note        :  1. Invoked by Flag_Check() using the function pointer "Flag_Region_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__FLAG_REGION"
//
// Parameter   :  i,j,k : Indices of the target element in the patch ptr[0][lv][PID]
//                lv    : Refinement level of the target patch
//                PID   : ID of the target patch
//
// Return      :  "true/false"  if the input cell "is/is not" within the region allowed for refinement
//-------------------------------------------------------------------------------------------------------
bool Flag_Region_CCSN( const int i, const int j, const int k, const int lv, const int PID )
{

   const double dh     = amr->dh[lv];
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   bool Within = true;


#  ifdef GREP
   const double dR [3] = { Pos[0]-GREP_Center[0],    Pos[1]-GREP_Center[1],    Pos[2]-GREP_Center[2]    };
#  else
   const double dR [3] = { Pos[0]-amr->BoxCenter[0], Pos[1]-amr->BoxCenter[1], Pos[2]-amr->BoxCenter[2] };
#  endif
   const double R      = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );


// (1) always allow refinement to the maximum level for cells with r <= CCSN_Ref_MaxLv_Rad
   if ( R <= CCSN_Ref_MaxLv_Rad_Code )
      return true;

// (2) check the maximum allowed refinement level based on angular resolution
   if ( CCSN_Ref_AngRes_Max_Code > 0.0  &&  2.0 * R * CCSN_Ref_AngRes_Max_Code > dh )
      Within = false;

   if ( !Within )   return Within;

// (3) check the maximum allowed refinement level based on the shock-radius-dependent angular resolution
   if ( CCSN_Ref_Shock_AngResMax_Code           >  0.0           &&
        CCSN_Ref_Shock_Rsh_Code                 <= CCSN_Rsh_Ave  &&
        CCSN_Ref_Shock_Rad_Code                 <= R             &&
        2.0 * R * CCSN_Ref_Shock_AngResMax_Code >  dh               )
      Within = false;

   if ( !Within )   return Within;

// (4) check allowed maximum refinement level based on distance to the center
   if ( CCSN_Ref_BaseRad_Code > 0.0 )
   {
      int dlv   = 0;
      int ratio = (int) ( R / CCSN_Ref_BaseRad_Code );

      while ( ratio )   { dlv += 1; ratio >>= 1; }

      if ( lv + dlv >= MAX_LEVEL )  Within = false;
   }


   return Within;

} // FUNCTION : Flag_Region_CCSN



#endif // #if ( MODEL == HYDRO )
