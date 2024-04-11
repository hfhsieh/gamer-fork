#include "GAMER.h"


extern bool   CCSN_CC_MaxRefine_Flag1;
extern bool   CCSN_CC_MaxRefine_Flag2;
extern int    CCSN_CC_MaxRefine_LV1;
extern int    CCSN_CC_MaxRefine_LV2;
extern double CCSN_CC_MaxRefine_Dens1;
extern double CCSN_CC_MaxRefine_Dens2;
extern double CCSN_CentralDens;

extern double CCSN_AngRes_Min;
extern double CCSN_AngRes_Max;



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_CoreCollapse
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
bool Flag_CoreCollapse( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   bool Flag      = false;
   bool MaxRefine = false;

   const double dh        = amr->dh[lv];
   const double Pos   [3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   const double Center[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   const double dR[3]     = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const double R         = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );

   const double CentralDens = CCSN_CentralDens / UNIT_D;


// (1) check if the allowed maximum level is reached
   if      ( CCSN_CC_MaxRefine_Flag1  &&  CentralDens < CCSN_CC_MaxRefine_Dens1 / UNIT_D )
      MaxRefine = lv >= CCSN_CC_MaxRefine_LV1;

   else if ( CCSN_CC_MaxRefine_Flag2  &&  CentralDens < CCSN_CC_MaxRefine_Dens2 / UNIT_D )
      MaxRefine = lv >= CCSN_CC_MaxRefine_LV2;


   if ( !MaxRefine ) {
//    (2) check if the minimum angular resolution is reached
      if ( CCSN_AngRes_Min > 0.0  &&  R * CCSN_AngRes_Min < dh )
         Flag = true;

//    (3) always refined to highest level in the region with r < 30 km
//    (3-a) always refine the innermost cells
      if ( R < amr->dh[lv] )
         Flag = true;

//    (3-b) refine the region with r < 30 km
      if ( R * UNIT_L < 3e6 )
         Flag = true;
   }


   return Flag;

} // FUNCTION : Flag_CoreCollapse



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Lightbulb
// Description :  Check if the element (i,j,k) of the input data satisfies the user-defined flag criteria
//
// Note        :  1. Invoked by "Flag_Check" using the function pointer "Flag_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this function will become useless
//                2. Enabled by the runtime option "OPT__FLAG_USER"
//                3. For lightbulb test problem
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
bool Flag_Lightbulb( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   bool Flag = false;

   const double dh        = amr->dh[lv];
   const double Pos   [3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   const double Center[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   const double dR[3]     = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const double R         = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );

   const real (*Rho )[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS];


// TODO: fine-tune the criteria
// (1) always refined to highest level in the region with r < 30 km
   if ( R * UNIT_L < 3e6 )
      Flag = true;

// (2) check if the minimum angular resolution is reached
   if ( !Flag  &&  CCSN_AngRes_Min > 0.0  &&  R * CCSN_AngRes_Min < dh )
      Flag = true;


   return Flag;

} // FUNCTION : Flag_Lightbulb



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Region_CCSN
// Description :  Check if the element (i,j,k) of the input patch is within the regions allowed to be refined
//                for the CCSN test problem
//
// Note        :  1. Invoked by Flag_Check() using the function pointer "Flag_Region_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__FLAG_REGION"
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[0][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
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


   const double Center[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   const double dR[3]     = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const double R         = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );


// check allowed maximum refine level based on angular resolution
   if (  CCSN_AngRes_Max > 0.0  &&  2.0 * R * CCSN_AngRes_Max >  dh  &&
      !( CCSN_AngRes_Min > 0.0  &&        R * CCSN_AngRes_Min >= dh )  )
      Within = False;


   return Within;

} // FUNCTION : Flag_Region_CCSN
