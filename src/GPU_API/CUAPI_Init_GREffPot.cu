#include "CUAPI.h"
#include "Profile.h"

#if ( defined GPU  &&  defined GRAVITY  &&  defined GREP )


#include "CUPOT.h"

extern int        GREP_LvUpdate;
extern int        GREPSg  [NLEVEL];
extern Profile_t *Phi_eff [NLEVEL][2];


// declare the GPU kernel requiring GREP_Data, GREP_EdgeL, GREP_Center, and r_max2
int CUPOT_SetConstMem_GREffPot( double h_GREP_Data[], double h_GREP_Radius[], double h_GREP_Center[],
                                int    h_GREP_NBin );



//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Init_GREffPot
// Description :  Set the auxiliary GPU constant-memory arrays for the GR effective potential
//
// Note        :  1. Invoked by Init_GREffPot()
//-------------------------------------------------------------------------------------------------------
void CUAPI_Init_GREffPot()
{

   Profile_t *Phi = Phi_eff[GREP_LvUpdate][ GREPSg[GREP_LvUpdate] ];

// check
   if ( Phi->NBin > GR_POT_NAUX_MAX )
      Aux_Error( ERROR_INFO, "Too many bins in average radial Profile %d !!\n", Phi->NBin );


   int Exitcode = CUPOT_SetConstMem_GREffPot( Phi->Data, Phi->Radius, Phi->Center, Phi->NBin );
   if (  Exitcode != 0  )
      Aux_Error( ERROR_INFO, "CUPOT_SetConstMem_GREffPot failed... Exitcode %d...\n", Exitcode );

} // FUNCTION : CUAPI_Init_GREffPot



#endif // #if ( defined GPU  &&  defined GRAVITY  &&  defined GREP )
