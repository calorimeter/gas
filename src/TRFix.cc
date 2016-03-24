/*
 * TRFix.cc
 *
 *  Created on: Feb 9, 2016
 *      Author: pepijn
 */

#include "TRFix.hh"

TRFix::
TRFix( const G4String& processName  )
   :        G4ForwardXrayTR(processName)
{
	fMatIndex1 = 1;
	fMatIndex2 = 2;
}


