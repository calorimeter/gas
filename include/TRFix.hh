/*
 * TRFix.hh
 *
 *  Created on: Feb 9, 2016
 *      Author: pepijn
 */


#ifndef TRFIX_H
#define TRFIX_H

#include "G4ForwardXrayTR.hh"

class TRFix : public G4ForwardXrayTR
{
  public:
	TRFix(  const G4String& processName="XrayTR"     );
};

#endif /* TRFIX_HH_ */
