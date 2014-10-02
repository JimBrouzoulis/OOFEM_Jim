/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2014   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "Contact/celnode2edge.h"
#include "floatmatrix.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "Contact/contactpair.h"


namespace oofem {
  
  
  
Node2EdgeContact :: Node2EdgeContact(int num, Domain *d, ContactDefinition *cDef, ContactPair *cPair) : StructuralContactElement(num, d, cDef, cPair)
{   
    this->numberOfDofMans = 3;
    this->cPair = static_cast< ContactPairNode2Edge * > ( cPair );
};   
 


void
Node2EdgeContact :: computeBmatrixAt(const FloatArray &lCoords, const FloatArray &traction, FloatMatrix &answer, TimeStep *tStep)
{
    // rename
    this->cPair->computeBmatrixAt(lCoords, traction, answer, tStep);
  
}


void
Node2EdgeContact :: setupIntegrationPoints()
{
    // Sets up the integration rule array which contains all the necessary integration points
    if ( this->integrationRule == NULL ) {
        //TODO sets a null pointer for the element in the iRule 
        this->integrationRule = new GaussIntegrationRule(1, this) ;
        this->integrationRule->SetUpPointsOnLine(1, _Unknown);
    }
    
  
}


ContactPair* 
Node2EdgeContact :: giveContactPair() 
{ 
    return dynamic_cast< ContactPair* > ( this->cPair ); 
  
}





    
}












