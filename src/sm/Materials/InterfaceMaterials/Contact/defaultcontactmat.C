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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#include "Materials/InterfaceMaterials/Contact/defaultcontactmat.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(DefaultContactMat);

DefaultContactMat :: DefaultContactMat(int n, Domain *d) : StructuralInterfaceMaterial(n, d) 
{ 
    this->epsN = 1.0e6;
}


void
DefaultContactMat :: giveEngTraction_3d( FloatArray &answer, GaussPoint *gp, const FloatArray &gap, TimeStep *tStep)
{
    // Returns the (engineering) traction vector in 3d based on the spatial gap.

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    this->initTempStatus(gp);
    
    answer.resize( 3 );
    answer.zero();
    double gN = gap.at( 3 );
    
    if ( gN < 0.0 ) {
        answer.at(3) = this->epsN * gN;
    }
    status->letTempJumpBe( gap );
    status->letTempTractionBe( answer );

}

  

void
DefaultContactMat :: give3dStiffnessMatrix_Eng( FloatMatrix &answer, MatResponseMode rMode,
                                            GaussPoint *gp, TimeStep *tStep )
{ 
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    FloatArray gap;
    gap = status->giveTempJump();
    answer.resize(3,3);
    answer.zero();
    if( gap.at(3) <= 0.0 ) {
        answer.at(3,3) =  this->epsN;
    }
    
}



IRResultType
DefaultContactMat :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->epsN = 1.0e6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->epsN, _IFT_DefaultContactMat_epsN);
    
    return StructuralInterfaceMaterial :: initializeFrom( ir );
}


void
DefaultContactMat :: giveInputRecord(DynamicInputRecord &input)
{
     StructuralInterfaceMaterial :: giveInputRecord(input);
     input.setField(this->epsN, _IFT_DefaultContactMat_epsN);
}







} // end namespace oofem
