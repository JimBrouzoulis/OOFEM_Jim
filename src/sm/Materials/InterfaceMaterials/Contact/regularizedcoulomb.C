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

#include "Materials/InterfaceMaterials/Contact/regularizedcoulomb.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(RegularizedCoulomb);

RegularizedCoulomb :: RegularizedCoulomb(int n, Domain *d) : StructuralInterfaceMaterial(n, d) { }

RegularizedCoulomb :: ~RegularizedCoulomb() { }


void
RegularizedCoulomb :: giveEngTraction_3d( FloatArray &answer, GaussPoint *gp, const FloatArray &gap, TimeStep *tStep)
{
    // Returns the (engineering) traction vector in 3d based on the spatial gap.

    RegularizedCoulombStatus *status = static_cast< RegularizedCoulombStatus * >( this->giveStatus( gp ) );
    this->initTempStatus(gp);
    
    computeEngTraction_3d( answer, gp, gap, tStep);
    
    // Update gp    
    status->letTempJumpBe( gap );
    status->letTempTractionBe( answer );

}


void
RegularizedCoulomb :: computeEngTraction_3d( FloatArray &answer, GaussPoint *gp, const FloatArray &gap, TimeStep *tStep)
{
    // Returns the (engineering) traction vector in 3d based on the spatial gap.

    RegularizedCoulombStatus *status = static_cast< RegularizedCoulombStatus * >( this->giveStatus( gp ) );
    this->initTempStatus(gp);
    
    // Traction vector
    answer.resize( 3 );
    answer.zero();
    double gN = gap.at( 3 );
    
    if ( gN < 0.0 ) {
        double pN = this->epsN * gN;
    
        
        FloatArray nT; 
        nT = { gap.at(1), gap.at(2) };
        double nTnorm = nT.computeNorm();
        if ( nTnorm < 1.0e-10 ) {
            nT.zero();
        } else {
            nT.times(1.0 / nTnorm);
        }
        
        FloatArray oldGap = status->giveJump(); // jump = gap
        FloatArray dgT; 
        dgT.beDifferenceOf( {gap.at(1), gap.at(2) }, { oldGap.at(1), oldGap.at(2) } );
        double gTdot = dgT.computeNorm() / tStep->giveTimeIncrement();
        double phi = this->computePhi(gTdot);
        
        answer.at( 1 ) = this->mu * abs(pN) * phi * nT.at( 1 );
        answer.at( 2 ) = this->mu * abs(pN) * phi * nT.at( 2 );
        answer.at( 3 ) = pN;
    }
 
}



double 
RegularizedCoulomb :: computePhi(const double gTdot)
{
 
    return tanh( gTdot / this->epsReg );
  
}  
  

void
RegularizedCoulomb :: give3dStiffnessMatrix_Eng( FloatMatrix &answer, MatResponseMode rMode,
                                            GaussPoint *gp, TimeStep *tStep )
{ 
    RegularizedCoulombStatus *status = static_cast< RegularizedCoulombStatus * > ( this->giveStatus(gp) );
    FloatArray gap0, t0;
    gap0 = status->giveTempJump();

    this->computeEngTraction_3d( t0, gp, gap0, tStep);
    
    const double gN = gap0.at( 3 );
    answer.resize(3,3);
    answer.zero();
    if( rMode == TangentStiffness ) {
        if( gN <= 0) {
            const double eps = 1.0e-8;
            FloatArray gap, t, dt;
            for ( int i = 1; i <= 3; i++ ) {
                gap = gap0;
                gap.at(i) += eps;
                computeEngTraction_3d( t, gp, gap, tStep);
                dt.beDifferenceOf(t,t0);
                answer.addSubVectorCol(dt, 1, i);
            }
            answer.times( 1.0 / eps );
        }
    } else {
        OOFEM_ERROR("only supports TangentStiffness");
    }
    
}





IRResultType
RegularizedCoulomb :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->mu, _IFT_RegularizedCoulomb_mu);

    // Optional
    this->epsN = 1.0e6;
    this->epsT = 1.0e6;
    this->epsReg = 1.0e-6;    
    IR_GIVE_OPTIONAL_FIELD(ir, this->epsN, _IFT_RegularizedCoulomb_epsN);
    IR_GIVE_OPTIONAL_FIELD(ir, this->epsT, _IFT_RegularizedCoulomb_epsT);
    IR_GIVE_OPTIONAL_FIELD(ir, this->epsReg, _IFT_RegularizedCoulomb_epsReg);
    
    return StructuralInterfaceMaterial :: initializeFrom( ir );
}


void
RegularizedCoulomb :: giveInputRecord(DynamicInputRecord &input)
{
     StructuralInterfaceMaterial :: giveInputRecord(input);
     input.setField(this->mu, _IFT_RegularizedCoulomb_mu);
     input.setField(this->epsN, _IFT_RegularizedCoulomb_epsN);
     input.setField(this->epsT, _IFT_RegularizedCoulomb_epsT);
     input.setField(this->epsReg, _IFT_RegularizedCoulomb_epsReg);
}


RegularizedCoulombStatus :: RegularizedCoulombStatus(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g)
{
}


RegularizedCoulombStatus :: ~RegularizedCoulombStatus()
{ }


void
RegularizedCoulombStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
     StructuralInterfaceMaterialStatus::printOutputAt( file, tStep );
}


void
RegularizedCoulombStatus :: initTempStatus()
{
    StructuralInterfaceMaterialStatus::initTempStatus( );
}




} // end namespace oofem
