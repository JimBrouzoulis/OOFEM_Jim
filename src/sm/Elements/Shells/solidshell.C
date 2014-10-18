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

#include "Elements/Shells/solidshell.h"
#include <Materials/structuralms.h>
#include "classfactory.h"
#include "fei3dhexalin.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "crosssection.h"

namespace oofem {
REGISTER_Element(SolidShell);

FEI3dHexaLin SolidShell :: interpolation;

SolidShell :: SolidShell(int n, Domain *aDomain) : LSpace(n, aDomain)
{
    numberOfDofMans  = 8;
    numberOfGaussPoints = 8;
    
}


void
SolidShell :: postInitialize() {
    
    LSpace :: postInitialize();
 
    
    int numEASparam = 0;
    
    switch ( this->EAS_type) { 
      case 1:   
        numEASparam = 1;
        break;
        
      case 2:
        numEASparam = 3;
        break;
    }
    
    this->u_k.resize(numberOfDofMans*3);
    this->u_k.zero();
    this->fE.resize(numEASparam);
    this->fE.zero();
    this->KEC.resize(numEASparam, numberOfDofMans*3);
    this->KEC.zero();
    this->invKEE.resize(numEASparam, numEASparam);
    this->invKEE.zero();
    this->alpha.resize(numEASparam);
    this->alpha.zero();
    
}

void
SolidShell :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 6);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

FEInterpolation *SolidShell :: giveInterpolation() const { return & interpolation; }


IRResultType
SolidShell :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 8;
    IRResultType result = this->NLStructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    
    // Check if EAS should be used
    this->EAS_type = 0;
    if ( ir->hasField(_IFT_SolidShell_EAS_type) ) {
      
        IR_GIVE_FIELD(ir, this->EAS_type,   _IFT_SolidShell_EAS_type);
      
    }
    

    return IRRT_OK;
}



void
SolidShell :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [ 6 x (8*3) ] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    FloatArray lCoords = * gp->giveNaturalCoordinates();
    interp->evaldNdx( dNdx, lCoords, FEIElementGeometryWrapper(this) );
    
    answer.resize(6, dNdx.giveNumberOfRows() * 3);
    answer.zero();
#if  0
    
    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);
        answer.at(3, 3 * i - 0) = dNdx.at(i, 3);

        answer.at(4, 3 * i - 1) = dNdx.at(i, 3);
        answer.at(4, 3 * i - 0) = dNdx.at(i, 2);
        answer.at(5, 3 * i - 2) = dNdx.at(i, 3);
        answer.at(5, 3 * i - 0) = dNdx.at(i, 1);
        answer.at(6, 3 * i - 2) = dNdx.at(i, 2);
        answer.at(6, 3 * i - 1) = dNdx.at(i, 1);
    }
    
#else
    
    // Do ANS
    // Need to evaluate strain at four points
    FloatArray A = { 0.0, -1.0, 0.0};
    FloatArray B = { 1.0,  0.0, 0.0};
    FloatArray C = { 0.0,  1.0, 0.0};
    FloatArray D = {-1.0,  0.0, 0.0};
    
    FloatMatrix dNdxA, dNdxB, dNdxC, dNdxD;
    interp->evaldNdx( dNdxA, A, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxB, B, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxC, C, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxD, D, FEIElementGeometryWrapper(this) );
    
    FloatMatrix dNdx0; 
    interp->evaldNdx( dNdx0, {lCoords.at(1), lCoords.at(2), 0.0}, FEIElementGeometryWrapper(this) );
        
    
    double NA = 0.5 * ( 1.0 - lCoords.at(2) );
    double NC = 0.5 * ( 1.0 + lCoords.at(2) );
    double NB = 0.5 * ( 1.0 - lCoords.at(1) );
    double ND = 0.5 * ( 1.0 + lCoords.at(1) );
    
    // E_xz = E(5) = N_A(gp) * B(A)(5,:) + N_C(gp) * B(C)(5,:)
    // E_yz = E(4) = N_B(gp) * B(B)(4,:) + N_D(gp) * B(D)(4,:)

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);
        answer.at(6, 3 * i - 2) = dNdx.at(i, 2);
        answer.at(6, 3 * i - 1) = dNdx.at(i, 1);
        
        // Evaluate thickness strain at the mid-surface
        //answer.at(3, 3 * i - 0) = dNdx.at(i, 3);
        answer.at(3, 3 * i - 0) = dNdx0.at(i, 3);
            
        
//         answer.at(4, 3 * i - 1) = NB * dNdxA.at(i, 3) + ND * dNdxC.at(i, 3);
//         answer.at(4, 3 * i - 0) = NB * dNdxA.at(i, 2) + ND * dNdxC.at(i, 2);

        answer.at(4, 3 * i - 1) = NB * dNdxB.at(i, 3) + ND * dNdxD.at(i, 3);
        answer.at(4, 3 * i - 0) = NB * dNdxB.at(i, 2) + ND * dNdxD.at(i, 2);

        
        answer.at(5, 3 * i - 2) = NA * dNdxA.at(i, 3) + NC * dNdxC.at(i, 3);
        answer.at(5, 3 * i - 0) = NA * dNdxA.at(i, 1) + NC * dNdxC.at(i, 1);
    }
#endif
}


void
SolidShell :: computeEASBmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
  
    FloatArray lCoords = *gp->giveNaturalCoordinates();
    double xi   = lCoords.at(1);
    double eta  = lCoords.at(2);
    double zeta = lCoords.at(3);
    FloatMatrix M;
    
     switch ( this->EAS_type) { 
       case 1: 
      
        // alt 1
        M.resize(6,1);
        M.zero();
        M.at(3,1) = zeta;
        break;
    
      case 2:
        // alt 2
        M.resize(6,3);
        M.zero();
        M.at(1,1) = xi;
        M.at(2,2) = eta;
        M.at(3,3) = zeta;
        break;
    }
    
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix J0mat; 
    FloatArray center = {0.0, 0.0, 0.0};

    interp->giveJacobianMatrixAt(J0mat, center, FEIElementGeometryWrapper(this) );
//     double J0 = interp->giveTransformationJacobian(center, FEIElementGeometryWrapper(this) );
    double J0 = J0mat.giveDeterminant();
    double J  = interp->giveTransformationJacobian(lCoords, FEIElementGeometryWrapper(this) );
  
    FloatMatrix T(6,6);
    this->computeBondTransformationMatrix(T, J0mat);
    //T.printYourself("T");
    
    //T.beUnitMatrix();
    answer.clear();
    answer.plusProductUnsym(T, M, J/J0);
    
    
}


 
void
SolidShell :: computeBondTransformationMatrix(FloatMatrix &answer, FloatMatrix &base)
{
  // given a bases {g1, g2, g3} compute the Bond (Voigt form) transformation matrix.
  // TODO is this with OOFEM ordering if components?
  FloatArray x, y, z;
  x.beColumnOf(base,1);
  y.beColumnOf(base,2);
  z.beColumnOf(base,3);
  answer = {
          { x(0) * x(0), x(1) * x(1), x(2) * x(2), x(1) * x(2), x(0) * x(2), x(0) * x(1) },
          { y(0) * y(0), y(1) * y(1), y(2) * y(2), y(1) * y(2), y(0) * y(2), y(0) * y(1) },
          { z(0) * z(0), z(1) * z(1), z(2) * z(2), z(1) * z(2), z(0) * z(2), z(0) * z(1) },
          { 2 * y(0) * z(0), 2 * y(1) * z(1), 2 * y(2) * z(2), y(2) * z(1) + y(1) * z(2), y(2) * z(0) + y(0) * z(2), y(1) * z(0) + y(0) * z(1) },
          { 2 * x(0) * z(0), 2 * x(1) * z(1), 2 * x(2) * z(2), x(2) * z(1) + x(1) * z(2), x(2) * z(0) + x(0) * z(2), x(1) * z(0) + x(0) * z(1) },
          { 2 * x(0) * y(0), 2 * x(1) * y(1), 2 * x(2) * y(2), x(2) * y(1) + x(1) * y(2), x(2) * y(0) + x(0) * y(2), x(1) * y(0) + x(0) * y(1) }
        };
  
}

void
SolidShell :: computeAlpha(FloatArray &answer, FloatArray &u)
{
    // compute alpha based on displacement update
    FloatArray deltaU;
    deltaU.beDifferenceOf(u,this->u_k);
   
    FloatMatrix KEE_inv = this->invKEE;
    FloatArray fE = this->fE;
    FloatMatrix KEC = this->KEC;
    FloatArray temp, deltaAlpha;
    temp.beProductOf(KEC, deltaU);
    temp.add(fE);
    deltaAlpha.beProductOf(KEE_inv,temp);
    
    FloatArray oldAlpha = this->alpha; // last converged values
    answer = this->alpha - deltaAlpha;
    
    
    // set current u-displcement
    this->u_k = u;
}

void
SolidShell :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    if ( this->EAS_type ) {
        FloatArray fC, fE, strain, u, vStrainC, vStrainE, vStress;
        FloatMatrix KEE, KBE, BC, BE, D;
        fE.clear();
        fC.clear();
        
        KEE.clear();
        
        this->computeVectorOf(VM_Total, tStep, u);
        FloatArray alpha;
        this->computeAlpha(alpha, u );
        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {

            computeBmatrixAt(gp, BC);
            this->computeEASBmatrixAt(gp, BE);
            
            vStrainC.beProductOf(BC, u);
            vStrainE.beProductOf(BE, alpha);
            vStrainC.add(vStrainE);
            this->computeStressVector(vStress, vStrainC, gp, tStep);
            double dV  = this->computeVolumeAround(gp);

            
            // Compute nodal internal forces at nodes as f = B^T*Stress dV
            fC.plusProduct(BC, vStress, dV);
            fE.plusProduct(BE, vStress, dV);
            
        }
        this->fE = fE;
        
        FloatMatrix KEE_inv, KEC;
        KEE_inv = this->invKEE;
        KEC = this->KEC;
        
        FloatArray temp, f;
        temp.beProductOf(KEE_inv,fE);
        answer = fC;
        answer.plusProduct(KEC, temp, -1.0);
        
    } else {
      
        LSpace :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
      
    }
    
    
}


void
SolidShell :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
  if ( this->EAS_type ) {
        FloatArray fC, fE, strain, u, vStrainC, vStrainE, vStress;
        FloatMatrix KEE, KEC, KCC, BC, BE, D, DBC, DBE;
        
        KEC.clear();
        KEE.clear();
        KCC.clear();
        
        
        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
            //LSpace :: computeBmatrixAt(gp, BC);
            this->computeBmatrixAt(gp, BC);
            this->computeEASBmatrixAt(gp, BE);
            double dV  = this->computeVolumeAround(gp);

            this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
            DBC.beProductOf(D, BC);
            DBE.beProductOf(D, BE);
            
            KCC.plusProductUnsym(BC, DBC, dV);
            KEC.plusProductUnsym(BE, DBC, dV);
            KEE.plusProductUnsym(BE, DBE, dV);
            
        }

        
        FloatMatrix KEE_inv;
        KEE_inv.beInverseOf(KEE);
        this->invKEE = KEE_inv;
        this->KEC = KEC;
        
        answer= KCC;
        
        FloatMatrix K, tempmat;
        tempmat.beProductOf(KEE_inv, KEC);
        
        answer.plusProductUnsym(KEC, tempmat, -1.0);
      
    } else {
      
        LSpace :: computeStiffnessMatrix(answer, rMode, tStep);
      
    }
}

  
  
  

} // end namespace oofem
