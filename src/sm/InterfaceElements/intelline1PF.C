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

#include "IntElLine1PF.h"
#include "node.h"
#include "structuralinterfacecrosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei2dlinelin.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(IntElLine1PF);

FEI2dLineLin IntElLine1PF :: interp(1, 1);


IntElLine1PF :: IntElLine1PF(int n, Domain *aDomain) :
    StructuralInterfaceElement(n, aDomain), PhaseFieldElement(n, aDomain)
{
    numberOfDofMans = 4;
}


void
IntElLine1PF :: computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{
    // Returns the modified N-matrix which multiplied with u give the spatial jump.

    FloatArray N;
    FEInterpolation *interp = this->giveInterpolation();
    interp->evalN( N, * ip->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 8);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = -N.at(1);
    answer.at(1, 3) = answer.at(2, 4) = -N.at(2);

    answer.at(1, 5) = answer.at(2, 6) = N.at(1);
    answer.at(1, 7) = answer.at(2, 8) = N.at(2);
}


void
IntElLine1PF :: computeGaussPoints()
{
    // Sets up the array of Gauss Points of the receiver.
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        //integrationRulesArray[0] = new LobattoIntegrationRule (1,domain, 1, 2);
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(4, _2dInterface); ///@todo - should be a parameter with num of ip
    }
}

void
IntElLine1PF :: computeCovarBaseVectorAt(IntegrationPoint *ip, FloatArray &G)
{
    FloatMatrix dNdxi;
    FEInterpolation *interp = this->giveInterpolation();
    interp->evaldNdxi( dNdxi, * ip->giveCoordinates(), FEIElementGeometryWrapper(this) );
    G.resize(2);
    G.zero();
    int numNodes = this->giveNumberOfNodes();
    for ( int i = 1; i <= dNdxi.giveNumberOfRows(); i++ ) {
        double X1_i = 0.5 * ( this->giveNode(i)->giveCoordinate(1) + this->giveNode(i + numNodes / 2)->giveCoordinate(1) ); // (mean) point on the fictious mid surface
        double X2_i = 0.5 * ( this->giveNode(i)->giveCoordinate(2) + this->giveNode(i + numNodes / 2)->giveCoordinate(2) );
        G.at(1) += dNdxi.at(i, 1) * X1_i;
        G.at(2) += dNdxi.at(i, 1) * X2_i;
    }
}

double
IntElLine1PF :: computeAreaAround(IntegrationPoint *ip)
{
    FloatArray G;
    this->computeCovarBaseVectorAt(ip, G);

    double weight  = ip->giveWeight();
    double ds = sqrt( G.dotProduct(G) ) * weight;

    double thickness  = this->giveCrossSection()->give(CS_Thickness, ip);
    return ds * thickness;
}


IRResultType
IntElLine1PF :: initializeFrom(InputRecord *ir)
{
    return StructuralInterfaceElement :: initializeFrom(ir);
}


void
IntElLine1PF :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(2, D_u, D_v);
}

void
IntElLine1PF :: computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // Transformation matrix to the local coordinate system
    FloatArray G;
    this->computeCovarBaseVectorAt(gp, G);
    G.normalize();

    answer.resize(2, 2);
    answer.at(1, 1) =  G.at(1);
    answer.at(2, 1) =  G.at(2);
    answer.at(1, 2) = -G.at(2);
    answer.at(2, 2) =  G.at(1);

}

FEInterpolation *
IntElLine1PF :: giveInterpolation() const
{
    return & interp;
}




void
IntElLine1PF :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    ///@todo this part is enough to do once
    IntArray IdMask_u, IdMask_d;
    this->giveDofManDofIDMask_u( IdMask_u );
    this->giveDofManDofIDMask_d( IdMask_d );
    this->computeLocationArrayOfDofIDs( IdMask_u, loc_u );
    this->computeLocationArrayOfDofIDs( IdMask_d, loc_d );

    int nDofs = this->computeNumberOfDofs();
    answer.resize( nDofs, nDofs );
    answer.zero();

    FloatMatrix answer1, answer2, answer3, answer4;
    this->computeStiffnessMatrix_uu(answer1, rMode, tStep);
    //this->computeStiffnessMatrix_ud(answer2, rMode, tStep);
    //this->computeStiffnessMatrix_du(answer3, rMode, tStep); //symmetric
    //answer3.beTranspositionOf( answer2 );
    this->computeStiffnessMatrix_dd(answer4, rMode, tStep);
    
    answer.assemble( answer1, loc_u, loc_u );
    //answer.assemble( answer2, loc_u, loc_d );
    //answer.assemble( answer3, loc_d, loc_u );
    answer.assemble( answer4, loc_d, loc_d );
}


void
IntElLine1PF :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // Computes the stiffness matrix of the receiver K_cohesive = int_A ( N^t * dT/dj * N ) dA
    FloatMatrix N, D, DN;
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.resize(0, 0);

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    FloatMatrix rotationMatGtoL;
    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(j);

        if ( this->nlGeometry == 0 ) {
            this->giveStiffnessMatrix_Eng(D, rMode, ip, tStep);
        } else if ( this->nlGeometry == 1 ) {
            this->giveStiffnessMatrix_dTdj(D, rMode, ip, tStep);
        } else {
            OOFEM_ERROR("nlgeometry must be 0 or 1!")
        }

        this->computeTransformationMatrixAt(ip, rotationMatGtoL);
        D.rotatedWith(rotationMatGtoL, 't');                      // transform stiffness to global coord system

        this->computeNmatrixAt(ip, N);
        DN.beProductOf(D, N);
        double dA = this->computeAreaAround(ip);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(N, DN, dA);
        } else {
            answer.plusProductUnsym(N, DN, dA);
        }
    }


    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}




void
IntElLine1PF :: computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) 
{    
    // Computation of tangent: K_dd = \int Nd^t * ( -kp*neg_Mac'(alpha_dot)/delta_t + g_c/l + G''*Psibar) * Nd + 
    //                                \int Bd^t * (  g_c * l * [G^1 G^2]^t * [G^1 G^2] ) * Bd 
    //                              = K_dd1 + K_dd2
    int ndofs   = 8 + 4;
	int ndofs_u = 8;
	int ndofs_d = 4;

    double l       = this->giveInternalLength();
    double g_c     = this->giveCriticalEnergy();
    double kp      = this->givePenaltyParameter();
    double Delta_t = tStep->giveTimeIncrement();

    answer.resize( ndofs_d, ndofs_d );
    answer.zero();


    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    FloatMatrix Nd, Bd, temp(2,2);

    temp.zero();
    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(j);

        computeBd_matrixAt(ip, Bd);
        computeNd_matrixAt(*ip->giveCoordinates(), Nd);
        double dA = this->computeAreaAround(ip);        

        //double Gbis   = this->computeGbis()
        double Gbis = 2.0;
        double Psibar  = this->computeFreeEnergy( ip, tStep );            
            

        // K_dd1 = ( -kp*neg_Mac'(d_dot) / Delta_t + g_c/ l + G'' * Psibar ) * N^t*N; 
        double Delta_d = computeDamageAt(ip, VM_Incremental, tStep);
        //double factorN = -kp * neg_MaCauleyPrime(Delta_d/Delta_t)/Delta_t +  g_c / l + Gbis * Psibar; 
        double factorN = g_c / l + Gbis * Psibar; 
        temp.plusProductSymmUpper(Nd, Nd, factorN * dA);
            

        // K_dd2 =  g_c * l * Bd^t * Bd;
        double factorB = g_c * l;
        temp.plusProductSymmUpper(Bd, Bd, factorB * dA);

    }

    IntArray indx1, indx2;
    indx1.setValues(2, 1, 2);
    indx2.setValues(2, 3, 4);
    answer.assemble(temp, indx1, indx1);
    answer.assemble(temp, indx2, indx2);

    answer.symmetrized();
    

}




double 
IntElLine1PF :: computeDamageAt(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // d = N_d * a_d
    NLStructuralElement *el = static_cast< NLStructuralElement* >(this->giveElement( ) );
    FloatArray dVec;
    computeDamageUnknowns(dVec, valueMode, stepN);
    dVec.resizeWithValues(2);
    FloatArray Nvec;
    el->giveInterpolation()->evalN(Nvec, *gp->giveCoordinates(), FEIElementGeometryWrapper(el));
    return Nvec.dotProduct(dVec);
}


//void
//IntElLine1PF :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
//{
//    // Computes the internal forces corresponding to the two fields u & d
//    IntArray IdMask_u, IdMask_d;
//    this->giveDofManDofIDMask_u( IdMask_u );
//    this->giveDofManDofIDMask_d( IdMask_d );
//    this->computeLocationArrayOfDofIDs( IdMask_u, loc_u );
//    this->computeLocationArrayOfDofIDs( IdMask_d, loc_d );
//    
//    int ndofs = this->computeNumberOfDofs();
//    answer.resize( ndofs);
//    answer.zero();
//
//    FloatArray answer_u(0);
//    FloatArray answer_d(0);
//
//    this->giveInternalForcesVector_u(answer_u, tStep, useUpdatedGpRecord);
//    this->giveInternalForcesVector_d(answer_d, tStep, useUpdatedGpRecord);
//    answer.assemble(answer_u, loc_u);
//    answer.assemble(answer_d, loc_d);
//}




void
IntElLine1PF :: giveInternalForcesVector(FloatArray &answer,
                                                       TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes internal forces
    // if useGpRecord == 1 then data stored in ip->giveStressVector() are used
    // instead computing stressVector through this->ComputeStressVector();
    // this must be done after you want internal forces after element->updateYourself()
    // has been called for the same time step.

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    FloatMatrix N, rotationMatGtoL, Nd, Bd;
    FloatArray u, traction, tractionTemp, jump, fu, fd;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);

    FloatArray a_d;
    this->computeDamageUnknowns( a_d, VM_Total, tStep );

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        this->computeNmatrixAt(ip, N);
        jump.beProductOf(N, u);
        this->computeTraction(traction, ip, jump, tStep);

        // compute internal cohesive forces as f = N^T*traction dA
        double dA = this->computeAreaAround(ip);
        
        fu.plusProduct(N, traction, dA);


        // damage field
        computeBd_matrixAt(ip, Bd);
        computeNd_matrixAt(*ip->giveCoordinates(), Nd);
        double kp      = this->givePenaltyParameter();
        double Delta_t = tStep->giveTimeIncrement();
        double d       = computeDamageAt( ip, VM_Total, tStep);
        double Delta_d = computeDamageAt(ip, VM_Incremental, tStep);
        double l       = this->giveInternalLength();
        double g_c     = this->giveCriticalEnergy();
        double Gprim   = computeGPrim(ip, VM_Total, tStep);
        double Psibar  = this->computeFreeEnergy( ip, tStep );
        
        // Dalpha/DX = Dalpha/Dxi * dxi/DX = B*a * 
        // Dalpha/Ds = Dalpha/Dxi * dxi/ds = B*a *  
        FloatArray temp;

        temp.beProductOf(Bd, a_d); 
        double gradd = temp.at(1);
        //double sectionalForcesScal = -kp * neg_MaCauley(Delta_d/Delta_t) + g_c / l * d + Gprim * Psibar;
        double sectionalForcesScal =  g_c / l * d + Gprim * Psibar;
   
	    double sectionalForcesVec = g_c * l * gradd;
        fd.plusProduct(Nd, sectionalForcesScal, dA);
        fd.plusProduct(Bd, sectionalForcesVec, dA);

    }

}



void
IntElLine1PF :: giveDofManDofIDMask_u(IntArray &answer)
{
	StructuralInterfaceElement :: giveDofManDofIDMask(-1, EID_MomentumBalance, answer); 
}

void
IntElLine1PF :: giveDofManDofIDMask_d(IntArray &answer)
{
    answer.setValues(1, T_f);
}




} // end namespace oofem
