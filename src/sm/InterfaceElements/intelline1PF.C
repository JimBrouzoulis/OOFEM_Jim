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
#include "IntElLine1.h"
#include "node.h"
#include "structuralinterfacecrosssection.h"
#include "structuralinterfacematerialstatus.h"
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
    answer.setValues(3, D_u, D_v, T_f);
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
    //answer2.printYourself();
    //this->computeStiffnessMatrix_du(answer3, rMode, tStep); //symmetric
    answer3.beTranspositionOf( answer2 );
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
        double g = this->computeG(ip, VM_Total, tStep);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(N, DN, dA*g);
        } else {
            answer.plusProductUnsym(N, DN, dA*g);
        }
    }


    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
IntElLine1PF :: computeStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
        IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    FloatMatrix N, rotationMatGtoL;
    FloatArray u, traction, tractionTemp, jump, fu, fd(2), fd4(4);

    IntArray IdMask_u;
    this->giveDofManDofIDMask_u( IdMask_u );
    this->computeVectorOfDofIDs(IdMask_u, VM_Total, tStep, u );

    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.resize(8,4);
    answer.zero();

    fd.zero();
    FloatArray a_d_temp, a_d, Bd, Nd;
    this->computeDamageUnknowns( a_d_temp, VM_Total, tStep );
    a_d.setValues(2, a_d_temp.at(1), a_d_temp.at(2) ); 

    FloatMatrix temp, Kdd(8,2);
    fu.zero();
    Kdd.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        this->computeNmatrixAt(ip, N);
        this->computeNd_vectorAt(*ip->giveCoordinates(), Nd);
        jump.beProductOf(N, u);
        this->computeTraction(traction, ip, jump, tStep);

        // compute internal cohesive forces as f = N^T*traction dA
        double Gprim   = computeGPrim(ip, VM_Total, tStep);
        //printf("%e\n",g);
        double dA = this->computeAreaAround(ip);
        
        fu.plusProduct(N, traction, dA*Gprim);
        temp.beDyadicProductOf(fu,Nd);
        Kdd.add(temp);

    }

    IntArray indxu, indx1, indx2;
    indxu.setValues(8, 1,2,3,4,5,6,7,8);
    indx1.setValues(2, 1, 2);
    indx2.setValues(2, 3, 4);
    
    answer.assemble(Kdd, indxu, indx1);
    
    answer.assemble(Kdd, indxu, indx2);

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
    FloatMatrix tempN(2,2), tempB(2,2), temp(2,2);
    FloatArray Nd, Bd;
    tempN.zero();
    tempB.zero();
    temp.zero();
    for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(j);

        computeBd_vectorAt(ip, Bd);
        computeNd_vectorAt(*ip->giveCoordinates(), Nd);
        double dA = this->computeAreaAround(ip);        

        //double Gbis   = this->computeGbis()
        double Gbis = 2.0;
        double Psibar  = this->computeFreeEnergy( ip, tStep );            
            

        // K_dd1 = ( -kp*neg_Mac'(d_dot) / Delta_t + g_c/ l + G'' * Psibar ) * N^t*N; 
        double Delta_d = computeDamageAt(ip, VM_Incremental, tStep);
        //double factorN = -kp * neg_MaCauleyPrime(Delta_d/Delta_t)/Delta_t +  g_c / l + Gbis * Psibar; 
        //double factorN = g_c / l + Gbis * Psibar; 

        double Psibar0 = this->givePsiBar0();
        double factorN = g_c / l + Gbis * this->MaCauley(Psibar-Psibar0);
        

        tempN.beDyadicProductOf(Nd, Nd);
        temp.add(factorN * dA, tempN);
        //tempN.printYourself();

        // K_dd2 =  g_c * l * Bd^t * Bd;
        double factorB = g_c * l;
        tempB.beDyadicProductOf(Bd, Bd);
        temp.add(factorB * dA, tempB);
    }

    IntArray indx1, indx2;
    indx1.setValues(2, 1, 2);
    indx2.setValues(2, 3, 4);
    
    answer.assemble(temp, indx1, indx1);
    
    answer.assemble(temp, indx2, indx2);

    //answer.symmetrized();
    //answer.printYourself();

}




double 
IntElLine1PF :: computeDamageAt(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // d = N_d * a_d
    //NLStructuralElement *el = static_cast< NLStructuralElement* >(this->giveElement( ) );
    StructuralInterfaceElement *el = this->giveElement( );
    FloatArray dVec;
    computeDamageUnknowns(dVec, valueMode, stepN);
    dVec.resizeWithValues(2);
    FloatArray Nvec;
    el->giveInterpolation()->evalN(Nvec, *gp->giveCoordinates(), FEIElementGeometryWrapper(el));
    return Nvec.dotProduct(dVec);
}




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

    FloatMatrix N, rotationMatGtoL;
    FloatArray u, traction, tractionTemp, jump, fu, fd(2), fd4(4);

    IntArray IdMask_u;
    this->giveDofManDofIDMask_u( IdMask_u );
    this->computeVectorOfDofIDs(IdMask_u, VM_Total, tStep, u );

    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.resize(0);
    fd.zero();
    FloatArray a_d_temp, a_d, Bd, Nd;
    this->computeDamageUnknowns( a_d_temp, VM_Total, tStep );
    a_d.setValues(2, a_d_temp.at(1), a_d_temp.at(2) ); 
    if( this->giveNumber() == 273) {
      //  printf("%e %e\n", a_d_temp.at(1), a_d_temp.at(2) );
    }
    fu.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRule->getIntegrationPoint(i);
        this->computeNmatrixAt(ip, N);
        jump.beProductOf(N, u);
        this->computeTraction(traction, ip, jump, tStep);

        // compute internal cohesive forces as f = N^T*traction dA
        double g = this->computeG(ip, VM_Total, tStep);
        //printf("%e\n",g);
        double dA = this->computeAreaAround(ip);
        
        fu.plusProduct(N, traction, dA*g);


        // damage field
        this->computeBd_vectorAt(ip, Bd);
        this->computeNd_vectorAt(*ip->giveCoordinates(), Nd);
        double kp      = this->givePenaltyParameter();
        double Delta_t = tStep->giveTimeIncrement();
        double d       = computeDamageAt( ip, VM_Total, tStep);
        double Delta_d = computeDamageAt(ip, VM_Incremental, tStep);
        double l       = this->giveInternalLength();
        double g_c     = this->giveCriticalEnergy();
        double Gprim   = computeGPrim(ip, VM_Total, tStep);
        double Psibar  = this->computeFreeEnergy( ip, tStep );
        //printf("psi %e \n",Psibar);

        double gradd = Bd.dotProduct(a_d);
        
        //double sectionalForcesScal = -kp * neg_MaCauley(Delta_d/Delta_t) + g_c / l * d + Gprim * Psibar;
        //double sectionalForcesScal =  g_c / l * d + Gprim * Psibar;
        double Psibar0 = this->givePsiBar0();
        double sectionalForcesScal =  g_c / l * d + Gprim * this->MaCauley(Psibar-Psibar0);
   


	    double sectionalForcesVec = g_c * l * gradd;
        fd = fd + ( Nd*sectionalForcesScal + Bd*sectionalForcesVec ) * dA;
     //   
        //double sectionalForcesScal =  -2.0 * Psibar;
        //fd = fd + ( Nd*sectionalForcesScal  ) * dA;
        
    }

    //FloatMatrix Kdd;
    //computeStiffnessMatrix_dd(Kdd, TangentStiffness, tStep); 
    //fd4.beProductOf(Kdd,a_d_temp);

    IntArray indx1, indx2;
    indx1.setValues(2, 1, 2);
    indx2.setValues(2, 3, 4);
    //fd4.zero();
    fd4.assemble(fd, indx1);
    fd4.assemble(fd, indx2);

    //fd4.printYourself();

// total vec
    IntArray IdMask_d;
    this->giveDofManDofIDMask_u( IdMask_u );
    this->giveDofManDofIDMask_d( IdMask_d );
    this->computeLocationArrayOfDofIDs( IdMask_u, loc_u );
    this->computeLocationArrayOfDofIDs( IdMask_d, loc_d );

    int nDofs = this->computeNumberOfDofs();
    answer.resize( nDofs );
    answer.zero();
    answer.assemble(fu, loc_u);
    answer.assemble(fd4, loc_d);
    //answer.printYourself();

}



double
IntElLine1PF :: computeFreeEnergy(GaussPoint *gp, TimeStep *tStep)
{
    
    StructuralInterfaceMaterialStatus *matStat = static_cast< StructuralInterfaceMaterialStatus * >( gp->giveMaterialStatus() );
    FloatArray strain, stress;
    //stress = matStat->giveTempFirstPKTraction();
	//strain = matStat->giveTempJump();
	stress = matStat->giveFirstPKTraction();
	strain = matStat->giveJump();
    //stress.printYourself();
    //strain.printYourself();
    if( this->giveNumber() == 273) {
        //printf("%e \n", 0.5 * stress.dotProduct( strain ) );
    }
    return 0.5 * stress.dotProduct( strain );
}




void
IntElLine1PF :: giveDofManDofIDMask_u(IntArray &answer)
{
	//StructuralInterfaceElement :: giveDofManDofIDMask(-1, EID_MomentumBalance, answer); 
    //IntElLine1 :: giveDofManDofIDMask(-1, EID_MomentumBalance, answer);
    answer.setValues(2, D_u, D_v);
}

void
IntElLine1PF :: giveDofManDofIDMask_d(IntArray &answer)
{
    answer.setValues(1, T_f);
}


void
IntElLine1PF :: computeLocationArrayOfDofIDs( const IntArray &dofIdArray, IntArray &answer )
{
    // Routine to extract compute the location array an element given an dofid array.
    answer.resize( 0 );
    //NLStructuralElement *el = this->giveElement();
    StructuralInterfaceElement *el = this->giveElement();
    int k = 0;
    for(int i = 1; i <= el->giveNumberOfDofManagers(); i++) {
        DofManager *dMan = el->giveDofManager( i );
        for(int j = 1; j <= dofIdArray.giveSize( ); j++) {

            if(dMan->hasDofID( (DofIDItem) dofIdArray.at( j ) )) {
                Dof *d = dMan->giveDofWithID( dofIdArray.at( j ) );
                answer.followedBy( k + d->giveNumber( ) );
            }
        }
        k += dMan->giveNumberOfDofs( );
    }
}

void
IntElLine1PF :: computeNd_vectorAt(const FloatArray &lCoords, FloatArray &N)
{
    StructuralInterfaceElement *el = this->giveElement();
    FloatArray Nvec;
    el->giveInterpolation( )->evalN( N, lCoords, FEIElementGeometryWrapper( el ) );

}

void
IntElLine1PF :: computeBd_vectorAt(GaussPoint *aGaussPoint, FloatArray &answer)
{
    // Returns the [numSpaceDim x nDofs] gradient matrix {B_d} of the receiver,
    // evaluated at gp.

    StructuralInterfaceElement *el = dynamic_cast< StructuralInterfaceElement* > ( this->giveElement() );
    FloatMatrix dNdxi;
    this->giveInterpolation()->evaldNdxi( dNdxi, *aGaussPoint->giveCoordinates( ), FEIElementGeometryWrapper( this ) );

    answer.resize(2);
    for (int i = 1; i <= dNdxi.giveNumberOfRows(); i++) {
        answer.at(i) = dNdxi.at(i,1);
    }
    FloatArray G;
    this->computeCovarBaseVectorAt(aGaussPoint, G);
    answer.times( 1.0 / sqrt(G.dotProduct(G)) );
}

} // end namespace oofem
