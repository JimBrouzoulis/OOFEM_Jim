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
#include "contact/contactdefinition.h"
#include "feinterpol.h"
#include "Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "Contact/contactpair.h"

// temporary
#include "crosssection.h"

namespace oofem {
  
  
  
Node2EdgeContact :: Node2EdgeContact(int num, Domain *d, ContactDefinition *cDef, ContactPair *cPair) : ContactElement(num, d, cDef)
{   
    this->numberOfDofMans = 3;
    this->cPair = static_cast< ContactPairNode2Edge * > ( cPair );
};   
  
int
Node2EdgeContact :: instanciateYourself(DataReader *dr)
{
    return 1;
}


void
Node2EdgeContact :: computeGap(FloatArray &answer, FloatArray &lCoords, TimeStep *tStep)
{
    this->cPair->computeGap(answer, lCoords, tStep);
    
    // Rotate gap to a local system
    FloatMatrix orthoBase;
    computeCurrentTransformationMatrixAt( lCoords, orthoBase, tStep );
    answer.rotatedWith(orthoBase, 't');

}

void
Node2EdgeContact :: performCPP(GaussPoint *gp, TimeStep *tStep) 
{
    this->cPair->performCPP(tStep);
    gp->setNaturalCoordinates( { this->cPair->giveCPPcoord() } );
}


void
Node2EdgeContact :: computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep)
{   
    StructuralInterfaceMaterial *mat = static_cast < StructuralInterfaceMaterial *> (this->giveContactMaterial() );
    mat->giveEngTraction_3d(t, gp, gap, tStep);
}





void
Node2EdgeContact :: computeContactForces(FloatArray &answer, TimeStep *tStep)
{
    answer.clear();
    FloatArray gap;
    
    for ( GaussPoint *gp : *this->integrationRule ) {
        this->performCPP(gp, tStep);
        FloatArray lCoords = *gp->giveNaturalCoordinates();
        
        this->computeGap(gap, lCoords, tStep); 
        if ( gap.at(3) < 0.0 ) {
            ;
            FloatArray t;
            this->computeContactTractionAt(gp, t, gap, tStep); // local system
            
            FloatMatrix N, globalSys;
            this->computeNmatrixAt(lCoords, N);
            this->computeCurrentTransformationMatrixAt( lCoords, globalSys, tStep );
            t.rotatedWith(globalSys, 't');                     // transform to global system
            double dA = this->computeCurrentAreaAround(gp, tStep);
            answer.plusProduct(N, t, dA);
        }
    }
}

  
void
Node2EdgeContact :: computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep)
{
    int ndofs = this->giveNumberOfDofManagers() * 3;
    answer.resize(ndofs, ndofs);
    answer.zero();
    
    for ( GaussPoint *gp : *this->integrationRule ) {
  
        this->performCPP(gp, tStep);
        FloatArray lCoords = *gp->giveNaturalCoordinates();
        FloatArray gap;
        this->computeGap(gap, lCoords, tStep);
    
        if( gap.at(3) < 0.0 ) {
            // D( delta(g) * t ) = delta(g) * D(t) + t * Delta( delta(g) ) = 1 + 2
          
            // Part 1: delta(g) * D(t) = delta(g) * dt/dg * D(g) 
            FloatMatrix N, D;
            
            this->computeNmatrixAt( lCoords, N); 
            
            StructuralInterfaceMaterial *mat = static_cast < StructuralInterfaceMaterial* > ( this->giveContactMaterial() );
            mat->give3dStiffnessMatrix_Eng(D, TangentStiffness, gp, tStep); 
            
            FloatMatrix globalSys, DN;
            this->computeCurrentTransformationMatrixAt( lCoords, globalSys, tStep );
            D.rotatedWith(globalSys, 't');                   // transform to global system
            DN.beProductOf(D,N);
            double dA = this->computeCurrentAreaAround(gp, tStep);
            answer.plusProductUnsym(N, DN, dA);  
            
            
            // Addition of non-symmetric and highly nonlinear term
            // part 2: t * dN/dxi * dxi/dx * Delta(x)
            // TODO doesn't give any effect on the convergence - wrong or unimportant!?
            FloatArray t;
            StructuralInterfaceMaterialStatus *matStat = static_cast < StructuralInterfaceMaterialStatus* > ( mat->giveStatus(gp) );
            t = matStat->giveTempTraction();
            FloatMatrix B;
            computeBmatrixAt(lCoords, t, B, tStep); //TODO rename and recheck
            //answer.add(dA, B);
        }
    }
}
  
  
void
Node2EdgeContact :: computeNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer)
{
    this->cPair->computeNmatrixAt(lCoords, answer);
}

void
Node2EdgeContact :: computeBmatrixAt(const FloatArray &lCoords, const FloatArray &traction, FloatMatrix &answer, TimeStep *tStep)
{
  
    this->cPair->computeBmatrixAt(lCoords, traction, answer, tStep);
  
}


void
Node2EdgeContact :: giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s)
{
    // should return location array for the master and the slave node
    // TODO this whole thing should be rewritten  
  
  
    answer.resize(9);
    answer.zero();
    IntArray dofIdArray = {D_u, D_v, D_w};

    // slave node
    for ( int i = 1; i <= dofIdArray.giveSize(); i++ ) {
        if ( this->giveDofManager(1)->hasDofID( (DofIDItem)dofIdArray.at(i) ) ) { // add corresponding number
            Dof *dof= this->giveDofManager(1)->giveDofWithID( (DofIDItem)dofIdArray.at(i) );
            answer.at(i) = s.giveDofEquationNumber( dof );
        } 
    }    
    
    // master node 1
    for ( int i = 1; i <= dofIdArray.giveSize(); i++ ) {
        if ( this->giveDofManager(2)->hasDofID( (DofIDItem)dofIdArray.at(i) ) ) { // add corresponding number
            Dof *dof= this->giveDofManager(2)->giveDofWithID( (DofIDItem)dofIdArray.at(i) );
            answer.at(3+i) = s.giveDofEquationNumber( dof );
        } 
    }    
    
    // master node 2
    for ( int i = 1; i <= dofIdArray.giveSize(); i++ ) {
        if ( this->giveDofManager(3)->hasDofID( (DofIDItem)dofIdArray.at(i) ) ) { // add corresponding number
            Dof *dof= this->giveDofManager(3)->giveDofWithID( (DofIDItem)dofIdArray.at(i) );
            answer.at(6+i) = s.giveDofEquationNumber( dof );
        } 
    }    
    
    //Element :: giveLocationArray(answer, s); will only give eqns for the dof ids I have not D_w for example. TODO 
}    


void
Node2EdgeContact :: setupIntegrationPoints()
{
    // Sets up the integration rule array which contains all the necessary integration points
    if ( this->integrationRule == NULL ) {
        //TODO sets a null pointer for the element in the iRule 
        this->integrationRule = new GaussIntegrationRule(1, NULL) ;
        this->integrationRule->SetUpPointsOnLine(1, _Unknown);
    }
    
  
}



double
Node2EdgeContact :: computeCurrentAreaAround(IntegrationPoint *ip, TimeStep *tStep)
{
    FloatArray g;
    this->computeCovarBaseVectorAt( *ip->giveNaturalCoordinates(), g, tStep);
    double weight  = ip->giveWeight();
    double ds = sqrt( g.dotProduct(g) ) * weight;
    double thickness  = this->cPair->giveMasterElement()->giveCrossSection()->give(CS_Thickness, ip);

    return ds * thickness;
    
}


void
Node2EdgeContact :: computeCovarBaseVectorAt(const FloatArray &lCoords, FloatArray &g, TimeStep *tStep)
{   
  this->cPair->computeCovarBaseVectorAt(lCoords, g, tStep); 
}




void
Node2EdgeContact :: computeCurrentNormalAt(const FloatArray &lCoords, FloatArray &normal, TimeStep *tStep)
{   
    FloatArray g;
    this->cPair->computeCovarBaseVectorAt(lCoords, g, tStep);
    normal.beVectorProductOf(g, {0, 0, 1});
    normal.normalize();
  
}




void
Node2EdgeContact :: computeCurrentTransformationMatrixAt(const FloatArray &lCoords, FloatMatrix &answer, TimeStep *tStep)
{
    // Transformation matrix to the local coordinate system

    FloatArray normal;
    this->computeCurrentNormalAt(lCoords, normal, tStep);
    answer.beLocalCoordSys( normal );
    
}





    
}












