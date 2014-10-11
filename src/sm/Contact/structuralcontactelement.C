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

#include "Contact/structuralcontactelement.h"


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
  
  
  
StructuralContactElement :: StructuralContactElement(int num, Domain *d, ContactDefinition *cDef, ContactPair *cPair) : ContactElement(num, d, cDef)
{   
    this->cPair = cPair;
};   
  
int
StructuralContactElement :: instanciateYourself(DataReader *dr)
{
    return 1;
}


void
StructuralContactElement :: computeGap(FloatArray &answer, FloatArray &lCoords, TimeStep *tStep)
{
    this->giveContactPair()->computeGap(answer, lCoords, tStep);
    // Rotate gap to a local system
    FloatMatrix orthoBase;
    this->giveContactPair()->computeCurrentTransformationMatrixAt(answer, orthoBase, tStep);
    answer.rotatedWith(orthoBase, 'n');
    if ( answer.at(3) <= 0.0 ) {
        this->setContactFlag();
        printf("    YES    ");
    } else { 
        printf("    NO     "); 
      
    }
}

void
StructuralContactElement :: performCPP(GaussPoint *gp, TimeStep *tStep) 
{
    this->cPair->performCPP(gp, tStep);
}


void
StructuralContactElement :: computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep)
{   
    StructuralInterfaceMaterial *mat = static_cast < StructuralInterfaceMaterial *> (this->giveContactMaterial() );
    mat->giveEngTraction_3d(t, gp, gap, tStep);
}





void
StructuralContactElement :: computeContactForces(FloatArray &answer, TimeStep *tStep)
{
    
    FloatArray gap;
    int ndofs = this->giveNumberOfDispDofs();
    answer.resize(ndofs);
    answer.zero();
    
    for ( GaussPoint *gp : *this->integrationRule ) {
        
        this->giveContactPair()->performCPP(gp, tStep);
        FloatArray lCoords = *gp->giveNaturalCoordinates();
        
        if ( !lCoords.giveSize() ) {
            this->resetContactFlag();
            continue;  
        }
        this->computeGap(gap, lCoords, tStep); 
        if ( this->isInContact() ) {
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
StructuralContactElement :: computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep)
{
    int ndofs = this->giveNumberOfDispDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    
    for ( GaussPoint *gp : *this->integrationRule ) {
  
        FloatArray lCoords = *gp->giveNaturalCoordinates();
        if ( this->isInContact() ) {
            // D( delta(g) * t ) = delta(g) * D(t) + t * Delta( delta(g) ) = 1 + 2
          
            // Part 1: delta(g) * D(t) = delta(g) * dt/dg * D(g) 
            FloatMatrix N, D;
            
            this->computeNmatrixAt( lCoords, N); 
            
            StructuralInterfaceMaterial *mat = static_cast < StructuralInterfaceMaterial* > ( this->giveContactMaterial() );
            mat->give3dStiffnessMatrix_Eng(D, TangentStiffness, gp, tStep); 
            
            FloatMatrix globalSys, DN;
            this->computeCurrentTransformationMatrixAt( lCoords, globalSys, tStep );
            D.rotatedWith(globalSys, 'n');                   // transform to global system
            DN.beProductOf(D,N);
            double dA = this->computeCurrentAreaAround(gp, tStep);
            answer.plusProductUnsym(N, DN, dA);  
            
            
            // Addition of non-symmetric and highly nonlinear term
            // part 2: t * dN/dxi * dxi/dx * Delta(x)
            // TODO doesn't give any effect on the convergence - wrong or unimportant!?
            //FloatArray t;
            //StructuralInterfaceMaterialStatus *matStat = static_cast < StructuralInterfaceMaterialStatus* > ( mat->giveStatus(gp) );
            //t = matStat->giveTempTraction();
            //FloatMatrix B;
            //computeBmatrixAt(lCoords, t, B, tStep); //TODO rename and recheck
            //answer.add(dA, B);
        }
    }
}


  
void
StructuralContactElement :: computeNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer)
{
    this->giveContactPair()->computeNmatrixAt(lCoords, answer);
}



void
StructuralContactElement :: giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s)
{
    // should return location array for the master and the slave node
  
    ContactPair *cPair = this->giveContactPair();
   
    int numDofs = ( cPair->giveNumberOfMasterNodes() + cPair->giveNumberOfSlaveNodes() ) * 3;
    answer.resize(numDofs);
    answer.zero();
    IntArray dofIdArray = {D_u, D_v, D_w};

    // First add eq numbers for the slaves
    int pos = 0;
    for ( int i = 1; i <= cPair->giveNumberOfSlaveNodes(); i++ ) {
        Node *node = cPair->giveSlaveNode(i); 
        for ( int i = 1; i <= dofIdArray.giveSize(); i++ ) {
            DofIDItem id = (DofIDItem) dofIdArray.at(i);
            if ( node->hasDofID( id ) ) { // add corresponding number
                Dof *dof= node->giveDofWithID( id );
                answer.at(pos + i) = s.giveDofEquationNumber( dof );
            } 
        }    
      
        pos += 3;
    }

    // Add eq. for master nodes
    for ( int i = 1; i <= cPair->giveNumberOfMasterNodes(); i++ ) {
        Node *node = cPair->giveMasterNode(i);
        for ( int i = 1; i <= dofIdArray.giveSize(); i++ ) {
            DofIDItem id = (DofIDItem) dofIdArray.at(i);
            if ( node->hasDofID( id ) ) { // add corresponding number
                Dof *dof= node->giveDofWithID( id );
                answer.at(pos + i) = s.giveDofEquationNumber( dof );
            } 
        }    
      
        pos += 3;
    }
   
}    


double
StructuralContactElement :: computeCurrentAreaAround(IntegrationPoint *ip, TimeStep *tStep)
{
    return this->giveContactPair()->computeCurrentAreaAround(ip, tStep);
}


void
StructuralContactElement :: computeCurrentTransformationMatrixAt(const FloatArray &lCoords, FloatMatrix &answer, TimeStep *tStep)
{
    this->cPair->computeCurrentTransformationMatrixAt(lCoords, answer, tStep);
}


int
StructuralContactElement :: giveNumberOfDispDofs() 
{ 
    return 3 * (this->giveContactPair()->giveNumberOfMasterNodes() + this->giveContactPair()->giveNumberOfSlaveNodes() );
}







//---------------------------------------
// Lagrange multiplier associated methods
//---------------------------------------


StructuralContactElementLagrange :: StructuralContactElementLagrange(int num, Domain *d, ContactDefinition *cDef, ContactPair *cPair) : 
       StructuralContactElement(num, d, cDef, cPair)
{

    
         
}

void
StructuralContactElementLagrange :: computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep)
{
    // Evaluates the contact traction vector, the gap should be in local system
    // TODO but g_N->0 => t_N = 0 => t_T = 0
    StructuralInterfaceMaterial *mat = static_cast < StructuralInterfaceMaterial *> (this->giveContactMaterial() );
    mat->giveEngTraction_3d(t, gp, gap, tStep);
          
    
    // The normal contact pressure is determined from the Lagarange multipliers 
    Dof *dof = this->giveDofManager(1)->giveDofWithID( this->giveDofIdArray().at(1) );
    double lambda = dof->giveUnknown(VM_Total, tStep);
    t.at(3) = lambda;

} 


void
StructuralContactElementLagrange :: giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s)
{
    // should return location array for the master and the slave node plus the lagrange multipliers
    StructuralContactElement ::giveLocationArray(answer, s);
    
    
    for ( int i = 1; i <= this->giveContactPair()->giveNumberOfSlaveNodes(); i++ ) {
      
        Node *node = this->giveContactPair()->giveSlaveNode(i);
        if ( node->hasDofID( (DofIDItem)this->giveDofIdArray().at(1) ) ) { 
      
            Dof *dof= node->giveDofWithID( (DofIDItem)this->giveDofIdArray().at(1) );
            answer.followedBy( s.giveDofEquationNumber(dof) );
        }
    }

   
}    




void
StructuralContactElementLagrange :: computeContactForces(FloatArray &answer, TimeStep *tStep)
{
    
    FloatArray gap, temp;
    int numDispDofs = this->giveNumberOfDispDofs();
    int ndofs = numDispDofs + this->integrationRule[0].giveNumberOfIntegrationPoints();
    
    answer.resize(ndofs);
    answer.zero();
    
    for ( GaussPoint *gp : *this->integrationRule ) {
        this->giveContactPair()->performCPP(gp, tStep);
        FloatArray lCoords = *gp->giveNaturalCoordinates();
        if ( !lCoords.giveSize() ) {
            this->resetContactFlag();
            continue;  
        }
        
        this->computeGap(gap, lCoords, tStep); // local system
        if ( this->isInContact() ) {
            FloatArray t;
            this->computeContactTractionAt(gp, t, gap, tStep); // local system
            
            FloatMatrix N, globalSys;
            this->computeNmatrixAt(lCoords, N);
            this->computeCurrentTransformationMatrixAt( lCoords, globalSys, tStep );
            t.rotatedWith(globalSys, 't');                     // transform to global system
            double dA = this->computeCurrentAreaAround(gp, tStep);
            temp.plusProduct(N, t, dA);
            
            answer.at( numDispDofs + gp->giveNumber() ) = gap.at(3)*dA; // The gap constraints are scaled with dA to give a symmetric tangent.
            
        }
    }    
    
    if ( temp.giveSize() ) {
        answer.addSubVector(temp, 1);
    }

    
}
  
  
  
void
StructuralContactElementLagrange :: computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep)
{

    
    int numDispDofs = this->giveNumberOfDispDofs();
    int ndofs = numDispDofs + this->integrationRule[0].giveNumberOfIntegrationPoints();
        
    answer.resize(ndofs, ndofs);
    answer.zero();

    FloatMatrix temp;
    FloatArray normal, C;
    FloatMatrix globalSys, DN;
    temp.clear();
    for ( GaussPoint *gp : *this->integrationRule ) {
  
         FloatArray lCoords = *gp->giveNaturalCoordinates();

         if ( this->isInContact() ) {
            // D( delta(g) * t ) = delta(g) * D(t) + t * Delta( delta(g) ) = 1 + 2
          
            // Part 1: delta(g) * D(t) = delta(g) * dt/dg * D(g) 
            FloatMatrix N, D;
            this->computeNmatrixAt( lCoords, N); 
            
            StructuralInterfaceMaterial *mat = static_cast < StructuralInterfaceMaterial* > ( this->giveContactMaterial() );
            mat->give3dStiffnessMatrix_Eng(D, TangentStiffness, gp, tStep); 
            D.at(3,3) = 0.0;
            
            
            this->computeCurrentTransformationMatrixAt( lCoords, globalSys, tStep );
            D.rotatedWith(globalSys, 'n');                   // transform to global system
            DN.beProductOf(D,N);
            double dA = this->computeCurrentAreaAround(gp, tStep);
            temp.plusProductUnsym(N, DN, dA);  
            
            
            this->giveContactPair()->computeNmatrixAt(lCoords, N);
            this->giveContactPair()->computeCurrentNormalAt(lCoords, normal, tStep);
            C.beTProductOf(N, normal);
            C.times(dA); // gives contact forces
            answer.addSubVectorRow(C, numDispDofs + gp->giveNumber(), 1 );
            answer.addSubVectorCol(C, 1, numDispDofs + gp->giveNumber() );
        }
    }
    
    if ( temp.giveNumberOfColumns() ) {
        IntArray loc;
        loc.enumerate(numDispDofs);
        answer.assemble(temp, loc, loc);
    }

    //TODO need to add a small number for the solver
    for ( int i = 1; i <= ndofs; i++ ) {
        answer.at(i,i) += 1.0e-9;
    }
}  
  

  
  
void
StructuralContactElementLagrange :: giveDofManagersToAppendTo(IntArray &answer)
{
    int numSlaveNodes = this->giveContactPair()->giveNumberOfSlaveNodes();
    answer.resize(numSlaveNodes);
    
    for ( int i = 1; i <= numSlaveNodes; i++ ) {
      Node *node = this->giveContactPair()->giveSlaveNode(i);
      answer.at(i) = node->giveNumber();  
    }
  
    
}
    
  
    
}












