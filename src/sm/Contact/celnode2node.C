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

#include "Contact/celnode2node.h"
#include "floatmatrix.h"
#include "masterdof.h"
//#include "unknownnumberingscheme.h"

#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "contact/contactdefinition.h"
#include "Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
namespace oofem {
  
  
  
Node2NodeContact :: Node2NodeContact(int num, Domain *d, ContactDefinition *cDef) : ContactElement(num, d, cDef)
{   
    this->area = 1.0;   // should be optional parameter
    this->epsN = 1.0e6; // penalty - should be given by 'contactmaterial'
    this->epsT = 1.0e6; 
    this->numberOfDofMans = 2;
};   
  
int
Node2NodeContact :: instanciateYourself(DataReader *dr)
{
    // compute normal as direction vector from master node to slave node
    FloatArray xs, xm, normal;
    xs = *this->giveDofManager(1)->giveCoordinates();
    xm = *this->giveDofManager(2)->giveCoordinates();
    normal = xs-xm;
    double norm = normal.computeNorm();
    if ( norm < 1.0e-8 ) {
        OOFEM_ERROR("Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.", 
          this->giveDofManager(2)->giveGlobalNumber(), this->giveDofManager(1)->giveGlobalNumber() )
    } else {
        this->normal = normal*(1.0/norm);
    }
    return 1;
}


void
Node2NodeContact :: computeGap(FloatArray &answer, TimeStep *tStep)
{
    // Computes the gap vector as the projection of the slave node penetration on a local coord system at the master node 
    // Normal penetration will be the third component
  
    // Compute slave penetration
    FloatArray xs, xm, uS, uM, ae;
    xs = *this->giveDofManager(1)->giveCoordinates();
    xm = *this->giveDofManager(2)->giveCoordinates();
    this->computeVectorOf( {D_u, D_v, D_w}, VM_Total, tStep, ae, true); // element solution vector
    uS = { ae.at(1), ae.at(2), ae.at(3) };
    uM = { ae.at(4), ae.at(5), ae.at(6) };
    xs.add(uS);
    xm.add(uM);
    FloatArray dx = xs-xm;
    
    // Compute gap as the projection on the local coordinate system at the master node
    FloatMatrix orthoBase;
    orthoBase.beLocalCoordSys( this->giveNormal() ); // create an othogonal base from the normal -> [t1^T; t2^T; n^T]
    answer.beProductOf(orthoBase, dx);
    //
}


void
Node2NodeContact :: computeCmatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *TimeStep)
{
    // TODO change name to vector/array
    // The normal is not updated for node2node which is for small deformations only
    // C = {n -n}
    FloatArray normal = this->giveNormal();
    answer = {  normal.at(1),  normal.at(2),  normal.at(3),
               -normal.at(1), -normal.at(2), -normal.at(3) };
    
}

void
Node2NodeContact :: computeTarraysAt(GaussPoint *gp, FloatArray &T1, FloatArray &T2, TimeStep *TimeStep)
{
    // Create T-arrays consisting of the tangent vectors e_T1 and e_T2
    // T1 = [e_T1, -e_T1] and T2 = [e_T2, -e_T2]
    FloatMatrix orthoBase;
    orthoBase.beLocalCoordSys( this->giveNormal() ); // create an othogonal base from the normal
    
    T1 = { orthoBase.at(1,1),  orthoBase.at(2,1),  orthoBase.at(3,1), 
          -orthoBase.at(1,1), -orthoBase.at(2,1), -orthoBase.at(3,1) };
    
    T2 = { orthoBase.at(1,2),  orthoBase.at(2,2),  orthoBase.at(3,2), 
          -orthoBase.at(1,2), -orthoBase.at(2,2), -orthoBase.at(3,2) };          
}


void
Node2NodeContact :: computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep)
{
    // TODO should be replaced with a call to constitutive model
    
    
    if ( gap.at(3) < 0.0 ) {
        StructuralInterfaceMaterial *mat = static_cast < StructuralInterfaceMaterial *> (this->giveContactMaterial() );
        mat->giveEngTraction_3d(t, gp, gap, tStep);
              
//         // Normal traction 
//         t = { 0.0, 0.0, this->epsN * gap.at(3) };
//         
//         // Add friction...
//         if ( this->giveContactDefinition()->giveContactMaterialNum() ) {
//             t.at(1) = this->epsT * gap.at(1);
//             t.at(2) = this->epsT * gap.at(2);
//             
//             
//         }
        
    } else {
        t = {0.0, 0.0, 0.0};
    }
  
}





void
Node2NodeContact :: computeContactForces(FloatArray &answer, TimeStep *tStep)
{
    answer.clear();
    FloatArray gap, C;
      
    this->computeGap(gap, tStep); // local system
    if ( gap.at(3) < 0.0 ) {
        GaussPoint *gp = this->integrationRule->getIntegrationPoint(0);
        FloatArray t;
        this->computeContactTractionAt(gp, t, gap, tStep); // local system
        
//         this->computeCmatrixAt(gp, C, tStep);
//         
//         // compute load vector
//         // fc = C^T * traction * A, Area - should be optional par
//         answer = t.at(3) * this->area * C;
//         
//         // Add friction forces
//         if ( this->giveContactDefinition()->giveContactMaterialNum() ) {
//             FloatArray T1, T2;
//             this->computeTarraysAt(gp, T1, T2, tStep);
//             answer.add( t.at(1) * this->area, T1 );
//             answer.add( t.at(2) * this->area, T2 );    
//         }
//         
        
        // new implementation
        FloatMatrix N;
        this->computeNmatrixAt(gp, N);
        
        FloatMatrix globalSys;
        globalSys.beLocalCoordSys( this->giveNormal() ); // should probably be stored so it is constant
        t.rotatedWith(globalSys, 't'); // transform to global system
        answer.beTProductOf(N, t);
        answer.times(this->area);
        t.printYourself("trac");
    }
  
}

  
void
Node2NodeContact :: computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep)
{
    GaussPoint *gp = this->integrationRule->getIntegrationPoint(0);
    FloatArray gap;
      
    this->computeGap(gap, tStep);
    if( gap.at(3) < 0.0 ) {
//         FloatArray C;
//         this->computeCmatrixAt(gp, C, tStep);
//         answer.beDyadicProductOf(C,C);
//         // this is the interface stiffness and should be obtained from that model
//         answer.times( this->epsN * this->area );
//         
        // add friction
    //     if ( this->giveContactDefinition()->giveContactMaterialNum() ) {
    //         FloatMatrix KT;
    //         this->computeFrictionTangent(KT, type, tStep);
    //         answer.add(KT);
    //     }
        
        FloatMatrix N, D;
        this->computeNmatrixAt(gp, N);
        
        
        StructuralInterfaceMaterial *mat = static_cast < StructuralInterfaceMaterial *> (this->giveContactMaterial() );
        mat->give3dStiffnessMatrix_Eng(D, TangentStiffness, gp, tStep); 
        
      
        
        FloatMatrix globalSys, DN;
        globalSys.beLocalCoordSys( this->giveNormal() ); // hould probably be stored so it is constant
        D.rotatedWith(globalSys, 't');                   // transform to global system
        DN.beProductOf(D,N);
        answer.clear();
        answer.plusProductUnsym(N, DN, this->area);
        
        //answer.printYourself("after");
                
        
    } else {
        answer.resize(6,6);
        answer.zero();
    }
    //answer.printYourself("after");
    
    
    
    if ( this->giveContactDefinition()->giveContactMaterialNum() ) {
        // new implementation
    }
}
  
  
void
Node2NodeContact :: computeFrictionTangent(FloatMatrix &answer, CharType type, TimeStep *tStep)
{
   // TODO change chartype
    GaussPoint *gp = this->integrationRule->getIntegrationPoint(0);
    StructuralInterfaceMaterial *mat = static_cast < StructuralInterfaceMaterial *> (this->giveContactMaterial() );
    mat->give3dStiffnessMatrix_Eng(answer, TangentStiffness, gp, tStep);
    /*
    // Stick case
    FloatArray T1, T2;
    this->computeTarraysAt(gp, T1, T2, tStep);
    answer.beDyadicProductOf(T1,T1);
    FloatMatrix temp;
    temp.beDyadicProductOf(T2,T2);
    answer.add(temp);
    answer.times( this->area * this->epsT );
    
    // Slip case*/
    
    
}

void
Node2NodeContact :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // N = [I, -I]
    answer.beNMatrixOf( {1, -1}, 3 );
}


void
Node2NodeContact :: giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s)
{
    // should return location array for the master and the slave node
    // TODO this whole thing should be rewritten  
  
  
    answer.resize(6);
    answer.zero();
    IntArray dofIdArray = {D_u, D_v, D_w};
    
    // master node
    for ( int i = 1; i <= dofIdArray.giveSize(); i++ ) {
        if ( this->giveDofManager(2)->hasDofID( (DofIDItem)dofIdArray.at(i) ) ) { // add corresponding number
            Dof *dof= this->giveDofManager(2)->giveDofWithID( (DofIDItem)dofIdArray.at(i) );
            answer.at(3+i) = s.giveDofEquationNumber( dof );
        } 
    }

    // slave node
    for ( int i = 1; i <= dofIdArray.giveSize(); i++ ) {
        if ( this->giveDofManager(1)->hasDofID( (DofIDItem)dofIdArray.at(i) ) ) { // add corresponding number
            Dof *dof= this->giveDofManager(1)->giveDofWithID( (DofIDItem)dofIdArray.at(i) );
            answer.at(i) = s.giveDofEquationNumber( dof );
        } 
    }    
}    


void
Node2NodeContact :: setupIntegrationPoints()
{
    // Sets up the integration rule array which contains all the necessary integration points
    if ( this->integrationRule == NULL ) {
        //TODO sets a null pointer for the element in the iRule 
        this->integrationRule = new GaussIntegrationRule(1, NULL) ;
        this->integrationRule->SetUpPoint(_Unknown);
    }
    
  
}













// node 2 node Lagrange


Node2NodeContactL :: Node2NodeContactL(int num, Domain *d, ContactDefinition *cDef) : Node2NodeContact(num, d, cDef)
{   

    this->area = 1.0;
};   



void
Node2NodeContactL :: giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s)
{
    //TODO maybe replace with element::giveLocationArray later but this requires to have separate methods for 1d,2d, 3d as
    // the dofid's must exist in the nodes 
    Node2NodeContact :: giveLocationArray(answer, s);
    
    // Add one lagrange dof
    if ( this->giveDofManager(1)->hasDofID( (DofIDItem)this->giveDofIdArray().at(1) ) ) { 
        Dof *dof= this->giveDofManager(1)->giveDofWithID( (DofIDItem)this->giveDofIdArray().at(1) );
        answer.followedBy( s.giveDofEquationNumber(dof) );
    }
    
}    



void
Node2NodeContactL :: computeContactForces(FloatArray &answer, TimeStep *tStep)
{
  
    //Loop through all the master objects and let them do their thing
    FloatArray gap, C, Fc;
    this->computeGap(gap, tStep);
    answer.resize( gap.giveSize() * 2 + 1);
    answer.zero();
        
    if( gap.at(3) < 0.0 ) {
    
        GaussPoint *gp = this->integrationRule->getIntegrationPoint(0);
        FloatArray t;
        this->computeContactTractionAt(gp, t ,gap, tStep);
        
        // new implementation
        FloatArray temp;
        FloatMatrix N;
        this->computeNmatrixAt(gp, N);
        
        FloatMatrix globalSys;
        globalSys.beLocalCoordSys( this->giveNormal() ); // should probably be stored so it is constant
        t.rotatedWith(globalSys, 't'); // transform to global system
        temp.beTProductOf(N, t);
        temp.times(this->area);

        answer.addSubVector(temp,1);
        answer.at( temp.giveSize() + 1 ) = gap.at(3);
    }
  
}
    


void
Node2NodeContactL :: computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep)
{
    answer.resize(7,7);
    answer.zero();
    
    FloatArray gap;
    this->computeGap(gap, tStep);
    
    if( gap.at(3) < 0.0 ) {
      
      GaussPoint *gp = this->integrationRule->getIntegrationPoint(0);
    
      FloatMatrix N, D;
      this->computeNmatrixAt(gp, N);
      
      StructuralInterfaceMaterial *mat = static_cast < StructuralInterfaceMaterial *> (this->giveContactMaterial() );
      mat->give3dStiffnessMatrix_Eng(D, TangentStiffness, gp, tStep); 
      D.at(3,3) = 0.0;
        
      FloatMatrix globalSys, DN, temp;
      globalSys.beLocalCoordSys( this->giveNormal() ); // should probably be stored so it is constant
      D.rotatedWith(globalSys, 't');                   // transform to global system
      DN.beProductOf(D,N);
      temp.clear();
      temp.plusProductUnsym(N, DN, this->area);
 
      FloatArray C;
      this->computeCmatrixAt(gp, C, tStep);
      int sz = C.giveSize();
      answer.resize(sz+1,sz+1);
      answer.zero();
      answer.assemble(temp, {1, 2, 3, 4, 5, 6}, {1, 2, 3, 4, 5, 6} );
      answer.addSubVectorCol(C, 1, sz + 1);
      answer.addSubVectorRow(C, sz + 1, 1);
        
      //answer.printYourself("after");
      
    }

    //TODO need to add a small number for the solver
    for ( int i = 1; i <= 7; i++ ) {
        answer.at(i,i) += 1.0e-6;
    }

    
    
    
}
  

  
void
Node2NodeContactL :: computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep)
{
    // should be replaced with a call to constitutive model
    // gap should be in a local system
    if ( gap.at(3) < 0.0 ) {
        StructuralInterfaceMaterial *mat = static_cast < StructuralInterfaceMaterial *> (this->giveContactMaterial() );
        mat->giveEngTraction_3d(t, gp, gap, tStep);
              
        Dof *dof = this->giveDofManager(1)->giveDofWithID( this->giveDofIdArray().at(1) );
        double lambda = dof->giveUnknown(VM_Total, tStep);
        t.at(3) = lambda;

    } else {
        t = {0.0, 0.0, 0.0};
    }
  
} 
  
    
void
Node2NodeContactL :: giveDofManagersToAppendTo(IntArray &answer)
{
    answer = {this->giveDofManager(1)->giveNumber()};
}
    
    
    
}












