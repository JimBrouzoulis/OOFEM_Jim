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

#include "Contact/contactpair.h"
#include "element.h"
#include "dofmanager.h"
#include "timestep.h"
#include "gausspoint.h"

// node2edge
#include "feinterpol2d.h"
#include "CrossSections/structuralcrosssection.h"

namespace oofem {

  
void
ContactPair :: computeCurrentNormalAt(const FloatArray &lCoords, FloatArray &normal, TimeStep *tStep)  
{
    FloatArray g1, g2;
    this->computeCovarTangentVectorsAt( lCoords, g1, g2, tStep);
    normal.beVectorProductOf(g1, g2);
    normal.normalize(); 
}
  
  
double 
ContactPair::computeCurrentAreaAround(GaussPoint* gp, TimeStep* tStep)
{
    double weight  = gp->giveWeight();
    FloatArray g1, g2, g3;
    this->computeCovarTangentVectorsAt( *gp->giveNaturalCoordinates(), g1, g2, tStep);
    g3.beVectorProductOf(g1, g2);
    return g3.computeNorm() * weight;
}

  
void 
ContactPair :: computeCurrentTransformationMatrixAt(const FloatArray &lCoords, FloatMatrix &answer, TimeStep *tStep) 
{
    // Transformation matrix to the local coordinate system
    FloatArray normal;
    this->computeCurrentNormalAt(lCoords, normal, tStep);
    answer.beLocalCoordSys( normal );
        
}

void
ContactPair :: computeNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer)
{
    // N = [N_slave, -N_master]
   
    FloatArray Nslave, Nmaster;
    giveSlaveNarray(lCoords, Nslave);
    giveMasterNarray(lCoords, Nmaster);
    Nmaster.negated();
    Nslave.append(Nmaster);

    answer.beNMatrixOf(Nslave, 3);

}
  
  

void
ContactPair :: giveCurrentCoordsArray(FloatArray &answer, TimeStep *tStep)
{
    // computes the array with updated node coords for slave and master nodes
    int numDofs = ( this->giveNumberOfMasterNodes() + this->giveNumberOfSlaveNodes() ) * 3;
    answer.resize(numDofs);
    answer.zero();

    FloatArray x;
    int pos = 0;
    for ( int i = 1; i <= this->giveNumberOfSlaveNodes(); i++ ) {
      
        this->giveSlaveNode(i)->giveUpdatedCoordinates(x, tStep);
        for ( int j = 1; j <= x.giveSize(); j++ ) {
            answer.at(pos + j) = x.at(j);
        }
        pos += 3;
        
    }
    
    for ( int i = 1; i <= this->giveNumberOfMasterNodes(); i++ ) {
      
        this->giveMasterNode(i)->giveUpdatedCoordinates(x, tStep);
        for ( int j = 1; j <= x.giveSize(); j++ ) {
            answer.at(pos + j) = x.at(j);
        }
        pos += 3;
        
    }
}
  
  
  
  
void
ContactPair :: performCPP(GaussPoint *gp, TimeStep *tStep)
{
    // should work on a local slave coord-> compute global slave coord 
    // (for when we only have slave gp's instead of actual slave nodes)
    // Compute updated coordinates for all the involved points

    FloatArray x, xibar;
    this->giveCurrentCoordsArray(x, tStep);
   
    this->computeCPP( xibar, x );   
    gp->setNaturalCoordinates( xibar );
}  
  
  
void
ContactPair :: computeGap(FloatArray &answer, FloatArray &lCoords, TimeStep *tStep)
{
    // Computes the gap vector from the CCP (closest point projection) of 
    // the slave node onto the master edge.
    // gap = xs(xibar) - xm(xibar) = [N_slave, -N_master] * xhat 
    FloatArray x;
    FloatMatrix N;
    this->giveCurrentCoordsArray(x, tStep);
    this->computeNmatrixAt(lCoords, N);
    answer.beProductOf(N, x);
}  
  
  
//--------------------------------  
  
  
  
ContactPairNode2Edge :: ContactPairNode2Edge(Element *el, int edge, Node *slave) : ContactPair()
{
    this->masterElement = el;
    this->masterElementEdgeNum = edge;
    this->slaveNodes.resize(1);
    this->slaveNodes[0] = slave;
    //this->xibar = -2.0; // outside element
}
    

int
ContactPairNode2Edge :: instanciateYourself(DataReader *dr)
{
  
    IntArray edgeNodes, globalNodeArray;
    FEInterpolation2d *interp = static_cast< FEInterpolation2d* > ( this->masterElement->giveInterpolation() );
    IntArray elNodes = this->masterElement->giveDofManArray(); 
    interp->computeEdgeMapping(edgeNodes, elNodes, this->masterElementEdgeNum );
    
    this->masterNodes.resize( edgeNodes.giveSize() );
    for ( int j = 1; j<= edgeNodes.giveSize(); j++ ) {        
        masterNodes[j-1] = this->masterElement->giveNode( edgeNodes.at(j) );
    }        
    
    
  return 1;
}    
    
    


void
ContactPairNode2Edge :: giveSlaveNarray(const FloatArray &lCoords, FloatArray &answer)
{
    answer = {1};
}


void
ContactPairNode2Edge :: giveMasterNarray(const FloatArray &lCoords, FloatArray &answer)
{
    FEInterpolation2d *interp = static_cast< FEInterpolation2d* > (this->masterElement->giveInterpolation() );
    interp->edgeEvalN(answer, this->masterElementEdgeNum, lCoords, FEIElementGeometryWrapper( this->masterElement ) );
}

 
void
ContactPairNode2Edge :: computeCovarTangentVectorsAt(const FloatArray &lCoords, FloatArray &g1, FloatArray &g2, TimeStep *tStep)
{
  
    
     // Computes the updated covariant tangent (base) vectors
     FloatArray dNdxi;
     FEInterpolation2d *interp = static_cast< FEInterpolation2d* > (this->masterElement->giveInterpolation() );
     interp->edgeEvaldNdxi( dNdxi, this->masterElementEdgeNum, lCoords, FEIElementGeometryWrapper(this->masterElement) );
     g1.resize(3); g1.zero();
     FloatArray x;
     for ( int i = 1; i <= dNdxi.giveSize(); i++ ) {
         x = *this->masterNodes[i-1]->giveCoordinates();
         g1.add( dNdxi.at(i), x);
     }
     StructuralCrossSection *cs = static_cast< StructuralCrossSection* > (  this->masterElement->giveCrossSection() );
     g2 = {0.0, 0.0, cs->give(CS_Thickness, NULL) }; // The second one point is in the out of plane dir. with t as length
}



void
ContactPairNode2Edge :: computeCPP(FloatArray &answer, const FloatArray &x)
{
    FloatArray xs  = { x.at(1), x.at(2), x.at(3) };
    FloatArray xm1 = { x.at(4), x.at(5), x.at(6) };
    FloatArray xm2 = { x.at(7), x.at(8), x.at(9) };
    
    // Perform CCP to find xibar - this is for a straight segment
    // TODO For now assume it to be a straight edge - generalize later
    // xibar = 1/l² * (xs-xm1).(xm2-xm1)  
    FloatArray a = xm2 - xm1;
    FloatArray dx = xs - xm1;
    double l2 = a.computeSquaredNorm();
    double xibar = dx.dotProduct(a) / l2;
    if ( xibar >= 0.0  &&  xibar <= 1.0) { 
        answer = { -1.0 + xibar * 2.0 }; // transform from [0,1] to [-1,1]
    } else {
        answer = {}; // outside
    }
      
}

void
ContactPairNode2Edge :: computeLinearizationOfCPP(const FloatArray &lCoords, FloatMatrix &answer, TimeStep *tStep)
{
    // Numerically compute the linearization of a CPP wrt contact element edgeNodes
    //TODO Add analytical version for straight segment
  
    // Numerical derivative
    FloatArray x0, x, xi0, xi, dxi;
    this->giveCurrentCoordsArray(x0, tStep);
    //x0 = xs; x0.append(xm1); x0.append(xm2);
    this->computeCPP( xi0, x0 );
  
    answer.resize( 1, x0.giveSize() );
    const double eps = 1.0e-8;
    for ( int i = 1; i <= x0.giveSize(); i++ ) {
        x = x0;
        x.at(i) += eps;
        this->computeCPP( xi, x );
        dxi.beDifferenceOf( xi, xi0 );
        answer.addSubVectorCol( dxi, 1, i );
    }
    answer.times( 1.0 / eps );
    
}


void
ContactPairNode2Edge :: computeBmatrixAt(const FloatArray &lCoords, const FloatArray &traction, FloatMatrix &answer, TimeStep *tStep)
{

    FloatArray dNdxi;
    FEInterpolation2d *interp = static_cast< FEInterpolation2d* > (this->masterElement->giveInterpolation() );
    interp->edgeEvaldNdxi( dNdxi, this->masterElementEdgeNum, lCoords, FEIElementGeometryWrapper(this->masterElement) );
     
    FloatMatrix dxidx;
    this->computeLinearizationOfCPP(lCoords, dxidx, tStep);
    FloatArray v;
    FloatMatrix temp;
    answer.resize(dxidx.giveNumberOfColumns(), dxidx.giveNumberOfColumns());
    answer.zero();
    int nno = 2;
    for ( int i = 1; i <= nno; i++ ) {
        v = { -dNdxi.at(i) * dxidx.at(1,i), -dNdxi.at(i) * dxidx.at(1,i+1), -dNdxi.at(i) * dxidx.at(1,i+2) };
        temp.beDyadicProductOf( traction, v);
        IntArray pos = {(i-1)*3 + 1, (i-1)*3 + 2 , (i-1)*3 + 3};
        pos.add(3);
        answer.assemble(temp, pos, pos );
     }
     
     
}    





//-------------------------------------------

ContactPairNode2Node :: ContactPairNode2Node(Node *master, Node *slave)
{
    /*
    this->masterNode = master;
    this->slaveNode  = slave;*/
    this->masterNodes.resize(1);
    this->masterNodes[0] = master;
    this->slaveNodes.resize(1);
    this->slaveNodes[0] = slave;
}


void
ContactPairNode2Node :: computeCovarTangentVectorsAt(const FloatArray &lCoords, FloatArray &g1, FloatArray &g2, TimeStep *tStep)
{
    FloatMatrix orthoSys;
    FloatArray normal;
    this->computeCurrentNormalAt(lCoords, normal, tStep);
    orthoSys.beLocalCoordSys(normal);
    g1.beColumnOf(orthoSys, 1);
    g2.beColumnOf(orthoSys, 2);
}



void
ContactPairNode2Node :: computeCurrentNormalAt(const FloatArray &lCoords, FloatArray &normal, TimeStep *tStep)
{
 
    // compute normal as direction vector from master node to slave node
    FloatArray xs, xm;
    xs = *this->giveSlaveNode(1)->giveCoordinates();
    xm = *this->giveMasterNode(1)->giveCoordinates();
    normal = xs-xm;
    double norm = normal.computeNorm();
    if ( norm < 1.0e-8 ) {
        OOFEM_ERROR("Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.", 
          this->giveMasterNode(1)->giveGlobalNumber(), this->giveSlaveNode(1)->giveGlobalNumber() )
    } else {
        normal.times( 1.0 / norm );
    }    
    
}


int
ContactPairNode2Node :: instanciateYourself(DataReader *dr)
{
  
//     IntArray edgeNodes, globalNodeArray;
//     FEInterpolation2d *interp = static_cast< FEInterpolation2d* > ( this->masterElement->giveInterpolation() );
//     IntArray elNodes = this->masterElement->giveDofManArray(); 
//     interp->computeEdgeMapping(edgeNodes, elNodes, this->masterElementEdgeNum );
//     
//     this->masterNodes.resize( edgeNodes.giveSize() );
//     for ( int j = 1; j<= edgeNodes.giveSize(); j++ ) {        
//         masterNodes[j-1] = this->masterElement->giveNode( edgeNodes.at(j) );
//     }        
//     
    
  return 1;
}    


}
