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

// node2edge
#include "feinterpol2d.h"

namespace oofem {

ContactPairNode2Edge :: ContactPairNode2Edge(Element *el, int edge, Node *slave) : ContactPair()
{
    this->masterElement = el;
    this->masterElementEdgeNum = edge;
    this->slaveNode = slave;
    this->xibar = -2.0; // outside element
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
ContactPairNode2Edge :: computeNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer)
{
    //NOTE General method to be reused
    // N = [I, -N_edge]
    FEInterpolation2d *interp = static_cast< FEInterpolation2d* > (this->masterElement->giveInterpolation() );
    FloatArray Nedge, N = {1};
    FloatMatrix NedgeMat;
    interp->edgeEvalN(Nedge, this->masterElementEdgeNum, lCoords, FEIElementGeometryWrapper(this->masterElement) );
    Nedge.negated();
    N.append(Nedge);
    answer.beNMatrixOf(N, 3);

}    


void
ContactPairNode2Edge :: computeCovarBaseVectorAt(const FloatArray &lCoords, FloatArray &g, TimeStep *tStep)
{
     //NOTE General method for 2D
     // Computes the updated covariant base vector (tangent vector) for a 2D edge
     FloatArray dNdxi;
     FEInterpolation2d *interp = static_cast< FEInterpolation2d* > (this->masterElement->giveInterpolation() );
     interp->edgeEvaldNdxi( dNdxi, this->masterElementEdgeNum, lCoords, FEIElementGeometryWrapper(this->masterElement) );
     g.resize(3);
     g.zero();
     FloatArray x;
     for ( int i = 1; i <= dNdxi.giveSize(); i++ ) {
         x = *this->masterNodes[i-1]->giveCoordinates();
         g.add( dNdxi.at(i), x);
     }
}

 
void
ContactPairNode2Edge :: performCPP(TimeStep *tStep)
{
    // should work on a local slave coord-> compute global slave coord 
    // (for when we only have slave gp's instead of actual slave nodes)
    // Compute updated coordinates for all the involved points
    FloatArray xs, xm1, xm2;
    this->slaveNode->giveUpdatedCoordinates(xs, tStep);
    this->masterNodes[0]->giveUpdatedCoordinates(xm1, tStep);
    this->masterNodes[1]->giveUpdatedCoordinates(xm2, tStep);
    
    // could x be created in a better way?
    FloatArray x, xibar;
    x = xs;
    x.append(xm1);
    x.append(xm2);
    this->computeCPP( xibar, x );
    this->xibar = xibar.at(1); // why not store the array?
  
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
    answer = { -1.0 + xibar * 2.0 }; // transform from [0,1] to [-1,1]
}

void
ContactPairNode2Edge :: computeLinearizationOfCPP(const FloatArray &lCoords, FloatMatrix &answer, TimeStep *tStep)
{
    // Numerically compute the linearization of a CPP wrt contact element edgeNodes
    //TODO Add analytical version for straight segment
  
    // Compute updated coordinates for all the involved points
    FloatArray xs, xm1, xm2;
    this->slaveNode->giveUpdatedCoordinates(xs, tStep);
    this->masterNodes[0]->giveUpdatedCoordinates(xm1, tStep);
    this->masterNodes[1]->giveUpdatedCoordinates(xm2, tStep);
    
    // Numerical derivative
    FloatArray x0, x, xi0, xi, dxi;
    x0 = xs; x0.append(xm1); x0.append(xm2);
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

void
ContactPairNode2Edge :: computeGap(FloatArray &answer, FloatArray &lCoords, TimeStep *tStep)
{
    // Computes the gap vector from the CCP (closest point projection) of 
    // the slave node onto the master edge.
    // gap = xs - xm(xibar) = xs - N_i(xibar) * xm_i
  
    this->slaveNode->giveUpdatedCoordinates(answer, tStep);
    
    FEInterpolation2d *interp = static_cast< FEInterpolation2d* > (this->masterElement->giveInterpolation()); 
    FloatArray N, xi;    
    interp->edgeEvalN(N, this->masterElementEdgeNum, lCoords, FEIElementGeometryWrapper(this->masterElement) );
    for ( int i = 1; i <= N.giveSize(); i++ ) {
        this->masterNodes[i-1]->giveUpdatedCoordinates(xi, tStep);
        answer.add( -N.at(i), xi );
    }
    
//     // Rotate gap to a local system
//     FloatMatrix orthoBase;
//     computeCurrentTransformationMatrixAt( lCoords, orthoBase, tStep );
//     answer.rotatedWith(orthoBase, 't');

}




}
