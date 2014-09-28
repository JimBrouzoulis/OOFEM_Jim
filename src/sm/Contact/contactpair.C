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

ContactPairNode2Edge :: ContactPairNode2Edge(Element *el, int edge, DofManager *slave) : ContactPair()
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
        masterNodes[j-1] = this->masterElement->giveDofManager(edgeNodes.at(j));
    }        

    
    // TODO TEMP - compute CPP
    this->xibar = -1.0;
    
    
    
  return 1;
}    
    
    
void
ContactPairNode2Edge :: computeNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer)
{
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
  
      // Updated covar base vector
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

      
    // Compute updaeted coordinates x = X + u
    FloatArray xs, xm1, xm2, us, um1, um2, ae;
    xs  = *this->slaveNode->giveCoordinates();
    xm1 = *this->masterNodes[0]->giveCoordinates();
    xm2 = *this->masterNodes[1]->giveCoordinates();
    
    this->masterElement->computeVectorOf( {D_u, D_v, D_w}, VM_Total, tStep, ae, true); // element solution vector
    us = { ae.at(1), ae.at(2), ae.at(3) };
    um1 = { ae.at(4), ae.at(5), ae.at(6) };
    um2 = { ae.at(7), ae.at(8), ae.at(9) };
    xs.add(us);
    xm1.add(um1);
    xm2.add(um2);
    
    // Perform CCP to find xibar
    // TODO For now assume it to be a straight edge - generalize later
    // xibar = 1/l² * (xs-xm1).(xm2-xm1)  
    FloatArray a = xm2 - xm1;
    FloatArray dx = xs - xm1;
    double l2 = a.computeSquaredNorm();
    double xibar = dx.dotProduct(a) / l2;
    this->xibar = -1.0 + xibar * 2.0; // transform from [0,1] to [-1,1]
}




}
