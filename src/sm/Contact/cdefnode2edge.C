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

#include "Contact/cdefnode2edge.h"
#include "domain.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "feinterpol2d.h"

#include "Contact/celnode2edge.h"
#include "Contact/contactpair.h"

namespace oofem {
REGISTER_ContactDefinition(ContactDefinitionNode2Edge)


ContactDefinitionNode2Edge :: ContactDefinitionNode2Edge(ContactManager *cMan) : ContactDefinition(cMan){}
    

IRResultType
ContactDefinitionNode2Edge :: initializeFrom(InputRecord *ir)
{
 
    IRResultType result; // Required by IR_GIVE_FIELD macro
    ContactDefinition :: initializeFrom(ir);
    
    IntArray masterElements;
    IntArray masterEdges;
    IntArray slaveNodes;
    IR_GIVE_FIELD(ir, masterElements, _IFT_ContactDefinitionNode2Edge_MasterElements);
    IR_GIVE_FIELD(ir, masterEdges, _IFT_ContactDefinitionNode2Edge_MasterEdges);
    IR_GIVE_FIELD(ir, slaveNodes, _IFT_ContactDefinitionNode2Edge_SlaveNodes);
    
    Domain *domain = this->giveContactManager()->giveDomain();
    for( int i = 1; i<= slaveNodes.giveSize(); i++ ) {
      
        // create a node2edge contact pair
        Element *el = domain->giveElement( masterElements.at(i) );
        ContactPair *cPair = new ContactPairNode2Edge(el, masterEdges.at(i), domain->giveDofManager( slaveNodes.at(i) ) );
        cPair->instanciateYourself(NULL); //TODO fix!
        
        IntArray edgeNodes, globalNodeArray;
        FEInterpolation2d *interp = static_cast< FEInterpolation2d* > ( el->giveInterpolation() );
        IntArray elNodes = el->giveDofManArray(); 
        interp->computeEdgeMapping(edgeNodes, elNodes, masterEdges.at(i) );
        
        
        globalNodeArray.resize(edgeNodes.giveSize() + 1);
        globalNodeArray.at(1) = slaveNodes.at(i);
        for ( int j = 1; j<= edgeNodes.giveSize(); j++ ) {        
            globalNodeArray.at(j+1) = el->giveDofManagerNumber(edgeNodes.at(j));
        }        
        
        
        ContactElement *cEl = new Node2EdgeContact(i, domain, this, cPair);
        
        
        // initialize contact element from dynamic input record - this is to be able to use regular element functions like giveLocArray
        DynamicInputRecord *dir;
        dir = CreateElementIR(i, _IFT_ContactDefinitionNode2Edge_Name, globalNodeArray );
        cEl->initializeFrom(dir);
        
        this->addContactElement(cEl);
    }
    
    return IRRT_OK;
}

int
ContactDefinitionNode2Edge :: instanciateYourself(DataReader *dr)
{
    ContactDefinition :: instanciateYourself(dr);
    return 1;
}


}
