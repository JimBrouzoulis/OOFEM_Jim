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

#ifndef contactdefinitionnode2node_h
#define contactdefinitionnode2node_h

#include "contact/contactdefinition.h"
#include "intarray.h"

///@name Input fields for _IFT_ContactDefinitionNode2Node
//@{
#define _IFT_ContactDefinitionNode2Edge_Name "cdef_node2edgep"
#define _IFT_ContactDefinitionNode2Edge_MasterElements "masterelements"
#define _IFT_ContactDefinitionNode2Edge_MasterEdges "masteredges"
#define _IFT_ContactDefinitionNode2Edge_SlaveNodes "slavenodes"
#define _IFT_ContactDefinitionNode2Edge_Lagrange "lagrange"
//@}

namespace oofem {
class Domain;
class ContactManager;
class ContactObject;
class ContactElement;


/**
 * This class manages a node 2 edge contact pair
 * @author Jim Brouzoulis
 */
class OOFEM_EXPORT ContactDefinitionNode2Edge : public ContactDefinition
{
protected:
  // Store them in a more flexible way - a ContactPair struct?
  IntArray masterEdges;
  IntArray masterElements;
  IntArray slavenodes;
  
public:

    /// Constructor.
    ContactDefinitionNode2Edge(ContactManager *cMan);
    /// Destructor.
    virtual ~ContactDefinitionNode2Edge(){};

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int instanciateYourself(DataReader *dr);
    
    virtual const char *giveClassName() const { return "ContactDefinitionNode2Edge"; };
    virtual const char *giveInputRecordName() const { return _IFT_ContactDefinitionNode2Edge_Name; };
 
};



} // end namespace oofem
#endif // contactdefinitionnode2node_h
