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

#ifndef celnode2node_h
#define celnode2node_h

#include "Contact/structuralcontactelement.h"

///@name Input fields for _IFT_ContactElement
//@{
#define _IFT_Node2NodeContactP_Name "node2nodecontactp"

#define _IFT_Node2NodeContactL_Name "node2nodecontactl"
#define _IFT_Node2NodeContactL2_Name "node2nodecontactl2"
//@}

namespace oofem {
class Domain;
class ContactDefinition;
class ContactPair;

class OOFEM_EXPORT Node2NodeContact : public StructuralContactElement
{
protected:
    ContactDefinition *cDef;
    
private:
    
    // should be set by input:
    double area; // The area associated with the node (default = 1)- in order to represent some physical dimension.  
    
    FloatArray normal;
    
public:
    
    /// Constructor.
    Node2NodeContact(int num, Domain *d, ContactDefinition *cDef, ContactPair *cPair);
    /// Destructor.
    virtual ~Node2NodeContact(){};
    virtual void setupIntegrationPoints();
    virtual const char *giveInputRecordName() const { return _IFT_Node2NodeContactP_Name; }
    
};



class OOFEM_EXPORT Node2NodeContactL2 : public StructuralContactElementLagrange
{
private:
    int lagrangeId; // dof Id associated with the Lagrange multiplier

public:

    /// Constructor.
    Node2NodeContactL2(int num, Domain *d, ContactDefinition *cDef, ContactPair *cPair);
    /// Destructor.
    virtual ~Node2NodeContactL2(){};
    virtual const char *giveInputRecordName() const { return _IFT_Node2NodeContactL2_Name; }
    
    virtual void setupIntegrationPoints();
    
};






} // end namespace oofem
#endif // celnode2node_h
