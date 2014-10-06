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

#ifndef celnode2node_h
#define celnode2node_h

#include "Contact/structuralcontactelement.h"

///@name Input fields for _IFT_Node2EdgeContact
//@{
#define _IFT_Node2EdgeContactP_Name "node2edgecontactp"

//@}

namespace oofem {
class Domain;
class ContactDefinition;
class ContactPair;

class OOFEM_EXPORT Node2EdgeContact : public StructuralContactElement
{
public:
    Node2EdgeContact(int num, Domain *d, ContactDefinition *cDef, ContactPair *cPair);
    virtual ~Node2EdgeContact(){};
    
    virtual void setupIntegrationPoints();

    virtual const char *giveInputRecordName() const { return _IFT_Node2EdgeContactP_Name; }
    
    //additional tangent from changes in contact position ->change name
    //virtual void computeBmatrixAt(const FloatArray &lCoords, const FloatArray &traction, FloatMatrix &answer, TimeStep *tStep);
    
};



class OOFEM_EXPORT Node2EdgeContactL : public StructuralContactElementLagrange
{
public:
    Node2EdgeContactL(int num, Domain *d, ContactDefinition *cDef, ContactPair *cPair);
    virtual ~Node2EdgeContactL(){};
    
    virtual void setupIntegrationPoints();

    virtual const char *giveInputRecordName() const { return _IFT_Node2EdgeContactP_Name; }
    
    //additional tangent from changes in contact position ->change name
    //virtual void computeBmatrixAt(const FloatArray &lCoords, const FloatArray &traction, FloatMatrix &answer, TimeStep *tStep);
    
};






} // end namespace oofem
#endif // celnode2node_h
