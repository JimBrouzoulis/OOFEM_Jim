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

#ifndef contactelement_h
#define contactelement_h

#include "chartype.h"
#include "valuemodetype.h"
#include "element.h"
#include "contact/contactdefinition.h"
#include <vector>


///@name Input fields for _IFT_ContactElement
//@{
//#define _IFT_ContactManager_Name "contactmanager"

//@}

namespace oofem {
class Domain;
class ContactManager;
class SparseMtrx;
class TimeStep;
class DofManager;
class GaussPoint;
class UnknownNumberingScheme;
class FloatMatrix;
class IntegrationRule;

// contact elements

class OOFEM_EXPORT ContactElement : public Element
{
private:
    // Pointer to its contact definition
    ContactDefinition *cDef;

    
    //std :: vector< ContactElement *> slaveObjectList; // remove?
    
    IntArray dofIdArray; // stores dofid's the element adds
    

    
protected:
    bool inContact;
    int petrubedEquation;
    
public:
    IntegrationRule *integrationRule;
    /// Constructor.
    ContactElement(int num, Domain *d, ContactDefinition *cDef);
    /// Destructor.
    virtual ~ContactElement(){};

    virtual void setupIntegrationPoints(){};
    
    ContactDefinition *giveContactDefinition() { return this->cDef; };
    
    Material *giveContactMaterial() { return this->cDef->giveContactMaterial(); }; 
    virtual int instanciateYourself(DataReader *dr){ return 1; };
    //virtual const char *giveClassName() const { return "ContactDefinition"; }
    bool isInContact() { return inContact; };
    virtual void giveDofManagersToAppendTo(IntArray &answer) { answer = {}; }; 
    
    
    virtual void computeContactForces(FloatArray &answer, TimeStep *tStep) = 0;   
    
    virtual void computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep) = 0;
    
    // numerical tangent
    virtual void computeNumContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep);
    void addNumericalPerturbationToEq(int eq) { this->petrubedEquation = eq; };
    void resetNumericalPerturbation() { this->petrubedEquation = 0; };
    virtual void computeVectorOf(const IntArray &dofIDMask, ValueModeType u, TimeStep *tStep, FloatArray &answer, bool padding = false);
    
    
    
    virtual void giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s) = 0;
    // set the dof id array if the contact element has its own dofs (Lagrange multipliers)
    virtual void setDofIdArray(IntArray &array){ this->dofIdArray = array; };
    virtual IntArray &giveDofIdArray(){ return this->dofIdArray; };
};


/*
class OOFEM_EXPORT Node2NodeContact : public ContactElement
{
protected:
    ContactDefinition *cDef;

private:
    DofManager *masterNode;
    DofManager *slaveNode;
    
    // should be set by input:
    double area; // The area associated with the node (default = 1)- in order to represent some physical dimension.  
    double epsN;  // penalty stiffness 
    
    FloatArray normal;
public:

    /// Constructor.
    Node2NodeContact(DofManager *master, DofManager *slave);
    /// Destructor.
    virtual ~Node2NodeContact(){};
    virtual int instanciateYourself(DataReader *dr);
    virtual void setupIntegrationPoints();
    
    virtual void computeGap(FloatArray &answer, TimeStep *tStep);
    virtual void computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep);
    virtual void computeCmatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *TimeStep);
    FloatArray &giveNormal() { return this->normal; };
    
    
    // Necessary methods - pure virtual in base class
    virtual void computeContactForces(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms);    
    
    virtual void computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep);
    
    virtual void giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s);
};




class OOFEM_EXPORT Node2NodeContactL : public Node2NodeContact
{
protected:
    ContactDefinition *cDef;

private:
    DofManager *masterNode;
    DofManager *slaveNode;
    int lagrangeId; // dof Id associated with the Lagrange multiplier
    
    // should be set by input:
    double area; // The area associated with the node (default = 1)- in order to represent some physical dimension.  
  
    
public:

    /// Constructor.
    Node2NodeContactL(DofManager *master, DofManager *slave);
    /// Destructor.
    virtual ~Node2NodeContactL(){};
    //virtual int instanciateYourself(DataReader *dr);
    virtual void giveDofManagersToAppendTo(IntArray &answer); 
    virtual void computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep);
    
    //virtual void computeGap(FloatArray &answer, TimeStep *tStep);
    //virtual void computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep);
    //virtual void computeCmatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *TimeStep);
    
    // Necessary methods - pure virtual in base class
    virtual void computeContactForces(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms);    
    
    virtual void computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep);
    
    virtual void giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s);
};*/



} // end namespace oofem
#endif // contactelement_h
