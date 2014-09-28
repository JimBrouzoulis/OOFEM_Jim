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

#ifndef contactpair_h
#define contactpair_h

#include "datareader.h"
#include "floatarray.h"

namespace oofem {
class Element;
class DofManager; 
class Node;
class TimeStep;

/**
 * Keeps track of the pairing node2node etc and should compute CPP, transformations etc
 *
 * @author Jim Brouzoulis
 */
class OOFEM_EXPORT ContactPair
{
public:

    ContactPair(){};
    virtual ~ContactPair(){};

    virtual int instanciateYourself(DataReader *dr){ return IRRT_OK; };
    virtual const char *giveClassName() const { return "ContactPair"; }
    
    
};


class OOFEM_EXPORT ContactPairNode2Edge : public ContactPair
{   
public:

    ContactPairNode2Edge(Element *el, int edge, Node *slave);
    virtual ~ContactPairNode2Edge() {};

    virtual int instanciateYourself(DataReader *dr);
    virtual const char *giveClassName() const { return "ContactPairNode2Edge"; }
    
    
    //element methods
    virtual void computeNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer);
    virtual void computeCovarBaseVectorAt(const FloatArray &lCoords, FloatArray &g, TimeStep *tStep);
    virtual void performCPP(TimeStep *tStep);
    virtual void computeGap(FloatArray &answer, FloatArray &lCoords, TimeStep *tStep);
    // gp method
    double giveCPPcoord() { return this->xibar; };
    Element *giveMasterElement() { return this->masterElement; };
private:
    Element *masterElement;   // the element to which the edge belongs TODO what to do with shared edges? 
    int masterElementEdgeNum; // local edge number on the master element 
    Node *slaveNode;
    
    std::vector< Node* > masterNodes;
    double xibar; // local coordinate from CPP
    
};
    





} // end namespace oofem
#endif // contactpair
