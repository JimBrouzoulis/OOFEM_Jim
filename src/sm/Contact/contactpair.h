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
#include "floatmatrix.h"

namespace oofem {
class Element;
class DofManager; 
class Node;
class TimeStep;
class GaussPoint;
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
    
    virtual Node *giveMasterNode(const int num) { return this->masterNodes[num-1]; };
    virtual Node *giveSlaveNode(const int num) { return this->slaveNodes[num-1]; };
    int giveNumberOfMasterNodes() { return this->masterNodes.size(); };
    int giveNumberOfSlaveNodes() { return this->slaveNodes.size(); };
    
    FloatArray &giveCPPcoords(const int num) { return this->CPPpoints[num-1]; };
    virtual void computeCovarTangentVectorsAt(const FloatArray &lCoords, FloatArray &g1, FloatArray &g2, TimeStep *tStep){};
    virtual void computeCurrentNormalAt(const FloatArray &lCoords, FloatArray &normal, TimeStep *tStep);
    virtual void computeCurrentTransformationMatrixAt(const FloatArray &lCoords, FloatMatrix &answer, TimeStep *tStep);
    virtual double computeCurrentAreaAround(GaussPoint *gp, TimeStep *tStep);
    
    virtual void performCPP(GaussPoint *gp, TimeStep *tStep);
    virtual void computeCPP(FloatArray &answer, const FloatArray &x){};
    
    virtual void computeGap(FloatArray &answer, FloatArray &lCoords, TimeStep *tStep);
    
    
    virtual void computeNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer);
    virtual void giveSlaveNarray(const FloatArray &lCoords, FloatArray &answer) { answer.clear(); };
    virtual void giveMasterNarray(const FloatArray &lCoords, FloatArray &answer) { answer.clear(); };;
    virtual void giveCurrentCoordsArray(FloatArray &answer, TimeStep *tStep);

protected:
    std :: vector< Node* > masterNodes;
    std :: vector< Node* > slaveNodes;
    
    std :: vector< FloatArray > CPPpoints; // list of all the CPP the cPair have
};





class OOFEM_EXPORT ContactPairNode2Edge : public ContactPair
{   
public:

    ContactPairNode2Edge(Element *el, int edge, Node *slave);
    virtual ~ContactPairNode2Edge() {};

    virtual int instanciateYourself(DataReader *dr);
    virtual const char *giveClassName() const { return "ContactPairNode2Edge"; }
    
    
    //element methods
    virtual void giveSlaveNarray(const FloatArray &lCoords, FloatArray &answer);
    virtual void giveMasterNarray(const FloatArray &lCoords, FloatArray &answer);
    
    
    virtual void computeBmatrixAt(const FloatArray &lCoords, const FloatArray &traction, FloatMatrix &answer, TimeStep *tStep);
    virtual void computeCovarTangentVectorsAt(const FloatArray &lCoords, FloatArray &g1, FloatArray &g2, TimeStep *tStep);
    
    virtual void computeCPP(FloatArray &answer, const FloatArray &x);
    virtual void computeLinearizationOfCPP(const FloatArray &lCoords, FloatMatrix &answer, TimeStep *tStep);
    

    Element *giveMasterElement() { return this->masterElement; };
private:
    Element *masterElement;   // the element to which the edge belongs TODO what to do with shared edges? 
    int masterElementEdgeNum; // local edge number on the master element 

    
    double xibar; // local coordinate from CPP - remove
    
};
    


class OOFEM_EXPORT ContactPairNode2Node : public ContactPair
{   
public:

    ContactPairNode2Node(Node *master, Node *slave);
    virtual ~ContactPairNode2Node() {};

    virtual int instanciateYourself(DataReader *dr);
    virtual const char *giveClassName() const { return "ContactPairNode2Node"; }
    
    
    //element methods
    virtual void giveSlaveNarray(const FloatArray &lCoords, FloatArray &answer) { answer = {1.0}; };
    virtual void giveMasterNarray(const FloatArray &lCoords, FloatArray &answer) { answer = {1.0}; };
    
    
    //virtual void computeBmatrixAt(const FloatArray &lCoords, const FloatArray &traction, FloatMatrix &answer, TimeStep *tStep);
    virtual void computeCovarTangentVectorsAt(const FloatArray &lCoords, FloatArray &g1, FloatArray &g2, TimeStep *tStep);
    virtual void computeCurrentNormalAt(const FloatArray &lCoords, FloatArray &normal, TimeStep *tStep);
    
    virtual void computeCPP(FloatArray &answer, const FloatArray &x) { answer = {0.0, 0.0, 0.0}; };
    //virtual void computeLinearizationOfCPP(const FloatArray &lCoords, FloatMatrix &answer, TimeStep *tStep);
  
    virtual double computeCurrentAreaAround(GaussPoint *gp, TimeStep *tStep) {return 1.0; }; //TODO
};



} // end namespace oofem
#endif // contactpair
