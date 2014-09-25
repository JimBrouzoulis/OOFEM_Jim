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

#include "contact/contactelement.h"
#include "dofmanager.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "valuemodetype.h"
#include "dofiditem.h"
#include "timestep.h"
#include "dof.h"
#include "masterdof.h"
#include "unknownnumberingscheme.h"
#include "domain.h"
#include "gaussintegrationrule.h"
#include "error.h"
#include "element.h"

namespace oofem {

ContactElement :: ContactElement(int num, Domain *d, ContactDefinition *cDef) : Element(num, d)
{ 
    this->dofIdArray.clear();
    this->integrationRule = NULL;
    this->petrubedEquation = 0;
    this->cDef = cDef;
};  
  


void
ContactElement :: computeVectorOf(const IntArray &dofIDMask, ValueModeType u, TimeStep *tStep, FloatArray &answer, bool padding)
{
    Element :: computeVectorOf( dofIDMask, u, tStep, answer, padding);
    
    // Add numerical perturbation to an equation if active
    const double eps = 1.0e-8;
    if ( this->petrubedEquation ) {
        
        answer.at(this->petrubedEquation) += eps;
      
    }
  
  
  
}



void
ContactElement :: computeNumContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep)
{
    FloatArray f0, fp, df;
    this->computeContactForces(f0, tStep); 
    const double eps = 1.0e-8;
    
    int sz = f0.giveSize();
    answer.resize(sz,sz);
    for ( int  i = 1; i <= sz; i++ ) {
        this->addNumericalPerturbationToEq(i);
        this->computeContactForces(fp, tStep);   
        df.beDifferenceOf(fp, f0);
        answer.addSubVectorCol(df, 1, i);
    }
    this->resetNumericalPerturbation();
    answer.times(1.0 / eps);
    
    
    
}



}












