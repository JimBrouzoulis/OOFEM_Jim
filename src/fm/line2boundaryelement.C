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

#include "line2boundaryelement.h"
#include "node.h"
#include "dof.h"
#include "fei2dlinequad.h"
#include "gaussintegrationrule.h"
#include "crosssection.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Line2BoundaryElement);

FEI2dLineQuad Line2BoundaryElement :: fei(1, 2);

Line2BoundaryElement :: Line2BoundaryElement(int n, Domain *aDomain) : FMElement(n, aDomain), SpatialLocalizerInterface(this)
{
    this->numberOfDofMans = 3;
    integrationRulesArray.resize( 0 );
}

Line2BoundaryElement :: ~Line2BoundaryElement()
{ }

FEInterpolation *Line2BoundaryElement :: giveInterpolation() const
{
    return & this->fei;
}

void Line2BoundaryElement :: giveDofManDofIDMask(int i, IntArray &nodeDofIDMask) const
{
    nodeDofIDMask = {V_u, V_v};
}

Interface *Line2BoundaryElement :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case SpatialLocalizerInterfaceType:
        return static_cast< SpatialLocalizerInterface * >(this);

    case EIPrimaryUnknownMapperInterfaceType:
        return static_cast< EIPrimaryUnknownMapperInterface * >(this);

    default:
        return FMElement :: giveInterface(it);
    }
}

void Line2BoundaryElement :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
                                                                                   TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n;
    this->fei.evalN( answer, lcoords, FEIElementGeometryWrapper(this) );

    IntArray dofIDs;
    this->giveElementDofIDMask(dofIDs);

    answer.resize( dofIDs.giveSize() );
    answer.zero();
    for ( int i = 1; i <= n.giveSize(); ++i ) {
        for ( int j = 1; j <= dofIDs.giveSize(); ++j ) {
            answer.at(j) += n.at(i) * this->giveNode(i)->giveDofWithID( dofIDs.at(j) )->giveUnknown(mode, tStep);
        }
    }
}
} // end namespace oofem
