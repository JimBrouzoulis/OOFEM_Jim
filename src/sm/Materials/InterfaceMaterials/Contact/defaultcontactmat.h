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

#ifndef defaultcontactmat_h
#define defaultcontactmat_h

#include "material.h"
#include "Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"

///@name Input fields for DefaultContactMat
//@{
#define _IFT_DefaultContactMat_Name "defaultcontactmat"
#define _IFT_DefaultContactMat_epsN "epsn"

//@}
//
///@author Jim Brouzoulis
//
namespace oofem {


/**
 * This class is the default constitutive behaviour for a contact interface. 
 * It implements a linear penalty spring in the normal direction only such 
 * that the normal pressure pN = epsN * gN. 
 */
class DefaultContactMat : public StructuralInterfaceMaterial
{
protected:
    double epsN;
    
public:
    /// Constructor
    DefaultContactMat( int n, Domain *d );
    /// Destructor
    virtual ~DefaultContactMat(){};

    virtual int hasNonLinearBehaviour() { return 0; }
    virtual const char *giveInputRecordName() const { return _IFT_DefaultContactMat_Name; }

    virtual void giveEngTraction_3d( FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &jump, TimeStep *);

    
    virtual void give3dStiffnessMatrix_Eng( FloatMatrix &answer, MatResponseMode rMode,
                                            GaussPoint *gp, TimeStep *tStep );

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual bool hasAnalyticalTangentStiffness( ) const { return true; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new StructuralInterfaceMaterialStatus(1, domain, gp); }
};

} // end namespace oofem
#endif // defaultcontactmat_h
