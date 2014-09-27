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

#ifndef regularizedcoulomb_h
#define regularizedcoulomb_h

#include "material.h"
#include "Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"

///@name Input fields for RegularizedCoulomb
//@{
#define _IFT_RegularizedCoulomb_Name "regularizedcoulomb"
#define _IFT_RegularizedCoulomb_mu "mu"
#define _IFT_RegularizedCoulomb_epsN "epsn"
#define _IFT_RegularizedCoulomb_epsT "epst"
#define _IFT_RegularizedCoulomb_epsReg "epsreg"

//@}
//
///@author Jim Brouzoulis
//
namespace oofem {
/**
 * This class implements the associated Material Status for RegularizedCoulomb.
 */
class RegularizedCoulombStatus : public StructuralInterfaceMaterialStatus
{
protected:
public:
    /// Constructor
    RegularizedCoulombStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~RegularizedCoulombStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();

};


/**
 * This class represents a regularized Coulomb friction model with a linear penalty formulation for the normal penetration. 
 * 
 */
class RegularizedCoulomb : public StructuralInterfaceMaterial
{
protected:
    double epsT;
    double epsN;
    double epsReg; // regularization parameter
    double mu;     // coefficient of friction
    
public:
    /// Constructor
    RegularizedCoulomb( int n, Domain *d );
    /// Destructor
    virtual ~RegularizedCoulomb();

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual const char *giveInputRecordName() const { return _IFT_RegularizedCoulomb_Name; }

    virtual void giveEngTraction_3d( FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &jump, TimeStep *);

    virtual void computeEngTraction_3d( FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &jump, TimeStep *);

    
    virtual void give3dStiffnessMatrix_Eng( FloatMatrix &answer, MatResponseMode rMode,
                                            GaussPoint *gp, TimeStep *tStep );
    

    void computeEngTraction(double &normalStress, FloatArray &shearStress, 
                             FloatArray &tempShearStressShift,
                             const double normalJump, const FloatArray &shearJump );

    double computePhi(const double gTdot);
    
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual bool hasAnalyticalTangentStiffness( ) const { return false; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new RegularizedCoulombStatus(1, domain, gp); }
};

} // end namespace oofem
#endif // regularizedcoulomb_h
