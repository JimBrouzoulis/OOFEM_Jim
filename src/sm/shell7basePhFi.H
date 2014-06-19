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

#ifndef Shell7BasePhFi_h
#define Shell7BasePhFi_h

#include "eleminterpmapperinterface.h"
#include "nodalaveragingrecoverymodel.h"
#include "layeredcrosssection.h"
#include "nlstructuralelement.h"
#include "vtkxmlexportmodule.h"
#include "zznodalrecoverymodel.h"
#include "fei3dwedgequad.h"
#include "fei3dtrquad.h"
#include "fracturemanager.h"
#include "cltypes.h"
#include <vector>
#include "phasefieldelement.h"
#include "shell7base.h"


namespace oofem {
class BoundaryLoad;

#define _ExportCZ

/**
 * This class represent a 7 parameter shell element with isotropic damage in each layer modelled by phase field theory.
 * Each node has 7 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 * @todo Add ref. to paper!
 * @author Martin Fagerström
 * @date 2014-06-03
 */
class Shell7BasePhFi : public Shell7Base, public PhaseFieldElement
{
public:
	int startIDdamage, endIDdamage;

    Shell7BasePhFi(int n, Domain *d); // constructor
    virtual ~Shell7BasePhFi() {}
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    virtual void giveDofManDofIDMask_u(IntArray &answer) const;
    virtual void giveDofManDofIDMask_d(IntArray &answer) const;
	//virtual void giveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep);
    
    double computeDamageAt(GaussPoint *gp, ValueModeType valueMode,  TimeStep *stepN);

	double computeGInLayer(int layer, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
	double computeGprimInLayer(int layer, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
	double computeDamageInLayerAt(int layer, GaussPoint *gp, ValueModeType valueMode,  TimeStep *stepN);
	IntArray computeDamageIndexArray(int layer);

	virtual int computeNumberOfDofs() { return this->giveNumberOfDofs(); }
    //virtual int checkConsistency();
    virtual void postInitialize();

    // Definition & identification
    virtual const char *giveClassName() const { return "Shell7BasePhFi"; }
   // virtual classType giveClassID() const { return Shell7BasePhFiClass; }


    // Element specific
    virtual int giveNumberOfDofs();		//@todo: should be protected?
	virtual int giveNumberOfuDofs();	//@todo: should be protected?
	virtual int giveNumberOfdDofs();	//@todo: should be protected?
    //virtual int giveNumberOfEdgeDofs() = 0; //todo: used for PhFi???
    //virtual int giveNumberOfEdgeDofManagers() = 0;

	virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

protected:

	enum SolutionFieldPhFi {
        Midplane,       ///< phi_bar 7 x_bar (3 dofs)
        Director,       ///< m (3 dofs)
        InhomStrain,    ///< gamma (1 dofs) - inhomogenious thickness strain
		Displacement,	///< phi_bar, m, gamma
		Damage,			///< d (n dofs where n is the number of layers in the layeredCS
        All,
        AllInv,
        EdgeInv,
    };

	int numberOfLayers;

    virtual const IntArray &giveOrdering(SolutionField fieldType) const;
	virtual const IntArray &giveOrderingPhFi(SolutionFieldPhFi fieldType) const = 0;

    virtual IRResultType initializeFrom(InputRecord *ir);

    // Tangent matrices
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void new_computeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, FloatArray &solVecI, FloatArray &solVecJ, MatResponseMode rMode, TimeStep *tStep);

    void computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode, TimeStep *);

    // Internal forces
	void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void computeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord = 0);
	void computeSectionalForcesAt_d(int ans_scalar, FloatArray &answer, IntegrationPoint *ip, Material *mat, TimeStep *tStep,  double zeta, int layer);


    // External forces

    // Transformation
	
	// Solution vectors
    //void giveSolutionVector_u(FloatArray &answer, const IntArray &dofIdArray, TimeStep *tStep);
    void giveSolutionVector_d(FloatArray &answer, const IntArray &dofIdArray, TimeStep *tStep);
    void computeVectorOfDofIDs(const IntArray &dofIdArray, ValueModeType u, TimeStep *stepN, FloatArray &answer);
//    void temp_computeBoundaryVectorOf(IntArray &dofIdArray, int boundary, ValueModeType u, TimeStep *stepN, FloatArray &answer);

//    virtual void edgeGiveUpdatedSolutionVector(FloatArray &answer, const int iedge, TimeStep *tStep);
//    void edgeGiveInitialSolutionVector(FloatArray &answer, const int iedge);

    //void giveInitialSolutionVector(FloatArray &answer);
    //void giveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep);
    void giveUpdatedSolutionVector_d(FloatArray &answer, TimeStep *tStep);
    //void giveUnknownsAt(FloatArray &lcoords, FloatArray &solVec, FloatArray &x, FloatArray &m, double gam, TimeStep *tStep);

    // Nodal averaging interface:
    //virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
    //virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep);

    // ZZ recovery
    //virtual void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatArray &answer, GaussPoint *gp, InternalStateType type);
    //void NodalRecoveryMI_computeNValProduct(FloatMatrix &answer, int layer, InternalStateType type, TimeStep *tStep);
    //void NodalRecoveryMI_computeNNMatrix(FloatArray &answer, int layer, InternalStateType type);
    //void NodalRecoveryMI_recoverValues(std::vector<FloatArray> &recoveredValues, int layer, InternalStateType type, TimeStep *tStep);

    // VTK interface
        
    //virtual void giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep );
    //virtual void giveShellExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep );
    
    
//    void recoverValuesFromIP(std::vector<FloatArray> &nodes, int layer, InternalStateType type, TimeStep *tStep);
    //void computeBmatrixForStressRecAt(FloatArray &lcoords, FloatMatrix &answer, int layer);


    // N and B matrices
	virtual void computeBdmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS){}; //@todo: define for damage

	virtual void computeBdmatrixAt(FloatArray &lCoords, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS); //@todo: define for damage

	virtual void computeNdvectorAt(const FloatArray &iLocCoords, FloatArray &answer){}; //@todo: define for damage

    // Misc
    //virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};
} // end namespace oofem
#endif
