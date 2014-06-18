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

#include "Shell7BasePhFi.h"
#include "node.h"
#include "load.h"
#include "structuralms.h"
#include "mathfem.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "feinterpol3d.h"
#include "fei3dtrquad.h"
#include "boundaryload.h"
#include "constantpressureload.h"
#include "constantsurfaceload.h"
#include "vtkxmlexportmodule.h"
#include "fracturemanager.h"
#include <fstream>

namespace oofem {

const int nLayers = 4;

Shell7BasePhFi :: Shell7BasePhFi(int n, Domain *aDomain) : Shell7Base(n, aDomain), PhaseFieldElement(n, aDomain){
	this->numberOfLayers = nLayers;
}

IRResultType Shell7BasePhFi :: initializeFrom(InputRecord *ir)
{
    Shell7Base :: initializeFrom(ir);
    return IRRT_OK;
}



void
Shell7BasePhFi :: postInitialize()
{
    
//    this->layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection()  );
//    this->fei       = dynamic_cast< FEInterpolation3d   * >( this->giveInterpolation() );
//    this->setupInitialNodeDirectors();

//    Element :: postInitialize();
    Shell7Base :: postInitialize();

	this->startIDdamage = this->domain->giveNextFreeDofID();
	
	if (!this->layeredCS->giveNumberOfLayers() == NULL) {
		this->endIDdamage = this->startIDdamage + this->layeredCS->giveNumberOfLayers() - 1;
	} else {
		this->endIDdamage = this->startIDdamage + numberOfLayers - 1;
	}
	
}

const IntArray &
Shell7BasePhFi :: giveOrdering(SolutionField fieldType) const
{
    OOFEM_ERROR("Shell7BasePhFi :: giveOrdering not implemented: Use Shell7BasePhFi :: giveOrderingPhFi instead");
	return 0;
}

//Interface *Shell7BasePhFi :: giveInterface(InterfaceType it)
//{
//    switch ( it ) {
//    case NodalAveragingRecoveryModelInterfaceType:
//        return static_cast< NodalAveragingRecoveryModelInterface * >( this );
//
//    case LayeredCrossSectionInterfaceType:
//        return static_cast< LayeredCrossSectionInterface * >( this );
//
//    case VTKXMLExportModuleElementInterfaceType:
//        return static_cast< VTKXMLExportModuleElementInterface * >( this );
//
//    case ZZNodalRecoveryModelInterfaceType:
//        return static_cast< ZZNodalRecoveryModelInterface * >( this );
//
//    case FailureModuleElementInterfaceType:
//        return static_cast< FailureModuleElementInterface * >( this );
//
//    default:
//        return StructuralElement :: giveInterface(it);
//    }
//}


void
Shell7BasePhFi :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer)
{
	IntArray answer_u, answer_d;

	this->giveDofManDofIDMask_u(answer_u);

    giveDofManDofIDMask_u(answer_u);
	giveDofManDofIDMask_d(answer_d);

	answer.resize(0);
	answer.followedBy(answer_u);
	answer.followedBy(answer_d);
	
	answer_u.printYourself();
	answer_d.printYourself();
	answer.printYourself();

}

void
Shell7BasePhFi :: giveDofManDofIDMask_u(IntArray &answer)
{
	answer.setValues(7, D_u, D_v, D_w, W_u, W_v, W_w, Gamma);
}

void
Shell7BasePhFi :: giveDofManDofIDMask_d(IntArray &answer)
{
	int sID, eID;

	sID = this->startIDdamage;
	eID = this->startIDdamage;

	answer.setValues(1, sID);		//@todo: modify!!!!
	for (int i = 2; i < eID - sID + 2; i++)
	{
		answer.followedBy(sID + i -1);
	}
}


double 
Shell7BasePhFi :: computeDamageAt( GaussPoint *gp, ValueModeType valueMode,  TimeStep *stepN)
{
    // d = N_d * a_d
    //NLStructuralElement *el = this->giveElement( );
    FloatArray dVec, dVecRed;
	
    this->computeDamageUnknowns(dVec, valueMode, stepN);		// should be a vector with all damage nodal values for this element
																// ordered such that all damage dofs associated with node 1 comes first
    FloatArray Nvec, lcoords;
	//int layer = this->layeredCS->giveLayer(gp);					
	int layer = this->layeredCS->give(CS_Layer, gp);					// @todo: possibly unnecessary expensive to each time loop over
	int numberOfNodes = this->giveNumberOfDofManagers();		// layers to determine the actual layer in order to extract the correct dofs
	
	IntArray indx(numberOfNodes);

	for (int i = 1; i <= numberOfNodes; i++)
	{
		indx.at(i) = layer + (i-1)*numberOfNodes;
	}
	indx.printYourself();

	dVecRed.beSubArrayOf(dVec, indx);
    this->giveInterpolation()->evalN(Nvec, *gp->giveCoordinates(), FEIElementGeometryWrapper(this));		// @todo anything strange here!!!!???
    return Nvec.dotProduct(dVec);
}


double 
Shell7BasePhFi :: computeDamageInLayerAt(int layer, GaussPoint *gp, ValueModeType valueMode,  TimeStep *stepN)
{
    // d = N_d * a_d
    //NLStructuralElement *el = this->giveElement( );
    FloatArray dVec, dVecRed;
	
    this->computeDamageUnknowns(dVec, valueMode, stepN);		// should be a vector with all damage nodal values for this element
																// ordered such that all damage dofs associated with node 1 comes first
    FloatArray Nvec, lcoords;
	//int layer = this->layeredCS->giveLayer(gp);					
	IntArray indx = computeDamageIndexArray(layer);
	
	indx.printYourself();

	dVecRed.beSubArrayOf(dVec, indx);
    this->giveInterpolation()->evalN(Nvec, *gp->giveCoordinates(), FEIElementGeometryWrapper(this));		// @todo anything strange here!!!!???
    return Nvec.dotProduct(dVec);
}

IntArray
Shell7BasePhFi :: computeDamageIndexArray(int layer)
{
	int numberOfNodes = this->giveNumberOfDofManagers();		
	
	IntArray indx(numberOfNodes);

	for (int i = 1; i <= numberOfNodes; i++)
	{
		indx.at(i) = layer + (i-1)*numberOfNodes;
	}
	return indx;

}

double
Shell7BasePhFi  :: computeGInLayer(int layer, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // computes Dg/Dd = (1-d)^2 + r0
    double d = this->computeDamageInLayerAt(layer, gp, valueMode, stepN);
    double r0 = 1.0e-10;
    return (1.0 - d) * (1.0 - d) + r0;
}

double 
Shell7BasePhFi  :: computeGprimInLayer(int layer, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN)
{
    // compute -2*(1-d)
    double d = this->computeDamageInLayerAt(layer, gp, valueMode, stepN);
    return -2.0 * (1.0 - d);
}
int
Shell7BasePhFi :: giveNumberOfDofs() 
{
    return this->giveNumberOfuDofs() + this->giveNumberOfdDofs();
}

int
Shell7BasePhFi :: giveNumberOfuDofs() 
{
    return 7 * this->giveNumberOfDofManagers();
}

int
Shell7BasePhFi :: giveNumberOfdDofs() 
{
    return this->endIDdamage - this->startIDdamage + 1;
}


// Tangent matrices

void
Shell7BasePhFi :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
	Shell7Base :: computeStiffnessMatrix(answer, rMode, tStep);
}

#if 0
void
Shell7BasePhFi :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    int ndofs = this->giveNumberOfuDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    
    FloatArray solVec, totalSolVec;
    this->giveUpdatedSolutionVector(totalSolVec, tStep);  
    this->giveUpdatedSolutionVector(solVec, tStep); // a
    
    // first 'solVec' corresponds to the point where the tangent i evaluated and solVecLeft, solVecRight
    // corresponds to the solution used for evaluation of the lambda matrices
    ///@todo rewrite this method since the XFEM part does not use this anymore
    this->new_computeBulkTangentMatrix(answer, solVec, solVec, solVec, rMode, tStep);		//@todo: Modify to include damage!!!!


    // Add contribution due to pressure load ///@todo should later be compted by the load
    int nLoads = this->boundaryLoadArray.giveSize() / 2;

    for ( int i = 1; i <= nLoads; i++ ) {     // For each pressure load that is applied
        int load_number = this->boundaryLoadArray.at(2 * i - 1);
        int iSurf = this->boundaryLoadArray.at(2 * i);         // load_id
        Load *load = this->domain->giveLoad(load_number);

        if ( dynamic_cast< ConstantPressureLoad * >( load ) ) {
            FloatMatrix K_pressure;
            this->computePressureTangentMatrix(K_pressure, load, iSurf, tStep);
            answer.add(K_pressure); // Should assemble with ordering
        }
    }

}
#endif

void
Shell7BasePhFi :: computeStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) {


}

void
Shell7BasePhFi :: computeStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) {


}

void
Shell7BasePhFi :: computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) {


}


void
Shell7BasePhFi :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) 
{
    //set displacement and nonlocal location array
    ///@todo this part is enough to do once
    //IntArray IdMask_u, IdMask_d;
    //this->giveDofManDofIDMask_u( IdMask_u );
    //this->giveDofManDofIDMask_d( IdMask_d );
    //this->computeLocationArrayOfDofIDs( IdMask_u, loc_u );		//@todo: needs to be defined??
    //this->computeLocationArrayOfDofIDs( IdMask_d, loc_d );		//@todo: needs to be defined??

	IntArray loc_u, loc_d;
	loc_u = giveOrderingPhFi(Displacement);
	loc_d = giveOrderingPhFi(Damage);

    int nDofs = this->computeNumberOfDofs();
    answer.resize( nDofs, nDofs );
    answer.zero();

    FloatMatrix answer1, answer2, answer3, answer4;
    this->computeStiffnessMatrix_uu(answer1, rMode, tStep);
    this->computeStiffnessMatrix_ud(answer2, rMode, tStep);
    this->computeStiffnessMatrix_du(answer3, rMode, tStep); //symmetric
    this->computeStiffnessMatrix_dd(answer4, rMode, tStep);
    
    answer.assemble( answer1, loc_u, loc_u );
    answer.assemble( answer2, loc_u, loc_d );
    answer.assemble( answer3, loc_d, loc_u );
    answer.assemble( answer4, loc_d, loc_d );
}


void
Shell7BasePhFi :: new_computeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, FloatArray &solVecI, FloatArray &solVecJ, MatResponseMode rMode, TimeStep *tStep)
{
    //@todo optimize method - remove I,J-stuff since XFEM does not need this anymore
    FloatMatrix A [ 3 ] [ 3 ], lambdaI [ 3 ], lambdaJ [ 3 ];
    FloatMatrix L(18,18), B;
    FloatMatrix K, tempAnswer;
	double g;

    int ndofs = Shell7BasePhFi :: giveNumberOfuDofs();
    answer.resize(ndofs, ndofs); tempAnswer.resize(ndofs, ndofs);
    answer.zero(); tempAnswer.zero();

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     
    FloatMatrix temp;
    FloatArray genEpsI, genEpsJ, genEps, lCoords;

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
        StructuralMaterial *mat = static_cast< StructuralMaterial* >( domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) ) );

        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);
            lCoords = *gp->giveCoordinates();

            this->computeBmatrixAt(lCoords, B, 0, 0);
            this->computeGeneralizedStrainVectorNew(genEpsI, solVecI, B);
            this->computeGeneralizedStrainVectorNew(genEpsJ, solVecJ, B);
            this->computeGeneralizedStrainVectorNew(genEps , solVec , B);

            // Material stiffness
            Shell7Base :: computeLinearizedStiffness(gp, mat, tStep, A, genEps);

            //double zeta = giveGlobalZcoord(gp->giveCoordinate(3));
            //double zeta = giveGlobalZcoord( gp->giveCoordinate( 3 ), *gp->giveCoordinates() );
            double zeta = giveGlobalZcoord( *gp->giveCoordinates() );
            this->computeLambdaGMatrices(lambdaI, genEpsI, zeta);
            this->computeLambdaGMatrices(lambdaJ, genEpsJ, zeta);

            // L = sum_{i,j} (lambdaI_i)^T * A^ij * lambdaJ_j
            // @todo Naive implementation - should be optimized 
            // note: L will only be symmetric if lambdaI = lambdaJ
            L.zero();
            for ( int j = 0; j < 3; j++ ) {
                for ( int k = 0; k < 3; k++ ) {
                    this->computeTripleProduct(temp, lambdaI [ j ], A [ j ][ k ], lambdaJ [ k ]);
                    L.add(temp);
                }
            }
            
            this->computeTripleProduct(K, B, L, B);
            double dV = this->computeVolumeAroundLayer(gp, layer);
			//damage = this->computeDamageAt(gp, VM_Total,  tStep);
			//g = this->computeGInLayer(gp, VM_Total,  tStep);				// Check so that the correct computeDamageAt method is used (in Shell7BasePhFi)
			g = this->computeGInLayer(layer, gp, VM_Total,  tStep);				
			K.times(g);

			//@todo: add scaling of stiffness matrix by damage!!!!
            tempAnswer.add(dV, K);


        }
    }

    const IntArray &ordering = this->giveOrderingPhFi(Displacement);		// @todo: defined in the actual element!!!!
    answer.assemble(tempAnswer, ordering, ordering);

}

#if 0
int 
	Shell7BasePhFi :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
	// Compute special IST quantities this element should support
	switch (type) {
	case IST_CauchyStressTensor:
		this->computeCauchyStressVector(answer, gp, tStep);			//@todo: add damage!!!
		return 1;
	default:
		return Element :: giveIPValue(answer, gp, type, tStep);
	}
}  
#endif // 0


// Internal forces

#if 0
void
Shell7BasePhFi :: giveUpdatedSolutionVector_d(FloatArray &answer, TimeStep *tStep)
{

	IntArray dofIdArray;
    
	//Shell7Base :: giveDofManDofIDMask(dummy, EID_MomentumBalance, dofIdArray);
    this->giveDofManDofIDMask_d(dofIdArray);

	this->computeVectorOfDofIDs(dofIdArray, VM_Total, tStep, temp);
    answer.assemble( temp, this->giveOrdering(AllInv) );
}
#endif

void
Shell7BasePhFi :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes internal forces as a summation of: sectional forces + convective mass force

    FloatArray solVec_u, solVec_d, solVec;
    this->giveUpdatedSolutionVector(solVec, tStep); // placement vector x
	//this->giveUpdatedSolutionVector_d(solVec_d, tStep); // damage vector, one component per layer
	
	//solVec.resize(solVec_u.giveSize() + solVec_d.giveSize());

    this->computeSectionalForces(answer, tStep, solVec, useUpdatedGpRecord);

    ///@todo How to treat the convective force? Only active during dynamic simulations
}


void
Shell7BasePhFi :: computeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord)
{
    //
    int ndofs = Shell7BasePhFi :: giveNumberOfDofs();
	int ndofs_u = Shell7BasePhFi :: giveNumberOfuDofs();
	int ndofs_d = ndofs - ndofs_u;
	double g;

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    int sectionalForces_ds;
	FloatArray fu(ndofs_u), fd(ndofs_d), ftemp_u, ftemp_d, sectionalForces, sectionalForces_dv;
    FloatArray genEps, genEpsD, totalSolVec, lCoords, Nd;
    FloatMatrix B, Bd;
    this->giveUpdatedSolutionVector(totalSolVec, tStep);    // => x, m, gam

    ///@todo For multiscale: This integration should be changed to a surface integration and 'sectionalForces' should be given by the RVE problem 

    fu.zero();
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );
		IntArray indx_d = computeDamageIndexArray(layer);

        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);
            lCoords = *gp->giveCoordinates();
            this->computeBmatrixAt(lCoords, B);
			this->computeBdmatrixAt(lCoords, Bd);
			this->computeNdvectorAt(lCoords, ftemp_d);


            this->computeGeneralizedStrainVectorNew(genEpsD, solVec, B);
            this->computeGeneralizedStrainVectorNew(genEps, totalSolVec, B); // used for computing the stress

            double zeta = giveGlobalZcoord( *gp->giveCoordinates() ); 
            this->computeSectionalForcesAt(sectionalForces, gp, mat, tStep, genEps, genEpsD, zeta); // these are per unit volume
			this->computeSectionalForcesAt_d(sectionalForces_ds, sectionalForces_dv, gp, mat, tStep, zeta, layer); // these are per unit volume

				
            // Computation of sectional forces: f = B^t*[N M T Ms Ts]^t
            ftemp_u.beTProductOf(B,sectionalForces);
            double dV = this->computeVolumeAroundLayer(gp, layer);
			g = this->computeGInLayer(layer, gp, VM_Total,  tStep);
			OOFEM_ERROR("Shell7BasePhFi :: computeSectionalForces, Need to implement sectional forces for damage evolution")
			fu.add(dV*g, ftemp_u);

			ftemp_d.times(sectionalForces_ds*dV);
			ftemp_d.plusProduct(Bd, sectionalForces_dv, dV);
			
			fd.assemble(ftemp_d, indx_d);

        }
    }

    answer.resize( ndofs );
    answer.zero();
    const IntArray &ordering_disp = this->giveOrderingPhFi(Displacement);
	const IntArray &ordering_damage = this->giveOrderingPhFi(Damage);
    answer.assemble(fu, ordering_disp);		//ordering_disp contains only displacement related dofs, not damage
	answer.assemble(fd, ordering_damage);

}



void
Shell7BasePhFi :: computeSectionalForcesAt_d(int sectionalForcesScal, FloatArray &sectionalForcesVec, IntegrationPoint *ip, Material *mat, TimeStep *tStep, double zeta, int layer)
{
}

// Computation of solution vectors


void
Shell7BasePhFi :: computeVectorOfDofIDs(const IntArray &dofIdArray, ValueModeType u, TimeStep *stepN, FloatArray &answer)
{
    // Routine to extract the solution vector for an element given an dofid array.
    // Size will be numberOfDofs and if a certain dofId does not exist a zero is used as value. 
    // This method may e.g. be used to obtain the enriched part of the solution vector
    ///@todo: generalize so it can be used by all XFEM elements
    answer.resize( Shell7BasePhFi ::giveNumberOfDofs());
    answer.zero();
    int k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);        
        for (int j = 1; j <= dofIdArray.giveSize(); j++ ) {
            
            if ( dMan->hasDofID( (DofIDItem) dofIdArray.at(j) ) ) {
                Dof *d = dMan->giveDofWithID( dofIdArray.at(j) );
                ///@todo: will fail if any other dof then gamma is excluded from enrichment 
                /// since I only add "j". Instead I should skip certain dof numbers when incrementing
                answer.at(k+j) = d->giveUnknown(u, stepN);	// Martin: modification to be used also for rates
            }
        }
        k += 7;
    }
}

// N and B matrices

void
computeBdmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui) //@todo: define for damage
{
	answer.resize(1,1);
	answer.zero();
}

void
computeBdmatrixAt(FloatArray &lCoords, FloatMatrix &answer, int li, int ui) //@todo: define for damage
{
	answer.resize(1,1);
	answer.zero();

} 

void
computeNdvectorAt(const FloatArray &iLocCoords, FloatArray &answer) //@todo: define for damage
{
	answer.resize(1);
	answer.zero();

}

#if 0
void
Shell7BasePhFi :: computeBmatrixAt(FloatArray &lcoords, FloatMatrix &answer, int li, int ui)
{
    // Returns the  matrix {B} of the receiver, evaluated at aGaussPoint. Such that
    // B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns
 
    int ndofs = Shell7BasePhFi :: giveNumberOfDofs();
    int ndofs_xm  = 3 * this->giveNumberOfDofManagers();
    answer.resize(18, ndofs);
    answer.zero();
    FloatArray N;
    FloatMatrix dNdxi;
    this->fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );
    this->fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

    /*    18   18   6
     * 6 [B_u   0   0
     * 6   0   B_w  0
     * 3   0   N_w  0
     * 2   0    0  B_gam
     * 1   0    0  N_gam]
     */
    int ndofman = this->giveNumberOfDofManagers();

    // First column
    for ( int i = 1, j = 0; i <= ndofman; i++, j += 3 ) {
        answer.at(1, 1 + j) = dNdxi.at(i, 1);
        answer.at(2, 2 + j) = dNdxi.at(i, 1);
        answer.at(3, 3 + j) = dNdxi.at(i, 1);
        answer.at(4, 1 + j) = dNdxi.at(i, 2);
        answer.at(5, 2 + j) = dNdxi.at(i, 2);
        answer.at(6, 3 + j) = dNdxi.at(i, 2);
    }

    // Second column
    for ( int i = 1, j = 0; i <= ndofman; i++, j += 3 ) {
        answer.at(7, ndofs_xm + 1 + j) = dNdxi.at(i, 1);
        answer.at(8, ndofs_xm + 2 + j) = dNdxi.at(i, 1);
        answer.at(9, ndofs_xm + 3 + j) = dNdxi.at(i, 1);
        answer.at(10, ndofs_xm + 1 + j) = dNdxi.at(i, 2);
        answer.at(11, ndofs_xm + 2 + j) = dNdxi.at(i, 2);
        answer.at(12, ndofs_xm + 3 + j) = dNdxi.at(i, 2);
        answer.at(13, ndofs_xm + 1 + j) = N.at(i);
        answer.at(14, ndofs_xm + 2 + j) = N.at(i);
        answer.at(15, ndofs_xm + 3 + j) = N.at(i);
    }

    // Third column
    for ( int i = 1, j = 0; i <= ndofman; i++, j += 1 ) {
        answer.at(16, ndofs_xm * 2 + 1 + j) = dNdxi.at(i, 1);
        answer.at(17, ndofs_xm * 2 + 1 + j) = dNdxi.at(i, 2);
        answer.at(18, ndofs_xm * 2 + 1 + j) = N.at(i);
    }
}


void
Shell7BasePhFi :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
{
    // Returns the displacement interpolation matrix {N} of the receiver,
    // evaluated at aGaussPoint.

    int ndofs = Shell7BasePhFi :: giveNumberOfDofs();
    int ndofs_xm  = 3 * this->giveNumberOfDofManagers();
    answer.resize(7, ndofs);
    answer.zero();
    FloatArray N;
    this->fei->evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );
    
    /*   nno*3 nno*3 nno
     * 3 [N_x   0    0
     * 3   0   N_m   0
     * 1   0    0  N_gmm ]
     */
    for ( int i = 1, j = 0; i <= this->giveNumberOfDofManagers(); i++, j += 3 ) {
        answer.at(1, 1 + j) = N.at(i);
        answer.at(2, 2 + j) = N.at(i);
        answer.at(3, 3 + j) = N.at(i);
        answer.at(4, ndofs_xm + 1 + j) = N.at(i);
        answer.at(5, ndofs_xm + 2 + j) = N.at(i);
        answer.at(6, ndofs_xm + 3 + j) = N.at(i);
        answer.at(7, ndofs_xm * 2 + i) = N.at(i);
    }
}

#endif



// VTK export
#if 0
void
Shell7BasePhFi :: vtkEvalInitialGlobalCoordinateAt(FloatArray &localCoords, int layer, FloatArray &globalCoords)
{
    double zeta = giveGlobalZcoordInLayer(localCoords.at(3), layer);
    FloatArray N;
    this->fei->evalN( N, localCoords, FEIElementGeometryWrapper(this) );

    globalCoords.resize(3);
    globalCoords.zero();
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        FloatArray &xbar = *this->giveNode(i)->giveCoordinates();
        FloatArray M = this->giveInitialNodeDirector(i);
        globalCoords += N.at(i) * ( xbar + zeta * M );
    }

}

void
Shell7BasePhFi :: vtkEvalInitialGlobalCZCoordinateAt(FloatArray &localCoords, int interface, FloatArray &globalCoords)
{
    double zeta = giveGlobalZcoordInLayer(1.0, interface);
    FloatArray N;
    this->fei->evalN( N, localCoords, FEIElementGeometryWrapper(this) );

    globalCoords.resize(3);
    globalCoords.zero();
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        FloatArray &xbar = *this->giveNode(i)->giveCoordinates();
        FloatArray M = this->giveInitialNodeDirector(i);
        globalCoords += N.at(i) * ( xbar + zeta * M );
    }

}

void
Shell7BasePhFi :: vtkEvalUpdatedGlobalCoordinateAt(FloatArray &localCoords, int layer, FloatArray &globalCoords, TimeStep *tStep)
{
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep);
    FloatArray x, m; double gam=0;
    this->giveUnknownsAt(localCoords, solVec, x, m, gam, tStep); 
    double zeta = giveGlobalZcoordInLayer(localCoords.at(3), layer);
    double fac = ( zeta + 0.5 * gam * zeta * zeta );
    globalCoords = x + fac*m;
}

void
Shell7BasePhFi :: giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{
    vtkPieces.resize(1);
    this->giveShellExportData(vtkPieces[0], primaryVarsToExport, internalVarsToExport, cellVarsToExport, tStep );
    //this->giveCZExportData(vtkPieces[1], primaryVarsToExport, internalVarsToExport, cellVarsToExport, tStep );
    
}

void 
Shell7BasePhFi :: giveShellExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )            
{   

    int numCells = this->layeredCS->giveNumberOfLayers();
    const int numCellNodes  = 15; // quadratic wedge
    int numTotalNodes = numCellNodes*numCells;

    vtkPiece.setNumberOfCells(numCells);
    vtkPiece.setNumberOfNodes(numTotalNodes);

    std::vector <FloatArray> nodeCoords;
    int val    = 1;
    int offset = 0;
    IntArray nodes(numCellNodes);

    // Compute fictious node coords
    int nodeNum = 1;
    for ( int layer = 1; layer <= numCells; layer++ ) {
        

        // Node coordinates
        this->giveFictiousNodeCoordsForExport(nodeCoords, layer);       
        
        for ( int node = 1; node <= numCellNodes; node++ ) {    
            vtkPiece.setNodeCoords(nodeNum, nodeCoords[node-1] );
            nodeNum += 1;
        }

        // Connectivity       
        for ( int i = 1; i <= numCellNodes; i++ ) {            
            nodes.at(i) = val++;
        }
        vtkPiece.setConnectivity(layer, nodes);
        
        // Offset
        offset += numCellNodes;
        vtkPiece.setOffset(layer, offset);

        // Cell types
        vtkPiece.setCellType(layer, 26); // Quadratic wedge
    }


    // Export nodal variables from primary fields        
    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport.giveSize(), numTotalNodes);

    std::vector<FloatArray> updatedNodeCoords;
    FloatArray u(3);
    std::vector<FloatArray> values;
    for ( int fieldNum = 1; fieldNum <= primaryVarsToExport.giveSize(); fieldNum++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(fieldNum);
        nodeNum = 1;
        for ( int layer = 1; layer <= numCells; layer++ ) {            
            
            if ( type == DisplacementVector ) { // compute displacement as u = x - X
                this->giveFictiousNodeCoordsForExport(nodeCoords, layer);
                this->giveFictiousUpdatedNodeCoordsForExport(updatedNodeCoords, layer, tStep);
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    u = updatedNodeCoords[j-1];
                    u.subtract(nodeCoords[j-1]);
                    vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, u);
                    nodeNum += 1;        
                }

            } else {
                NodalRecoveryMI_recoverValues(values, layer, ( InternalStateType ) 1, tStep); // does not work well - fix
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, values[j-1]);
                    nodeNum += 1;
                }
            }
        }
    }

    // Export nodal variables from internal fields
    
    vtkPiece.setNumberOfInternalVarsToExport( internalVarsToExport.giveSize(), numTotalNodes );
    for ( int fieldNum = 1; fieldNum <= internalVarsToExport.giveSize(); fieldNum++ ) {
        InternalStateType type = ( InternalStateType ) internalVarsToExport.at(fieldNum);
        nodeNum = 1;
        //this->recoverShearStress(tStep);
        for ( int layer = 1; layer <= numCells; layer++ ) {            
            recoverValuesFromIP(values, layer, type, tStep);        
            for ( int j = 1; j <= numCellNodes; j++ ) {
                vtkPiece.setInternalVarInNode( fieldNum, nodeNum, values[j-1] );
                //ZZNodalRecoveryMI_recoverValues(el.nodeVars[fieldNum], layer, type, tStep);          
                nodeNum += 1;        
            }                                
        }  
    }


    // Export cell variables
    FloatArray average;
    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), numCells);
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);;
      
        for ( int layer = 1; layer <= numCells; layer++ ) {     
            IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
            VTKXMLExportModule::computeIPAverage(average, iRuleL, this, type, tStep);
            
            if ( average.giveSize() == 6 ) {
                vtkPiece.setCellVar(i, layer, convV6ToV9Stress(average) );
            } else {
                vtkPiece.setCellVar(i, layer, average );
            }

        }

    }




}


void 
Shell7BasePhFi :: recoverValuesFromIP(std::vector<FloatArray> &recoveredValues, int layer, InternalStateType type, TimeStep *tStep)
{
    // recover nodal values by coosing the ip closest to the node

    //FEInterpolation *interpol = static_cast< FEInterpolation * >( &this->interpolationForExport );

    // composite element interpolator
    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);

    int numNodes = localNodeCoords.giveNumberOfColumns();
    recoveredValues.resize(numNodes);
    
    IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
    IntegrationPoint *ip;

    // Find closest ip to the nodes
    IntArray closestIPArray(numNodes);
    FloatArray nodeCoords, ipCoords, ipValues;

    for ( int i = 1; i <= numNodes; i++ ) {
        nodeCoords.beColumnOf(localNodeCoords, i);
        double distOld = 3.0; // should not be larger
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            ip = iRule->getIntegrationPoint(j);
            ipCoords = *ip->giveCoordinates();
            double dist = nodeCoords.distance(ipCoords);
            if ( dist < distOld ) {
                closestIPArray.at(i) = j;
                distOld = dist;
            }
        }
    }

   InternalStateValueType valueType =  giveInternalStateValueType(type);

    // recover ip values
    for ( int i = 1; i <= numNodes; i++ ) {
        ip = iRule->getIntegrationPoint( closestIPArray.at(i) );
        this->giveIPValue(ipValues, ip, type, tStep);
        if ( valueType == ISVT_TENSOR_S3 ) {
            recoveredValues[i-1].resize(9);
            recoveredValues[i-1] = convV6ToV9Stress(ipValues);
        } else if ( ipValues.giveSize() == 0 && type == IST_AbaqusStateVector) {
            recoveredValues[i-1].resize(23);
            recoveredValues[i-1].zero();
        } else if ( ipValues.giveSize() == 0 ) {
            recoveredValues[i-1].resize(giveInternalStateTypeSize(valueType));
            recoveredValues[i-1].zero();

        } else {
            recoveredValues[i-1] = ipValues;
        }
    }

}


void 
Shell7BasePhFi :: recoverShearStress(TimeStep *tStep)
{
    // Recover shear stresses at ip by numerical integration of the momentum balance through the thickness
    std::vector<FloatArray> recoveredValues;
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    IntegrationRule *iRuleThickness = specialIntegrationRulesArray[ 0 ];
    FloatArray dS, Sold;
    FloatMatrix B, Smat(2,6); // 2 stress components * num of in plane ip ///@todo generalize
    Smat.zero();
    FloatArray Tcon(6), Trec(6);  Tcon.zero(); Trec.zero();
    
     for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        this->recoverValuesFromIP(recoveredValues, layer, IST_StressTensor, tStep);
        //this->ZZNodalRecoveryMI_recoverValues(recoveredValues, layer, IST_StressTensor, tStep);
        double thickness = this->layeredCS->giveLayerThickness(layer);

        //set up vector of stresses in the ip's = [S1_xx, S1_yy, S1_xy, ..., Sn_xx, Sn_yy, Sn_xy]
        int numNodes = 15;
        FloatArray aS(numNodes*3); 
        for ( int j = 1, pos = 0; j <= numNodes; j++, pos+=3 ) {
            aS.at(pos + 1) = recoveredValues[j-1].at(1);   // S_xx
            aS.at(pos + 2) = recoveredValues[j-1].at(2);   // S_yy
            aS.at(pos + 3) = recoveredValues[j-1].at(6);   // S_xy
        }
        int numInPlaneIP = 6;

        for ( int i = 0; i < iRuleThickness->giveNumberOfIntegrationPoints(); i++ ) { 
            double  dz = thickness * iRuleThickness->getIntegrationPoint(i)->giveWeight();
            
            for ( int j = 0; j < numInPlaneIP; j++ ) { 

                int point = i*numInPlaneIP + j; // integration point number
                GaussPoint *gp = iRuleL->getIntegrationPoint(point);

                this->computeBmatrixForStressRecAt(*gp->giveCoordinates(), B, layer);
                dS.beProductOf(B,aS*(-dz)); // stress increment

                StructuralMaterialStatus* status = dynamic_cast< StructuralMaterialStatus* > ( gp->giveMaterialStatus() );
                Sold = status->giveStressVector();
                
                Smat.at(1,j+1) += dS.at(1); // add increment from each level
                Smat.at(2,j+1) += dS.at(2);

                //Tcon.at(j+1) += Sold.at(5)*dz;

                // Replace old stresses with  - this should probably not be done as it may affect the convergence in a nonlinear case
                Sold.at(5) = Smat.at(1,j+1); // S_xz
                Sold.at(4) = Smat.at(2,j+1); // S_yz


                status->letStressVectorBe(Sold);
                //Trec.at(j+1) += Sold.at(5)*dz;
            }
        }


    }

}


void
Shell7BasePhFi :: computeBmatrixForStressRecAt(FloatArray &lcoords, FloatMatrix &answer, int layer)
{
    // Returns the  special matrix {B} of the receiver, evaluated at aGaussPoint. Such that
    // B*a = [dS_xx/dx + dS_xy/dy, dS_yx/dx + dS_yy/dy ]^T, where a is the vector of in plane 
    // stresses [S_xx, S_yy, S_xy]
 
    // set up virtual cell geometry for an qwedge
    const int numNodes = 15;
    std::vector<FloatArray> nodes;
    giveFictiousNodeCoordsForExport(nodes, layer);

    int VTKWedge2EL [] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    FloatArray *coords[numNodes];


    for ( int i = 1; i <= numNodes; i++ ) {
        int pos = VTKWedge2EL[ i-1 ];
        coords[ i - 1 ] = &nodes[ pos - 1];
        
    }
    
    FEInterpolation *interpol = static_cast< FEInterpolation * >( &this->interpolationForExport );
    FloatMatrix dNdx;
    interpol->evaldNdx( dNdx, lcoords, FEIVertexListGeometryWrapper(numNodes, (const FloatArray **)coords ) );
    
    /*    
     * 1 [d/dx  0   d/dy
     * 1   0   d/dy d/dx]
     */
    int ndofs = numNodes*3;
    answer.resize(2, ndofs);
    for ( int i = 1, j = 0; i <= numNodes; i++, j += 3 ) {
        answer.at(1, j + 1) = dNdx.at(i, 1);
        answer.at(1, j + 3) = dNdx.at(i, 2);
        answer.at(2, j + 2) = dNdx.at(i, 2);
        answer.at(2, j + 3) = dNdx.at(i, 1);
    }
    
}





void 
Shell7BasePhFi :: giveFictiousNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer)
{
    // compute fictious node coords
    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
    
    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        this->vtkEvalInitialGlobalCoordinateAt(localCoords, layer, coords);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }

}


void 
Shell7BasePhFi :: giveFictiousCZNodeCoordsForExport(std::vector<FloatArray> &nodes, int interface)
{
    // compute fictious node coords
    FloatMatrix localNodeCoords;
    this->interpolationForCZExport.giveLocalNodeCoords(localNodeCoords);
    
    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        localCoords.at(3) = 1.0;
        this->vtkEvalInitialGlobalCoordinateAt(localCoords, interface, coords);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }

}

void 
Shell7BasePhFi :: giveFictiousUpdatedNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, TimeStep *tStep)
{
    // compute fictious node coords

    FloatMatrix localNodeCoords;
    this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, layer, coords, tStep);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }
}


//void
//Shell7BasePhFi :: giveLocalNodeCoordsForExport(FloatArray &nodeLocalXi1Coords, FloatArray &nodeLocalXi2Coords, FloatArray &nodeLocalXi3Coords) {
//    // Local coords for a quadratic wedge element (VTK cell type 26)
//    double z = 0.999;
//    nodeLocalXi1Coords.setValues(15, 1., 0., 0., 1., 0., 0., .5, 0., .5, .5, 0., .5, 1., 0., 0.);      
//    nodeLocalXi2Coords.setValues(15, 0., 1., 0., 0., 1., 0., .5, .5, 0., .5, .5, 0., 0., 1., 0.);
//    nodeLocalXi3Coords.setValues(15, -z, -z, -z,  z,  z,  z, -z, -z, -z,  z,  z,  z, 0., 0., 0.);
//}

#endif



// Misc functions



} // end namespace oofem
