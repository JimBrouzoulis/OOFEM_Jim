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

#include "misesmat.h"
#include "Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "stressvector.h"
#include "strainvector.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(MisesMat);

// constructor
MisesMat :: MisesMat(int n, Domain *d) : StructuralMaterial(n, d)
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    H = 0.;
    sig0 = 0.;
    G = 0.;
    K = 0.;
}

// destructor
MisesMat :: ~MisesMat()
{
    delete linearElasticMaterial;
}

// specifies whether a given material mode is supported by this model
int
MisesMat :: hasMaterialModeCapability(MaterialMode mode)
{
    return mode == _3dMat || mode == _1dMat || mode == _PlaneStrain;
}

// reads the model parameters from the input file
IRResultType
MisesMat :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    StructuralMaterial :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants
    G = static_cast< IsotropicLinearElasticMaterial * >(linearElasticMaterial)->giveShearModulus();
    K = static_cast< IsotropicLinearElasticMaterial * >(linearElasticMaterial)->giveBulkModulus();

    IR_GIVE_FIELD(ir, sig0, _IFT_MisesMat_sig0); // uniaxial yield stress

    H = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, H, _IFT_MisesMat_h); // hardening modulus
    /*********************************************************************************************************/
    IR_GIVE_FIELD(ir, omega_crit, _IFT_MisesMat_omega_crit); // critical damage

    IR_GIVE_OPTIONAL_FIELD(ir, a, _IFT_MisesMat_a); // exponent in damage law
    /********************************************************************************************************/

    return IRRT_OK;
}

// creates a new material status  corresponding to this class
MaterialStatus *
MisesMat :: CreateStatus(GaussPoint *gp) const
{
    MisesMatStatus *status;
    status = new MisesMatStatus(1, this->giveDomain(), gp);
    return status;
}


// returns the stress vector in 3d stress space
void
MisesMat :: giveRealStressVector(FloatArray &answer,
                                 GaussPoint *gp,
                                 const FloatArray &totalStrain,
                                 TimeStep *tStep)
{
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    //this->initGpForNewStep(gp);
    this->performPlasticityReturn(gp, totalStrain);
    double omega = computeDamage(gp, tStep);
    FloatArray totalStress = status->giveTempEffectiveStress();
    totalStress.times(1 - omega);
    answer = totalStress;
    status->setTempDamage(omega);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}


void
MisesMat :: giveFirstPKStressVector_3d(FloatArray &answer,
                                       GaussPoint *gp,
                                       const FloatArray &totalDefGradOOFEM,
                                       TimeStep *tStep)
{
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);
    //this->initGpForNewStep(gp);

    double kappa, dKappa, yieldValue, mi;
    FloatMatrix F, oldF, invOldF;
    FloatArray totalStrain;
    FloatArray s(6);
    F.beMatrixForm(totalDefGradOOFEM); //(method assumes full 3D)

    kappa = status->giveCumulativePlasticStrain();
    oldF.beMatrixForm( status->giveFVector() );
    invOldF.beInverseOf(oldF);
    //relative deformation radient
    FloatMatrix f;
    f.beProductOf(F, invOldF);
    //compute elastic predictor
    FloatMatrix trialLeftCauchyGreen, help;

    f.times( 1./cbrt(f.giveDeterminant()) );

    help.beProductOf(f, status->giveTempLeftCauchyGreen());
    trialLeftCauchyGreen.beProductTOf(help, f);
    FloatMatrix E;
    E.beTProductOf(F, F);
    E.at(1, 1) -= 1.0;
    E.at(2, 2) -= 1.0;
    E.at(3, 3) -= 1.0;
    E.times(0.5);

    FloatArray e;
    e.beSymVectorFormOfStrain(E);

    StrainVector leftCauchyGreen(_3dMat);
    StrainVector leftCauchyGreenDev(_3dMat);
    double leftCauchyGreenVol;

    leftCauchyGreen.beSymVectorFormOfStrain(trialLeftCauchyGreen);

    leftCauchyGreen.computeDeviatoricVolumetricSplit(leftCauchyGreenDev, leftCauchyGreenVol);
    StressVector trialStressDev(_3dMat);
    leftCauchyGreenDev.applyDeviatoricElasticStiffness(trialStressDev, G / 2.);
    s = trialStressDev;

    //check for plastic loading
    double trialS = trialStressDev.computeStressNorm();
    double sigmaY = sig0 + H * kappa;
    //yieldValue = sqrt(3./2.)*trialS-sigmaY;
    yieldValue = trialS - sqrt(2. / 3.) * sigmaY;


    //store deviatoric trial stress(reused by algorithmic stiffness)
    status->letTrialStressDevBe(trialStressDev);
    //the return-mapping algorithm
    double J = F.giveDeterminant();
    mi = leftCauchyGreenVol * G;
    if ( yieldValue > 0 ) {
        //dKappa =sqrt(3./2.)* yieldValue/(H + 3.*mi);
        //kappa = kappa + dKappa;
        //trialStressDev.times(1-sqrt(6.)*mi*dKappa/trialS);
        dKappa = ( yieldValue / ( 2 * mi ) ) / ( 1 + H / ( 3 * mi ) );
        FloatArray n = trialStressDev;
        n.times(2 * mi * dKappa / trialS);
        ////return map
        s.beDifferenceOf(trialStressDev, n);
        kappa += sqrt(2. / 3.) * dKappa;


        //update of intermediate configuration
        trialLeftCauchyGreen.beMatrixForm(s);
        trialLeftCauchyGreen.times(1.0 / G);
        trialLeftCauchyGreen.at(1, 1) += leftCauchyGreenVol;
        trialLeftCauchyGreen.at(2, 2) += leftCauchyGreenVol;
        trialLeftCauchyGreen.at(2, 2) += leftCauchyGreenVol;
        trialLeftCauchyGreen.times(J * J);
    }

    //addition of the elastic mean stress
    FloatMatrix kirchhoffStress;
    kirchhoffStress.beMatrixForm(s);
    kirchhoffStress.at(1, 1) += 1. / 2. *  K * ( J * J - 1 );
    kirchhoffStress.at(2, 2) += 1. / 2. *  K * ( J * J - 1 );
    kirchhoffStress.at(3, 3) += 1. / 2. *  K * ( J * J - 1 );


    FloatMatrix iF, Ep(3, 3), S;
    FloatArray vF, vS, ep;

    //transform Kirchhoff stress into Second Piola - Kirchhoff stress
    iF.beInverseOf(F);
    help.beProductOf(iF, kirchhoffStress);
    S.beProductTOf(help, iF);

    this->computeGLPlasticStrain(F, Ep, trialLeftCauchyGreen, J);

    ep.beSymVectorFormOfStrain(Ep);
    vS.beSymVectorForm(S);
    vF.beVectorForm(F);
    answer.beVectorForm(kirchhoffStress);

    status->setTrialStressVol(mi);
    status->letTempLeftCauchyGreenBe(trialLeftCauchyGreen);
    status->setTempCumulativePlasticStrain(kappa);
    status->letTempStressVectorBe(answer);
    status->letTempStrainVectorBe(e);
    status->letTempPlasticStrainBe(ep);
    status->letTempPVectorBe(answer);
    status->letTempFVectorBe(vF);
}


void
MisesMat :: computeGLPlasticStrain(const FloatMatrix &F, FloatMatrix &Ep, FloatMatrix b, double J)
{
    FloatMatrix I, invB, help;
    I.resize(3, 3);
    I.beUnitMatrix();
    invB.beInverseOf(b);
    help.beTProductOf(F, invB);
    Ep.beProductOf(help, F);
    Ep.times( pow(J, -2. / 3.) );
    Ep.subtract(I);
    Ep.times(1. / 2.);
}


// returns the stress vector in 3d stress space
// computed from the previous plastic strain and current total strain
void
MisesMat :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain)
{
    double kappa, yieldValue, dKappa;
    FloatArray reducedStress;
    FloatArray strain, plStrain;
    MaterialMode mode = gp->giveMaterialMode();
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    StressVector fullStress(mode);
    this->initTempStatus(gp);
    //this->initGpForNewStep(gp);
    // get the initial plastic strain and initial kappa from the status
    plStrain = status->givePlasticStrain();
    kappa = status->giveCumulativePlasticStrain();

    // === radial return algorithm ===
    if ( mode == _1dMat ) {
        LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
        double E = lmat->give('E', gp);
        /*trial stress*/
        fullStress.at(1) = E * ( totalStrain.at(1) - plStrain.at(1) );
        double sigmaY = sig0 + H * kappa;
        double trialS = fullStress.at(1);
        trialS = fabs(trialS);
        /*yield function*/
        yieldValue = trialS - sigmaY;
        // === radial return algorithm ===
        if ( yieldValue > 0 ) {
            dKappa = yieldValue / ( H + E );
            kappa += dKappa;
            fullStress.at(1) = fullStress.at(1) - dKappa *E *signum( fullStress.at(1) );
            plStrain.at(1) = plStrain.at(1) + dKappa *signum( fullStress.at(1) );
        }
    } else if ( ( mode == _PlaneStrain ) || ( mode == _3dMat ) ) {
        // elastic predictor
        StrainVector elStrain(totalStrain, mode);
        elStrain.subtract(plStrain);
        StrainVector elStrainDev(mode);
        double elStrainVol;
        elStrain.computeDeviatoricVolumetricSplit(elStrainDev, elStrainVol);
        StressVector trialStressDev(mode);
        elStrainDev.applyDeviatoricElasticStiffness(trialStressDev, G);
        /**************************************************************/
        double trialStressVol;
        trialStressVol = 3 * K * elStrainVol;
        /**************************************************************/

        // store the deviatoric and trial stress (reused by algorithmic stiffness)
        status->letTrialStressDevBe(trialStressDev);
        status->setTrialStressVol(trialStressVol);
        // check the yield condition at the trial state
        double trialS = trialStressDev.computeStressNorm();
        double sigmaY = sig0 + H * kappa;
        yieldValue = sqrt(3. / 2.) * trialS - sigmaY;
        if ( yieldValue > 0. ) {
            // increment of cumulative plastic strain
            dKappa = yieldValue / ( H + 3. * G );
            kappa += dKappa;
            StrainVector dPlStrain(mode);
            // the following line is equivalent to multiplication by scaling matrix P
            trialStressDev.applyDeviatoricElasticCompliance(dPlStrain, 0.5);
            // increment of plastic strain
            dPlStrain.times(sqrt(3. / 2.) * dKappa / trialS);
            plStrain.add(dPlStrain);
            // scaling of deviatoric trial stress
            trialStressDev.times(1. - sqrt(6.) * G * dKappa / trialS);
        }

        // assemble the stress from the elastically computed volumetric part
        // and scaled deviatoric part

        double stressVol = 3. * K * elStrainVol;
        trialStressDev.computeDeviatoricVolumetricSum(fullStress, stressVol);
    }

    // store the effective stress in status
    status->letTempEffectiveStressBe(fullStress);
    // store the plastic strain and cumulative plastic strain
    status->letTempPlasticStrainBe(plStrain);
    status->setTempCumulativePlasticStrain(kappa);
}


double
MisesMat :: computeDamageParam(double tempKappa)
{
    double tempDam;
    if ( tempKappa > 0. ) {
        tempDam = omega_crit * ( 1.0 - exp(-a * tempKappa) );
    } else {
        tempDam = 0.;
    }

    return tempDam;
}

double
MisesMat :: computeDamageParamPrime(double tempKappa)
{
    double tempDam;
    if ( tempKappa >= 0. ) {
        tempDam = omega_crit * a * exp(-a * tempKappa);
    } else {
        tempDam = 0.;
    }

    return tempDam;
}



double
MisesMat :: computeDamage(GaussPoint *gp,  TimeStep *tStep)
{
    double tempKappa, dam;
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    dam = status->giveDamage();
    computeCumPlastStrain(tempKappa, gp, tStep);
    double tempDam = computeDamageParam(tempKappa);
    if ( dam > tempDam ) {
        tempDam = dam;
    }

    return tempDam;
}


void MisesMat :: computeCumPlastStrain(double &tempKappa, GaussPoint *gp, TimeStep *tStep)
{
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    tempKappa = status->giveTempCumulativePlasticStrain();
}

void
MisesMat :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    give3dSSMaterialStiffnessMatrix(answer, mode, gp, tStep);
}

void
MisesMat :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    ///@todo Directly compute dPdF instead.
    FloatMatrix dSdE;
    this->give3dLSMaterialStiffnessMatrix(dSdE, mode, gp, tStep);
    this->give_dPdF_from(dSdE, answer, gp);
}

// returns the consistent (algorithmic) tangent stiffness matrix
void
MisesMat :: give3dSSMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
{
    // start from the elastic stiffness
    this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
    if ( mode != TangentStiffness ) {
        return;
    }

    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    double kappa = status->giveCumulativePlasticStrain();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double dKappa = tempKappa - kappa;

    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return;
    }

    // === plastic loading ===

    // yield stress at the beginning of the step
    double sigmaY = sig0 + H * kappa;

    // trial deviatoric stress and its norm
    StressVector trialStressDev(status->giveTrialStressDev(), _3dMat);
    //double trialStressVol = status->giveTrialStressVol();
    double trialS = trialStressDev.computeStressNorm();

    // one correction term
    FloatMatrix stiffnessCorrection(6, 6);
    stiffnessCorrection.beDyadicProductOf(trialStressDev, trialStressDev);
    double factor = -2. * sqrt(6.) * G * G / trialS;
    double factor1 = factor * sigmaY / ( ( H + 3. * G ) * trialS * trialS );
    stiffnessCorrection.times(factor1);
    answer.add(stiffnessCorrection);

    // another correction term
    stiffnessCorrection.bePinvID();
    double factor2 = factor * dKappa;
    stiffnessCorrection.times(factor2);
    answer.add(stiffnessCorrection);

    //influence of damage
    //    double omega = computeDamageParam(tempKappa);
    double omega = status->giveTempDamage();
    answer.times(1. - omega);
    const FloatArray &effStress = status->giveTempEffectiveStress();
    double omegaPrime = computeDamageParamPrime(tempKappa);
    double scalar = -omegaPrime *sqrt(6.) * G / ( 3. * G + H ) / trialS;
    stiffnessCorrection.beDyadicProductOf(effStress, trialStressDev);
    stiffnessCorrection.times(scalar);
    answer.add(stiffnessCorrection);
}

void
MisesMat :: give1dStressStiffMtrx(FloatMatrix &answer,
                                  MatResponseMode mode,
                                  GaussPoint *gp,
                                  TimeStep *tStep)
{
    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    double kappa = status->giveCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double omega = status->giveTempDamage();
    double E = answer.at(1, 1);
    if ( mode != TangentStiffness ) {
        return;
    }


    if ( ( tempKappa - kappa ) <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        answer.times(1 - omega);
        return;
    }


    // === plastic loading ===
    const FloatArray &stressVector = status->giveTempEffectiveStress();
    double stress = stressVector.at(1);
    answer.resize(1, 1);
    answer.at(1, 1) = ( 1 - omega ) * E * H / ( E + H ) - computeDamageParamPrime(tempKappa) * E / ( E + H ) * stress * signum(stress);
}

void
MisesMat :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                     MatResponseMode mode,
                                     GaussPoint *gp,
                                     TimeStep *tStep)
{
    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    if ( mode != TangentStiffness ) {
        return;
    }

    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double kappa = status->giveCumulativePlasticStrain();
    double dKappa = tempKappa - kappa;

    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return;
    }

    // === plastic loading ===
    // yield stress at the beginning of the step
    double sigmaY = sig0 + H * kappa;

    // trial deviatoric stress and its norm
    StressVector trialStressDev(status->giveTrialStressDev(), _PlaneStrain);
    //double trialStressVol = status->giveTrialStressVol();
    double trialS = trialStressDev.computeStressNorm();

    // one correction term
    FloatMatrix stiffnessCorrection(4, 4);
    stiffnessCorrection.beDyadicProductOf(trialStressDev, trialStressDev);
    double factor = -2. * sqrt(6.) * G * G / trialS;
    double factor1 = factor * sigmaY / ( ( H + 3. * G ) * trialS * trialS );
    stiffnessCorrection.times(factor1);
    answer.add(stiffnessCorrection);

    // another correction term
    stiffnessCorrection.resize(4, 4);
    stiffnessCorrection.zero();
    stiffnessCorrection.at(1, 1) = stiffnessCorrection.at(2, 2) = stiffnessCorrection.at(3, 3) = 2. / 3.;
    stiffnessCorrection.at(1, 2) = stiffnessCorrection.at(1, 3) = stiffnessCorrection.at(2, 1) = -1. / 3.;
    stiffnessCorrection.at(2, 3) = stiffnessCorrection.at(3, 1) = stiffnessCorrection.at(3, 2) = -1. / 3.;
    stiffnessCorrection.at(4, 4) = 0.5;
    double factor2 = factor * dKappa;
    stiffnessCorrection.times(factor2);
    answer.add(stiffnessCorrection);

    //influence of damage
    double omega = status->giveTempDamage();
    answer.times(1. - omega);
    const FloatArray &effStress = status->giveTempEffectiveStress();
    double omegaPrime = computeDamageParamPrime(tempKappa);
    double scalar = -omegaPrime *sqrt(6.) * G / ( 3. * G + H ) / trialS;
    stiffnessCorrection.beDyadicProductOf(effStress, trialStressDev);
    stiffnessCorrection.times(scalar);
    answer.add(stiffnessCorrection);
}

void
MisesMat :: give3dLSMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    // start from the elastic stiffness

    FloatMatrix I(6, 6);
    I.at(1, 1) = I.at(2, 2) = I.at(3, 3) = 1;
    I.at(4, 4) = I.at(5, 5) = I.at(6, 6) = 0.5;
    FloatArray delta(6);
    delta.at(1) = delta.at(2) = delta.at(3) = 1;

    FloatMatrix F, F_Tr;
    F.beMatrixForm( status->giveTempFVector() );
    double J;
    J = F.giveDeterminant();

    StressVector trialStressDev(status->giveTrialStressDev(), _3dMat);
    double trialStressVol = status->giveTrialStressVol();
    double trialS = trialStressDev.computeStressNorm();
    FloatArray n = trialStressDev;
    if ( trialS == 0 ) {
        n.resize(6);
    } else {
        n.times(1 / trialS);
    }


    FloatMatrix Cdev(6, 6);
    FloatMatrix C(6, 6);
    FloatMatrix help(6, 6);
    help.beDyadicProductOf(delta, delta);
    C = help;
    help.times(-1. / 3.);
    FloatMatrix C1 = I;
    C1.add(help);
    C1.times(2 * trialStressVol);

    FloatMatrix n1(6, 6), n2(6, 6);
    n1.beDyadicProductOf(n, delta);
    n2.beDyadicProductOf(delta, n);
    help = n1;
    help.add(n2);
    help.times(-2. / 3. * trialS);
    C1.add(help);
    Cdev = C1;
    C.times(K * J * J);

    help = I;
    help.times( -K * ( J * J - 1 ) );
    C.add(help);
    FloatMatrix Cvol = C;
    C.add(C1);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    FloatMatrix invF(3, 3);
    FloatMatrix T(6, 6), tT(6, 6);

    invF.beInverseOf(F);
    //////////////////////////////////////////////////
    //first row of pull back transformation matrix
    T.at(1, 1) = invF.at(1, 1) * invF.at(1, 1);
    T.at(1, 2) = invF.at(1, 2) * invF.at(1, 2);
    T.at(1, 3) = invF.at(1, 3) * invF.at(1, 3);
    T.at(1, 4) = 2. * invF.at(1, 2) * invF.at(1, 3);
    T.at(1, 5) = 2. * invF.at(1, 1) * invF.at(1, 3);
    T.at(1, 6) = 2. * invF.at(1, 1) * invF.at(1, 2);
    //second row of pull back transformation matrix
    T.at(2, 1) = invF.at(2, 1) * invF.at(2, 1);
    T.at(2, 2) = invF.at(2, 2) * invF.at(2, 2);
    T.at(2, 3) = invF.at(2, 3) * invF.at(2, 3);
    T.at(2, 4) = 2. * invF.at(2, 2) * invF.at(2, 3);
    T.at(2, 5) = 2. * invF.at(2, 1) * invF.at(2, 3);
    T.at(2, 6) = 2. * invF.at(2, 1) * invF.at(2, 2);
    //third row of pull back transformation matrix
    T.at(3, 1) = invF.at(3, 1) * invF.at(3, 1);
    T.at(3, 2) = invF.at(3, 2) * invF.at(3, 2);
    T.at(3, 3) = invF.at(3, 3) * invF.at(3, 3);
    T.at(3, 4) = 2. * invF.at(3, 2) * invF.at(3, 3);
    T.at(3, 5) = 2. * invF.at(3, 1) * invF.at(3, 3);
    T.at(3, 6) = 2. * invF.at(3, 1) * invF.at(3, 2);
    //fourth row of pull back transformation matrix
    T.at(4, 1) = invF.at(2, 1) * invF.at(3, 1);
    T.at(4, 2) = invF.at(2, 2) * invF.at(3, 2);
    T.at(4, 3) = invF.at(2, 3) * invF.at(3, 3);
    T.at(4, 4) = ( invF.at(2, 2) * invF.at(3, 3) + invF.at(2, 3) * invF.at(3, 2) );
    T.at(4, 5) = ( invF.at(2, 1) * invF.at(3, 3) + invF.at(2, 3) * invF.at(3, 1) );
    T.at(4, 6) = ( invF.at(2, 1) * invF.at(3, 2) + invF.at(2, 2) * invF.at(3, 1) );
    //fifth row of pull back transformation matrix
    T.at(5, 1) = invF.at(1, 1) * invF.at(3, 1);
    T.at(5, 2) = invF.at(1, 2) * invF.at(3, 2);
    T.at(5, 3) = invF.at(1, 3) * invF.at(3, 3);
    T.at(5, 4) = ( invF.at(1, 2) * invF.at(3, 3) + invF.at(1, 3) * invF.at(3, 2) );
    T.at(5, 5) = ( invF.at(1, 1) * invF.at(3, 3) + invF.at(1, 3) * invF.at(3, 1) );
    T.at(5, 6) = ( invF.at(1, 1) * invF.at(3, 2) + invF.at(1, 2) * invF.at(3, 1) );
    //sixth row of pull back transformation matrix
    T.at(6, 1) = invF.at(1, 1) * invF.at(2, 1);
    T.at(6, 2) = invF.at(1, 2) * invF.at(2, 2);
    T.at(6, 3) = invF.at(1, 3) * invF.at(2, 3);
    T.at(6, 4) = ( invF.at(1, 2) * invF.at(2, 3) + invF.at(1, 3) * invF.at(2, 2) );
    T.at(6, 5) = ( invF.at(1, 1) * invF.at(2, 3) + invF.at(1, 3) * invF.at(2, 1) );
    T.at(6, 6) = ( invF.at(1, 1) * invF.at(2, 2) + invF.at(1, 2) * invF.at(2, 1) );
    ///////////////////////////////////////////

    if ( mode != TangentStiffness ) {
        help.beProductTOf(C, T);
        answer.beProductOf(T, help);
        return;
    }


    double kappa = status->giveCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double dKappa = sqrt(3. / 2.) * ( status->giveTempCumulativePlasticStrain() - kappa );
    //double dKappa = ( status->giveTempCumulativePlasticStrain() - kappa);
    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        help.beProductTOf(C, T);
        answer.beProductOf(T, help);
        return;
    }

    // === plastic loading ===
    //dKappa = dKappa*sqrt(3./2.);
    // trial deviatoric stress and its norm


    double beta0, beta1, beta2, beta3, beta4;
    if ( trialS == 0 ) {
        beta1 = 0;
    } else {
        beta1 = 2 * trialStressVol * dKappa / trialS;
    }

    if ( trialStressVol == 0 ) {
        beta0 = 0;
        beta2 = 0;
        beta3 = beta1;
        beta4 = 0;
    } else {
        beta0 = 1 + H / 3 / trialStressVol;
        beta2 = ( 1 - 1 / beta0 ) * 2. / 3. * trialS * dKappa / trialStressVol;
        beta3 = 1 / beta0 - beta1 + beta2;
        beta4 = ( 1 / beta0 - beta1 ) * trialS / trialStressVol;
    }

    FloatMatrix N;
    N.beDyadicProductOf(n, n);
    N.times(-2 * trialStressVol * beta3);
    answer.resize(6, 6);

    C1.times(-beta1);
    FloatMatrix mN(3, 3);
    mN.at(1, 1) = n.at(1);
    mN.at(1, 2) = n.at(6);
    mN.at(1, 3) = n.at(5);
    mN.at(2, 1) = n.at(6);
    mN.at(2, 2) = n.at(2);
    mN.at(2, 3) = n.at(4);
    mN.at(3, 1) = n.at(5);
    mN.at(3, 2) = n.at(4);
    mN.at(3, 3) = n.at(3);
    FloatMatrix mN2(3, 3);
    mN2.beProductOf(mN, mN);

    double volN2 = 1. / 3. * ( mN2.at(1, 1) + mN2.at(2, 2) + mN2.at(3, 3) );
    FloatArray devN2(6);
    devN2.at(1) = mN2.at(1, 1) - volN2;
    devN2.at(2) = mN2.at(2, 2) - volN2;

    devN2.at(3) = mN2.at(3, 3) - volN2;
    devN2.at(4) = mN2.at(2, 3);
    devN2.at(5) = mN2.at(1, 3);
    devN2.at(6) = mN2.at(1, 2);
    FloatMatrix nonSymPart;
    nonSymPart.beDyadicProductOf(n, devN2);
    //symP.beTranspositionOf(nonSymPart);
    //symP.add(nonSymPart);
    //symP.times(1./2.);
    nonSymPart.times(-2 * trialStressVol * beta4);

    C.add(C1);
    C.add(N);
    C.add(nonSymPart);
    help.beProductTOf(C, T);
    answer.beProductOf(T, help);
}


#ifdef __OOFEG
#endif

int
MisesMat :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        const FloatArray &ep = status->givePlasDef();
        ///@todo Fix this so that it doesn't just fill in zeros for plane stress:
        StructuralMaterial :: giveFullSymVectorForm( answer, ep, gp->giveMaterialMode() );
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveCumulativePlasticStrain();
        return 1;
    } else if ( ( type == IST_DamageScalar ) || ( type == IST_DamageTensor ) ) {
        answer.resize(1);
        answer.at(1) = status->giveDamage();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

//=============================================================================

MisesMatStatus :: MisesMatStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), plasticStrain(), tempPlasticStrain(), trialStressD()
{
    damage = tempDamage = 0.;
    kappa = tempKappa = 0.;
    effStress.resize(6);
    tempEffStress.resize(6);
    leftCauchyGreen.resize(3, 3);
    leftCauchyGreen.at(1, 1) = leftCauchyGreen.at(2, 2) = leftCauchyGreen.at(3, 3) = 1;
}

MisesMatStatus :: ~MisesMatStatus()
{ }

void
MisesMatStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    //int i, n;

    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "status { ");
    /*
     * // this would not be correct, since the status is already updated and kappa is now the "final" value
     * if ( tempKappa > kappa ) {
     * fprintf(file, " Yielding, ");
     * } else {
     * fprintf(file, " Unloading, ");
     * }
     */

    /*********************************************************************************/

    fprintf(file, "damage %.4e", tempDamage);

    /********************************************************************************/
    /*
     * // print the plastic strain
     *  n = plasticStrain.giveSize();
     *  fprintf(file, " plastic_strains ");
     *  for ( i = 1; i <= n; i++ ) {
     *      fprintf( file, " % .4e", plasticStrain.at(i) );
     *  }
     */
    // print the cumulative plastic strain
    fprintf(file, ", kappa ");
    fprintf(file, " % .4e", kappa);

    fprintf(file, "}\n");
    /*
     * //print Left Cauchy - Green deformation tensor
     * fprintf(file," left Cauchy Green");
     * fprintf( file, " % .4e",tempLeftCauchyGreen.at(1,1) );
     * fprintf( file, " % .4e",tempLeftCauchyGreen.at(2,2) );
     * fprintf( file, " % .4e",tempLeftCauchyGreen.at(3,3) );
     * fprintf( file, " % .4e",tempLeftCauchyGreen.at(2,3) );
     * fprintf( file, " % .4e",tempLeftCauchyGreen.at(1,3) );
     * fprintf( file, " % .4e",tempLeftCauchyGreen.at(1,2) );
     *
     * //print deformation gradient
     * fprintf(file," Deformation Gradient");
     * fprintf( file, " % .4e",tempDefGrad.at(1,1) );
     * fprintf( file, " % .4e",tempDefGrad.at(1,2) );
     * fprintf( file, " % .4e",tempDefGrad.at(1,3) );
     * fprintf( file, " % .4e",tempDefGrad.at(2,1) );
     * fprintf( file, " % .4e",tempDefGrad.at(2,2) );
     * fprintf( file, " % .4e",tempDefGrad.at(2,3) );
     * fprintf( file, " % .4e",tempDefGrad.at(3,1) );
     * fprintf( file, " % .4e",tempDefGrad.at(3,2) );
     * fprintf( file, " % .4e",tempDefGrad.at(3,3) );
     */

    fprintf(file, "}\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void MisesMatStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrain.giveSize() == 0 ) {
        plasticStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        plasticStrain.zero();
    }

    tempDamage = damage;
    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
    trialStressD.clear(); // to indicate that it is not defined yet
    tempLeftCauchyGreen = leftCauchyGreen;
}


// updates internal variables when equilibrium is reached
void
MisesMatStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    kappa = tempKappa;
    damage = tempDamage;
    trialStressD.clear(); // to indicate that it is not defined any more
    leftCauchyGreen = tempLeftCauchyGreen;
}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
MisesMatStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data

    // write plastic strain (vector)
    if ( ( iores = plasticStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write cumulative plastic strain (scalar)
    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write damage (scalar)
    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}



contextIOResultType
MisesMatStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read plastic strain (vector)
    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read cumulative plastic strain (scalar)
    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read damage (scalar)
    if ( !stream.read(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem
