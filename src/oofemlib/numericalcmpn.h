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
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//
// FILE: numericalcmpn.h
//

#ifndef numericalcmpn_h
#define numericalcmpn_h

/**
 * Type representing numerical component. The components of characteristic equations are mapped
 * to their corresponding numrical counterparts using these common component types.
 * All numerical methods solving the same problem have to use the same and compulsory
 * NumericalCmpn values. This allows to use generally any numerical method instance (even added in future)
 * without changing any code.
 */
enum NumericalCmpn {
    LinearEquationLhs,
    LinearEquationRhs,   /* for linear equation problem in the form Ax=b */
    LinearEquationSolution,
    EigenValues,
    EigenVectors,
    AEigvMtrx,
    BEigvMtrx,
    NumberOfEigenValues,
    PrescribedTolerancy,
    ReachedLevel,
    IncrementOfSolution,
    RequiredIterations,
    NonLinearLhs,
    NonLinearRhs_Total,
    NonLinearRhs_Incremental,
    InitialNonLinearRhs,
    TotalNonLinearSolution,
    StepLength,
    CurrentLevel,
    InternalRhs,
    IncrementOfNonlinearSolution,
    InternalState
};

#endif // numericalcmpn_h
