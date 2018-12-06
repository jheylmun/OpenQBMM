/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Jeff Heylmun
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "phaseFluxIntegrator.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseFluxIntegrator, 0);
    defineRunTimeSelectionTable(phaseFluxIntegrator, dictionary);
}

void Foam::phaseFluxIntegrator::setCoeffs
(
    boolList& storeFields,
    boolList& storeDeltas
)
{
    List<scalarList> c = coeffs();
    List<scalarList> f = Fcoeffs();

    for (label i = 0; i < coeffs().size(); i++)
    {
        for (label j = 0; j < c[i].size() - 1; j++)
        {
            if (mag(c[i][j]) > SMALL)
            {
                storeFields[j] = true;
            }
            if (mag(f[i][j]) > SMALL)
            {
                storeDeltas[j] = true;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxIntegrator::phaseFluxIntegrator
(
    phaseModel& phase1,
    phaseModel& phase2
)
:
    phase1_
    (
        (phase2.granular() || phase2.slavePressure())
      ? phase2 : phase1
    ),
    phase2_
    (
        (phase2.granular() || phase2.slavePressure())
      ? phase1 : phase2
    ),
    gradAlpha_(phase1_.gradAlpha())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxIntegrator::~phaseFluxIntegrator()
{}


// * * * * * * * * * * * * * Public Member Fucntions * * * * * * * * * * * * //


void Foam::phaseFluxIntegrator::integrateFluxes
(
    const dimensionedVector& g,
    volVectorField& Ui,
    volScalarField& pi
)
{
    List<scalarList> c = coeffs();
    List<scalarList> f = Fcoeffs();

    const dimensionedScalar& deltaT = Ui.mesh().time().deltaT();

    const volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;

    label sign = 1;
    if (&phase1_ != &(phase1_.fluid().phase1()))
    {
        sign = -1;
    }

    for (label stepi = 0; stepi < nSteps(); stepi++)
    {
        volVectorField Fh(phase1_.fluid().F()*sign);

        phase1_.updateFluxes();
        phase2_.updateFluxes(1.0 - phase1_.alphaf());

        phase1_.advect
        (
            stepi,
            c[stepi],
            f[stepi],
            deltaT,
            g,
            Fh,
            Ui,
            pi
        );
        phase1_.decode();
        //phase1_.encode();

        phase2_.advect
        (
            stepi,
            c[stepi],
            f[stepi],
            deltaT,
            g,
            -Fh,
            Ui,
            pi
        );
        alpha2 = 1.0 - alpha1;
        alpha2.correctBoundaryConditions();
        phase2_.decode();
        //phase2_.encode();


        Ui = phase1_.fluid().mixtureU();
        pi = phase1_.fluid().mixturep();
    }
}

// ************************************************************************* //
