/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

    #include "RK3SSPPhase.H"
#include "twoPhaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxIntegrators
{

    defineTypeNameAndDebug(RK3SSPPhase, 0);
    addToRunTimeSelectionTable(phaseFluxIntegrator, RK3SSPPhase, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxIntegrators::RK3SSPPhase::RK3SSPPhase
(
    phaseModel& phase1,
    phaseModel& phase2
)
:
    phaseFluxIntegrator(phase1, phase2)
{
    boolList storeFields(nSteps(), false);
    boolList storeDeltas(nSteps(), false);
    setCoeffs(storeFields, storeDeltas);

    phase1_.setNSteps(storeFields, storeDeltas);
    phase2_.setNSteps(storeFields, storeDeltas);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxIntegrators::RK3SSPPhase::~RK3SSPPhase()
{}


// * * * * * * * * * * *  Protected Member Fucntions * * * * * * * * * * * * //

Foam::List<Foam::scalarList>
Foam::phaseFluxIntegrators::RK3SSPPhase::coeffs() const
{
    return {{1.0}, {3.0/4.0, 1.0/4.0}, {1.0/3.0, 0.0, 2.0/3.0}};
}

Foam::List<Foam::scalarList>
Foam::phaseFluxIntegrators::RK3SSPPhase::Fcoeffs() const
{
    return {{1.0}, {0.0, 1.0/4.0}, {0.0, 0.0, 2.0/3.0}};
}

// ************************************************************************* //
