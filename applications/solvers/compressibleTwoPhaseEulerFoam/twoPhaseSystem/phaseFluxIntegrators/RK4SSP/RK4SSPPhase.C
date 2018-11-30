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

#include "RK4SSPPhase.H"
#include "twoPhaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxIntegrators
{

    defineTypeNameAndDebug(RK4SSPPhase, 0);
    addToRunTimeSelectionTable(phaseFluxIntegrator, RK4SSPPhase, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxIntegrators::RK4SSPPhase::RK4SSPPhase
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

Foam::phaseFluxIntegrators::RK4SSPPhase::~RK4SSPPhase()
{}


// * * * * * * * * * * *  Protected Member Fucntions * * * * * * * * * * * * //

Foam::List<Foam::scalarList>
Foam::phaseFluxIntegrators::RK4SSPPhase::coeffs() const
{
    List<scalarList> c(10);
    c[0] = {1.0};
    for (label i = 1; i < 4; i++)
    {
        c[i] = scalarList(i + 1, 0.0);
        c[i][i] = 1.0;
    }
    c[4] = {3.0/5.0, 0.0, 0.0, 0.0, 2.0/5.0};
    for (label i = 5; i < 9; i++)
    {
        c[i] = scalarList(i + 1, 0.0);
        c[i][i] = 1.0;
    }
    c[9] = {1.0/25.0, 0, 0, 0, 9.0/25.0, 0, 0, 0, 0, 3.0/5.0};
    return c;
}

Foam::List<Foam::scalarList>
Foam::phaseFluxIntegrators::RK4SSPPhase::Fcoeffs() const
{
    List<scalarList> f(10);
    f[0] = {1.0/6.0};
    for (label i = 1; i < 4; i++)
    {
        f[i] = scalarList(i + 1, 0.0);
        f[i][i] = 1.0/6.0;
    }
    f[4] = {0.0, 0.0, 0.0, 0.0, 1.0/15.0};
    for (label i = 5; i < 9; i++)
    {
        f[i] = scalarList(i + 1, 0.0);
        f[i][i] = 1.0/6.0;
    }
    f[9] = {0, 0, 0, 0, 3.0/50.0, 0, 0 , 0, 0, 1.0/10.0};
    return f;
}

// ************************************************************************* //
