/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2018 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "noAdvection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentAdvectionSchemes
{
    defineTypeNameAndDebug(noAdvection, 0);

    addToRunTimeSelectionTable
    (
        momentAdvection,
        noAdvection,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentAdvectionSchemes::noAdvection::noAdvection
(
    const dictionary& dict,
    const quadratureApproximation& quadrature,
    const surfaceScalarField& phi,
    const word& support
)
:
    momentAdvection(dict, quadrature, phi, support)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentAdvectionSchemes::noAdvection::~noAdvection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::momentAdvectionSchemes::noAdvection::realizableCo() const
{
    return scalar(1);
}

Foam::scalar Foam::momentAdvectionSchemes::noAdvection::CoNum() const
{
    return scalar(0);
}

void Foam::momentAdvectionSchemes::noAdvection::update
(
    const bool localPhi,
    const bool wallCollisions
)
{
    return;
}

void Foam::momentAdvectionSchemes::noAdvection::update
(
    const mappedPtrList<volVectorField>& Us,
    const bool wallCollisions
)
{
    return;
}

// ************************************************************************* //
