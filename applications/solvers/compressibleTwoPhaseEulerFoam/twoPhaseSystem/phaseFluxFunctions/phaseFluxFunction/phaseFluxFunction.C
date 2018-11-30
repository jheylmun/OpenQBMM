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

#include "phaseFluxFunction.H"
#include "surfaceInterpolate.H"
#include "fvMatrix.H"
#include "subCycle.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseFluxFunction, 0);
    defineRunTimeSelectionTable(phaseFluxFunction, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxFunction::phaseFluxFunction
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    mesh_(mesh),
    dict_
    (
        mesh.schemesDict().subDict
        (
            IOobject::groupName("compressible", phaseName)
        )
    ),
    own_
    (
        IOobject
        (
            "own",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("own", dimless, 1.0)
    ),
    nei_
    (
        IOobject
        (
            "nei",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("nei", dimless, -1.0)
    ),
    residualAlpha_(dict_.lookupOrDefault("residualAlpha", 0)),
    epsilon_(dict_.lookupOrDefault("cut)ffMa", 1e-10)),
    alphaf_
    (
        IOobject
        (
            IOobject::groupName("alphaf", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    Uf_
    (
        IOobject
        (
            IOobject::groupName("Uf", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, Zero)
    ),
    phi_
    (
        IOobject
        (
            IOobject::groupName("phi", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimVelocity*dimArea, 0.0)
    ),
    pf_
    (
        IOobject
        (
            IOobject::groupName("pf", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxFunction::~phaseFluxFunction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::surfaceScalarField& Foam::phaseFluxFunction::phi() const
{
    return phi_;
}

Foam::tmp<Foam::surfaceScalarField> Foam::phaseFluxFunction::alphaPhi() const
{
    return phi_*alphaf_;
}

Foam::tmp<Foam::volVectorField> Foam::phaseFluxFunction::gradAlpha() const
{
    return fvc::surfaceIntegrate(mesh_.Sf()*alphaf_);
}

Foam::tmp<Foam::volTensorField> Foam::phaseFluxFunction::gradU() const
{
    return fvc::surfaceIntegrate(mesh_.Sf()*Uf_);
}

Foam::tmp<Foam::volVectorField> Foam::phaseFluxFunction::gradp() const
{
    return fvc::surfaceIntegrate(mesh_.Sf()*pf_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
