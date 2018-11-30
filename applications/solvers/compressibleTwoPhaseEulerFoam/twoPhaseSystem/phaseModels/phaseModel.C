/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "phaseModel.H"
#include "twoPhaseSystem.H"
#include "diameterModel.H"
#include "fvMatrix.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "dragModel.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "fixedValueFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseModel, 0);
    defineRunTimeSelectionTable(phaseModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("alpha", dimless, 0)
    ),
    fluid_(fluid),
    name_(phaseName),
    phaseDict_
    (
        phaseProperties.subDict(name_)
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        fluid.subDict(phaseName).lookup("residualAlpha")
    ),
    residualRho_
    (
        "residualRho",
        dimDensity,
        fluid.subDict(phaseName).lookup("residualRho")
    ),
    thermoPtr_
    (
        rhoThermo::New(fluid.mesh(), phaseName)
    ),
    rho_(thermoPtr_->rho()),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    he_(thermoPtr_->he()),
    E_
    (
        IOobject
        (
            IOobject::groupName("E", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        he_ + 0.5*magSqr(U_),
        thermoPtr_->T().boundaryField().types()
    ),
    gradAlpha_
    (
        IOobject
        (
            IOobject::groupName("gradAlpha", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fvc::grad(*this)
    )
{
    thermoPtr_->validate("phaseModel " + name_, "h", "e");
    rho_.writeOpt() = IOobject::AUTO_WRITE;

    const word phiName = IOobject::groupName("phi", name_);

    IOobject phiHeader
    (
        phiName,
        fluid_.mesh().time().timeName(),
        fluid_.mesh(),
        IOobject::NO_READ
    );

    dPtr_ = diameterModel::New
    (
        phaseDict_,
        *this
    );
}

void Foam::phaseModel::setTurbulenceModel()
{
    turbulence_ =
        PhaseCompressibleTurbulenceModel<phaseModel>::New
        (
            *this,
            this->rho(),
            this->U(),
            this->massFlux(),
            this->phi(),
            *this
        );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseModel& Foam::phaseModel::otherPhase() const
{
    return fluid_.otherPhase(*this);
}


const Foam::PhaseCompressibleTurbulenceModel<Foam::phaseModel>&
Foam::phaseModel::turbulence() const
{
    return turbulence_();
}

Foam::PhaseCompressibleTurbulenceModel<Foam::phaseModel>&
Foam::phaseModel::turbulence()
{
    return turbulence_();
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::d() const
{
    return dPtr_().d();
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::nuEff() const
{
    return turbulence_->nuEff();
}


Foam::tmp<Foam::scalarField> Foam::phaseModel::nuEff(const label patchi) const
{
    return turbulence_->nuEff(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::k() const
{
    return turbulence_->k();
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::epsilon() const
{
    return turbulence_->epsilon();
}


Foam::tmp<Foam::volSymmTensorField> Foam::phaseModel::R() const
{
    return turbulence_->R();
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::pPrime() const
{
    return turbulence_->pPrime();
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseModel::pPrimef() const
{
    return turbulence_->pPrimef();
}


Foam::tmp<Foam::volSymmTensorField> Foam::phaseModel::devRhoReff() const
{
    return turbulence_->devRhoReff();
}


Foam::tmp<Foam::fvVectorMatrix> Foam::phaseModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return turbulence_->divDevRhoReff(U);
}


void Foam::phaseModel::correct()
{
    return dPtr_->correct();
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::dissipationCoeff()const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ProductionCoeff",
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
            dimensionedScalar("zero", dimensionSet(1, -2, -2, 0, 0), 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::productionCoeff()const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "productionCoeff",
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
            dimensionedScalar("zero", dimensionSet(1, -2, -2, 0, 0), 0.0)
        )
    );
}


bool Foam::phaseModel::read(const dictionary& phaseProperties)
{
    phaseDict_ = phaseProperties.subDict(name_);
    return dPtr_->read(phaseDict_);
}


void Foam::phaseModel::correctInflowOutflow(surfaceScalarField& alphaPhi) const
{
    surfaceScalarField::Boundary& alphaPhiBf = alphaPhi.boundaryFieldRef();
    const volScalarField::Boundary& alphaBf = boundaryField();
    const surfaceScalarField::Boundary& phiBf = phi().boundaryField();

    forAll(alphaPhiBf, patchi)
    {
        fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];

        if (!alphaPhip.coupled())
        {
            alphaPhip = phiBf[patchi]*alphaBf[patchi];
        }
    }
}


// ************************************************************************* //
