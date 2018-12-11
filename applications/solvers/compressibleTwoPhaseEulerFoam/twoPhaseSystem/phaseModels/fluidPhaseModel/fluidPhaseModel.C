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

#include "fluidPhaseModel.H"
#include "twoPhaseSystem.H"
#include "diameterModel.H"
#include "fvMatrix.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        fluidPhaseModel,
        dictionary,
        fluid
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidPhaseModel::fluidPhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    phaseModel(fluid, phaseProperties, phaseName),
    p_
    (
        IOobject
        (
            IOobject::groupName("p", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->thermoPtr_->p(),
        this->thermoPtr_->p().boundaryField()
    ),
    alphaRho_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        (*this)*rho(),
        this->boundaryField()
    ),
    alphaRhoU_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoU", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        alphaRho_*U_,
        U_.boundaryField()
    ),
    alphaRhoE_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoE", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        alphaRho_*E_,
        E_.boundaryField()
    ),
    massFlux_
    (
        IOobject
        (
            IOobject::groupName("massFlux", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fvc::flux(U_)*fvc::interpolate((*this)*rho())
    ),
    momentumFlux_
    (
        IOobject
        (
            IOobject::groupName("momentumFlux", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        massFlux_*fvc::interpolate(U_)
    ),
    energyFlux_
    (
        IOobject
        (
            IOobject::groupName("energyFlux", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        massFlux_*fvc::interpolate(E_)
    ),
    fluxFunction_
    (
        fluidFluxFunction::New(fluid.mesh(), phaseName)
    )
{
    setTurbulenceModel();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidPhaseModel::~fluidPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluidPhaseModel::setNSteps
(
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    label nSteps = storeFields.size();
    this->storedFieldIndexes_ = labelList(nSteps, -1);
    this->storedDeltaIndexes_ = labelList(nSteps, -1);
    bool setAlpha = !otherPhase().granular();

    label currFieldIndex = 0;
    label currDeltaIndex = 0;
    for (label stepi = 0; stepi < nSteps; stepi++)
    {
        if (storeFields[stepi])
        {
            this->storedFieldIndexes_[stepi] = currFieldIndex;
            currFieldIndex++;

            if (setAlpha)
            {
                alphas_.append
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "alpha"+Foam::name(stepi),
                                name()
                            ),
                            fluid_.mesh().time().timeName(),
                            fluid_.mesh()
                        ),
                        fluid_.mesh(),
                        dimensionedScalar("zero", dimless, 0.0)
                    )
                );
            }
            alphaRhos_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "alphaRho" + Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedScalar("zero", alphaRho_.dimensions(), 0.0)
                )
            );
            alphaRhoUs_.append
            (
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "alphaRhoU" + Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedVector("zero", alphaRhoU_.dimensions(), Zero)
                )
            );
            alphaRhoEs_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "alphaRhoE" + Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedScalar("zero", alphaRhoE_.dimensions(), 0.0)
                )
            );

        }

        if (storeDeltas[stepi])
        {
            this->storedDeltaIndexes_[stepi] = currDeltaIndex;
            currDeltaIndex++;

            if (setAlpha)
            {

                deltaAlphas_.append
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "deltaAlpha"+Foam::name(stepi),
                                name()
                            ),
                            fluid_.mesh().time().timeName(),
                            fluid_.mesh()
                        ),
                        fluid_.mesh(),
                        dimensionedScalar("zero", dimless/dimTime, 0.0)
                    )
                );
            }
            deltaAlphaRhos_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "deltaAlphaRho" + Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedScalar("zero", alphaRho_.dimensions()/dimTime, 0.0)
                )
            );
            deltaAlphaRhoUs_.append
            (
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "deltaAlphaRhoU" + Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedVector("zero", alphaRhoU_.dimensions()/dimTime, Zero)
                )
            );
            deltaAlphaRhoEs_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "deltaAlphaRhoE" + Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedScalar("zero", alphaRhoE_.dimensions()/dimTime, 0.0)
                )
            );

        }
    }
}

void Foam::fluidPhaseModel::correctP()
{
    volScalarField R = constant::thermodynamic::RR/thermoPtr_->W();
    R.dimensions().reset
    (
        dimPressure/dimDensity/dimTemperature
    );

    if (thermoPtr_->he().name()[0] == 'e')
    {
        p_ = rho_*R*thermoPtr_->he()/thermoPtr_->Cv();
    }
    else if (thermoPtr_->he().name()[0] == 'h')
    {
        p_ = rho_*R*thermoPtr_->he()/thermoPtr_->Cp();
    }

}

Foam::tmp<Foam::volScalarField> Foam::fluidPhaseModel::c() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("c", name()),
            sqrt(thermoPtr_->gamma()*p_/rho_)
        )
    );
}

void Foam::fluidPhaseModel::advect
(
    const label stepi,
    const scalarList& coeffs,
    const scalarList& Fcoeffs,
    const dimensionedScalar& deltaT,
    const dimensionedVector& g,
    const volVectorField& F,
    const volVectorField& Ui,
    const volScalarField& pi
)
{
    label fi = this->storedFieldIndexes_[stepi];
    label di = this->storedDeltaIndexes_[stepi];

    if (fi != -1)
    {
        if (!otherPhase().granular())
        {
            alphas_[fi] = *this;
        }
        alphaRhos_[fi] = alphaRho_;
        alphaRhoUs_[fi] = alphaRhoU_;
        alphaRhoEs_[fi] = alphaRhoE_;
    }

    tmp<volScalarField> alpha;
    tmp<volScalarField> deltaAlpha;
    if (!otherPhase().granular())
    {
        alpha = (*this)*coeffs[stepi];
        deltaAlpha = -(Ui & gradAlpha_);
    }

    volScalarField deltaAlphaRho(-fvc::div(massFlux_));
    volVectorField deltaAlphaRhoU
    (
       - fvc::div(momentumFlux_)
       + pi*gradAlpha_
       + F
       + alphaRho_*g
    );
    volScalarField deltaAlphaRhoE
    (
        -fvc::div(energyFlux_) + ((*this)*rho_*g & U_)
    );

    if (otherPhase().granular())
    {
        deltaAlphaRhoE +=
            (F & otherPhase().U())
          - pi*fvc::div(otherPhase().alphaPhi());
    }
    else
    {
        deltaAlphaRhoE += (pi*Ui) & gradAlpha_;
    }

    if (di != -1)
    {
        if (!otherPhase().granular())
        {
            deltaAlphas_[di] = deltaAlpha;
        }
        deltaAlphaRhos_[di] = deltaAlphaRho;
        deltaAlphaRhoUs_[di] = deltaAlphaRhoU;
        deltaAlphaRhoEs_[di] = deltaAlphaRhoE;
    }

    if (!otherPhase().granular())
    {
        alpha.ref() *= coeffs[stepi];
        deltaAlpha.ref() *= Fcoeffs[stepi];
    }
    volScalarField alphaRho(alphaRho_*coeffs[stepi]);
    volVectorField alphaRhoU(alphaRhoU_*coeffs[stepi]);
    volScalarField alphaRhoE(alphaRhoE_*coeffs[stepi]);

    deltaAlphaRho *= Fcoeffs[stepi];
    deltaAlphaRhoU *= Fcoeffs[stepi];
    deltaAlphaRhoE *= Fcoeffs[stepi];

    label fieldi = 0;
    label deltai = 0;
    for (label i = 0; i < stepi; i++)
    {
        if (this->storedFieldIndexes_[i] != -1)
        {
            if (!otherPhase().granular())
            {
                alpha.ref() += alphas_[fieldi]*coeffs[i];
            }
            alphaRho += alphaRhos_[fieldi]*coeffs[i];
            alphaRhoU += alphaRhoUs_[fieldi]*coeffs[i];
            alphaRhoE += alphaRhoEs_[fieldi]*coeffs[i];
            fieldi++;
        }

        if (this->storedDeltaIndexes_[i] != -1)
        {
            if (!otherPhase().granular())
            {
                deltaAlpha.ref() += deltaAlphas_[deltai]*Fcoeffs[i];
            }
            deltaAlphaRho += deltaAlphaRhos_[deltai]*Fcoeffs[i];
            deltaAlphaRhoU += deltaAlphaRhoUs_[deltai]*Fcoeffs[i];
            deltaAlphaRhoE += deltaAlphaRhoEs_[deltai]*Fcoeffs[i];
            deltai++;
        }
    }
    if (!otherPhase().granular())
    {
        refCast<volScalarField>(*this) += deltaT*deltaAlpha;
        this->correctBoundaryConditions();
    }

    alphaRho_ = alphaRho + deltaT*deltaAlphaRho;
    alphaRho_.correctBoundaryConditions();

    alphaRhoU_ = alphaRhoU + deltaT*deltaAlphaRhoU;
    alphaRhoU_.correctBoundaryConditions();

    alphaRhoE_ = alphaRhoE + deltaT*deltaAlphaRhoE;
    alphaRhoE_.correctBoundaryConditions();
}


void Foam::fluidPhaseModel::updateFluxes()
{
    volScalarField H(IOobject::groupName("H", name()), E_ + p_/rho_);

    this->fluxFunction_->updateFluxes
    (
        massFlux_,
        momentumFlux_,
        energyFlux_,
        *this,
        rho_,
        U_,
        E_,
        p_,
        c(),
        fluid_.U(),
        fluid_.p()
    );
    gradAlpha_ = fluxFunction_->gradAlpha();
}

void Foam::fluidPhaseModel::updateFluxes(const surfaceScalarField& alphaf)
{
    volScalarField H(IOobject::groupName("H", name()), E_ + p_/rho_);
    this->fluxFunction_->updateFluxes
    (
        massFlux_,
        momentumFlux_,
        energyFlux_,
        alphaf,
        rho_,
        U_,
        E_,
        p_,
        c(),
        fluid_.U(),
        fluid_.p()
    );
    gradAlpha_ = fluxFunction_->gradAlpha();
}


void Foam::fluidPhaseModel::decode()
{
    this->max(residualAlpha_);
    rho_ = alphaRho_/(*this);
    rho_.correctBoundaryConditions();

    U_ = alphaRhoU_/alphaRho_;
    U_.correctBoundaryConditions();

    E_ = alphaRhoE_/alphaRho_;
    E_.correctBoundaryConditions();

    he_ = E_ - 0.5*magSqr(U_);
    he_.correctBoundaryConditions();

    correctP();
    p_.correctBoundaryConditions();
    thermoPtr_->p() = p_;
    thermoPtr_->correct();
}


void Foam::fluidPhaseModel::encode()
{
    (*this).max(residualAlpha_);

    alphaRho_ = (*this)*rho_;
    alphaRho_.correctBoundaryConditions();

    alphaRhoU_ = alphaRho_*U_;
    alphaRhoU_.correctBoundaryConditions();

    E_ = he_ + 0.5*magSqr(U_);
    alphaRhoE_ = alphaRho_*E_;
    alphaRhoE_.correctBoundaryConditions();
}


void Foam::fluidPhaseModel::correctThermo()
{
    E_ = he_ + 0.5*magSqr(U_);
    E_.correctBoundaryConditions();

    alphaRhoE_ = alphaRho_*E_;
    alphaRhoE_.correctBoundaryConditions();

    correctP();
    p_.correctBoundaryConditions();
    thermoPtr_->p() = p_;
    thermoPtr_->correct();
}

// ************************************************************************* //
