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

#include "granularPhaseModel.H"
#include "twoPhaseSystem.H"
#include "diameterModel.H"
#include "fvMatrix.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        granularPhaseModel,
        dictionary,
        granular
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularPhaseModel::granularPhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    phaseModel(fluid, phaseProperties, phaseName),
    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            this->phaseDict_
        )
    ),
    conductivityModel_
    (
        kineticTheoryModels::conductivityModel::New
        (
            this->phaseDict_
        )
    ),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            this->phaseDict_
        )
    ),
    granularPressureModel_
    (
        kineticTheoryModels::granularPressureModel::New
        (
            this->phaseDict_
        )
    ),
    frictionalStressModel_
    (
        kineticTheoryModels::frictionalStressModel::New
        (
            this->phaseDict_
        )
    ),

    e_("e", dimless, this->phaseDict_),
    alphaMax_("alphaMax", dimless, this->phaseDict_),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        this->phaseDict_
    ),

    maxNut_
    (
        "maxNut",
        dimensionSet(0,2,-1,0,0),
        this->phaseDict_.lookupOrDefault<scalar>("maxNut",1000)
    ),

    Theta_
    (
        IOobject
        (
            IOobject::groupName("Theta", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),

    Ps_
    (
        IOobject
        (
            IOobject::groupName("Ps", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (
            Theta_
           *granularPressureModel_->granularPressureCoeff
            (
                (*this),
                radialModel_->g0(*this, alphaMinFriction_, alphaMax_),
                this->rho_,
                e_
            )
        )
    ),
    Pfric_
    (
        IOobject
        (
            IOobject::groupName("Pfric", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        frictionalStressModel_->frictionalPressure
        (
            *this,
            alphaMinFriction_,
            alphaMax_
        )
    ),

    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phaseName),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),

    gs0_
    (
        IOobject
        (
            IOobject::groupName("gs0", phaseName),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phaseName),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),

    nut_
    (
        IOobject
        (
            IOobject::groupName("nut", phaseName),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),

    nuFric_
    (
        IOobject
        (
            IOobject::groupName("nuFric", phaseName),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
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
        this->boundaryField().types()
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
        U_.boundaryField().types()
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
        E_.boundaryField().types()
    ),
    alphaRhoPTE_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPTE", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        1.5*alphaRho_*Theta_,
        Theta_.boundaryField().types()
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
    PTEFlux_
    (
        IOobject
        (
            IOobject::groupName("PTEFlux", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        1.5*massFlux_*fvc::interpolate(Theta_)
    ),
    fluxFunction_
    (
        granularFluxFunction::New(fluid.mesh(), phaseName)
    )
{
    setTurbulenceModel();

    // Kinetic energy is not included in E
    E_ = he_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularPhaseModel::~granularPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::granularPhaseModel::setNSteps
(
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    label nSteps = storeFields.size();
    this->storedFieldIndexes_ = labelList(nSteps, -1);
    this->storedDeltaIndexes_ = labelList(nSteps, -1);

    label currFieldIndex = 0;
    label currDeltaIndex = 0;
    for (label stepi = 0; stepi < nSteps; stepi++)
    {
        if (storeFields[stepi])
        {
            this->storedFieldIndexes_[stepi] = currFieldIndex;
            currFieldIndex++;

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
            alphaRhoPTEs_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "alphaRhoPTE" + Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedScalar("zero", alphaRhoPTE_.dimensions(), 0.0)
                )
            );

        }

        if (storeDeltas[stepi])
        {
            this->storedDeltaIndexes_[stepi] = currDeltaIndex;
            currDeltaIndex++;

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

            deltaAlphaRhoPTEs_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "deltaAlphaRhoPTE" + Foam::name(stepi),
                            name()
                        ),
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh()
                    ),
                    fluid_.mesh(),
                    dimensionedScalar("zero", alphaRhoPTE_.dimensions()/dimTime, 0.0)
                )
            );

        }
    }
}

Foam::tmp<Foam::volScalarField> Foam::granularPhaseModel::c() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("a", name()),
            sqrt
            (
                pPrime()/rho_
            + 2.0/3.0*sqr
                (
                    granularPressureModel_->granularPressureCoeff
                    (
                        *this,
                        gs0_,
                        rho_,
                        e_
                    )
                )*Theta_/sqr(Foam::max(alphaRho_, residualAlphaRho()))
            )
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::granularPhaseModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField> Foam::granularPhaseModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField> Foam::granularPhaseModel::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", name()),
                this->fluid_.mesh().time().timeName(),
                this->fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (nut_)*dev(twoSymm(fvc::grad(U())))
          - (lambda_*fvc::div(U()))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::granularPhaseModel::pPrime() const
{
    tmp<volScalarField> tpPrime
    (
        Theta_
       *granularPressureModel_->granularPressureCoeffPrime
        (
            *this,
            radialModel_->g0(*this, alphaMinFriction_, alphaMax_),
            radialModel_->g0prime(*this, alphaMinFriction_, alphaMax_),
            this->rho_,
            e_
        )
     +  frictionalStressModel_->frictionalPressurePrime
        (
            *this,
            alphaMinFriction_,
            alphaMax_
        )
    );

    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField> Foam::granularPhaseModel::pPrimef() const
{
    return fvc::interpolate(pPrime());
}


Foam::tmp<Foam::volSymmTensorField>
Foam::granularPhaseModel::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", this->name()),
                this->fluid_.mesh().time().timeName(),
                this->fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (this->rho_*nut_)*dev(twoSymm(fvc::grad(this->U_)))
          - ((this->rho_*lambda_)*fvc::div(phi()))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix> Foam::granularPhaseModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(rho()*nut_, U)
      - fvc::div
        (
            (this->rho_*nut_)*dev2(::Foam::T(fvc::grad(U)))
          + ((this->rho_*lambda_)*fvc::div(U))
           *dimensioned<symmTensor>("I", dimless, symmTensor::I)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::dissipationCoeff() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            12.0*(1.0 - sqr(e_))*gs0_*sqr(*this)*rho_
           /(sqrt(Foam::constant::mathematical::pi)*d())
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::granularPhaseModel::productionCoeff() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "dissipationCoeff",
            81.0*(*this)*sqr(otherPhase().mu())*magSqr(U_ - otherPhase().U())
           /(
                gs0_*pow3(d())*rho_*sqrt(Foam::constant::mathematical::pi)
              + dimensionedScalar("small", dimMass, vSmall)
            )
        )
    );
}


void Foam::granularPhaseModel::advect
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
        alphaRhos_[fi] = alphaRho_;
        alphaRhoUs_[fi] = alphaRhoU_;
        alphaRhoEs_[fi] = alphaRhoE_;
        alphaRhoPTEs_[fi] = alphaRhoPTE_;
    }

    volScalarField deltaAlphaRho(-fvc::div(massFlux_));
    volVectorField deltaAlphaRhoU
    (
      - fvc::div(momentumFlux_)
      + F
      - (*this)*otherPhase().gradp()
      + alphaRho_*g
    );
    volScalarField deltaAlphaRhoE(-fvc::div(energyFlux_));
    volScalarField deltaAlphaRhoPTE(-fvc::div(PTEFlux_) - Ps_*fvc::div(phi()));

    if (di != -1)
    {
        deltaAlphaRhos_[di] = deltaAlphaRho;
        deltaAlphaRhoUs_[di] = deltaAlphaRhoU;
        deltaAlphaRhoEs_[di] = deltaAlphaRhoE;
        deltaAlphaRhoPTEs_[di] = deltaAlphaRhoPTE;
    }

    volScalarField alphaRho(alphaRho_*coeffs[stepi]);
    volVectorField alphaRhoU(alphaRhoU_*coeffs[stepi]);
    volScalarField alphaRhoE(alphaRhoE_*coeffs[stepi]);
    volScalarField alphaRhoPTE(alphaRhoPTE_*coeffs[stepi]);

    deltaAlphaRho *= Fcoeffs[stepi];
    deltaAlphaRhoU *= Fcoeffs[stepi];
    deltaAlphaRhoE *= Fcoeffs[stepi];
    deltaAlphaRhoPTE *= Fcoeffs[stepi];

    label fieldi = 0;
    label deltai = 0;
    for (label i = 0; i < stepi; i++)
    {
        if (this->storedFieldIndexes_[i] != -1)
        {
            alphaRho += alphaRhos_[fieldi]*coeffs[i];
            alphaRhoU += alphaRhoUs_[fieldi]*coeffs[i];
            alphaRhoE += alphaRhoEs_[fieldi]*coeffs[i];
            alphaRhoPTE += alphaRhoPTEs_[fieldi]*coeffs[i];
            fieldi++;
        }

        if (this->storedDeltaIndexes_[i] != -1)
        {
            deltaAlphaRho += deltaAlphaRhos_[deltai]*Fcoeffs[i];
            deltaAlphaRhoU += deltaAlphaRhoUs_[deltai]*Fcoeffs[i];
            deltaAlphaRhoE += deltaAlphaRhoEs_[deltai]*Fcoeffs[i];
            deltaAlphaRhoPTE += deltaAlphaRhoPTEs_[deltai]*Fcoeffs[i];
            deltai++;
        }
    }

    alphaRho_ = alphaRho + deltaT*deltaAlphaRho;
    alphaRho_.correctBoundaryConditions();

    alphaRhoU_ = alphaRhoU + deltaT*deltaAlphaRhoU;
    alphaRhoU_.correctBoundaryConditions();

    alphaRhoE_ = alphaRhoE + deltaT*deltaAlphaRhoE;
    alphaRhoE_.correctBoundaryConditions();

    alphaRhoPTE_ = alphaRhoPTE + deltaT*deltaAlphaRhoPTE;
    alphaRhoPTE_.correctBoundaryConditions();
}


void Foam::granularPhaseModel::updateFluxes()
{
    // calculate fluxes with
    volScalarField Ptot
    (
        IOobject::groupName("Ptot", name()),
        Ps_ + Pfric_
    );
    this->fluxFunction_->updateFluxes
    (
        massFlux_,
        momentumFlux_,
        energyFlux_,
        PTEFlux_,
        *this,
        rho_,
        U_,
        E_,
        Theta_,
        Ptot,
        c(),
        fluid_.p()
    );
    gradAlpha_ = fluxFunction_->gradAlpha();
}

void Foam::granularPhaseModel::updateFluxes(const surfaceScalarField& alphaf)
{
    NotImplemented;
}


void Foam::granularPhaseModel::decode()
{
    volScalarField& alpha = *this;
    alpha = alphaRho_/rho_;
    alpha.max(residualAlpha_);

    U_ = alphaRhoU_/alphaRho_;
    E_ = alphaRhoE_/this->alphaRho_;
    he_ = E_;
    he_.correctBoundaryConditions();
    thermoPtr_->correct();

    // Update kinetic theory quantities
    Theta_ = 2.0/3.0*alphaRhoPTE_/alphaRho_;
    Theta_.max(0);
    Theta_.min(100);
    Theta_.correctBoundaryConditions();

    volScalarField ThetaSqrt(sqrt(Theta_));
    scalar sqrtPi = sqrt(Foam::constant::mathematical::pi);
    volScalarField da(dPtr_->d());

    gs0_ = radialModel_->g0(*this, alphaMinFriction_, alphaMax_);
    kappa_ = conductivityModel_->kappa(*this, Theta_, gs0_, rho_, da, e_);
    nut_ = viscosityModel_->nu(*this, Theta_, gs0_, rho_, da, e_);
    lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

    // Frictional pressure
    Pfric_ =
        frictionalStressModel_->frictionalPressure
        (
            *this,
            alphaMinFriction_,
            alphaMax_
        );
    Pfric_.correctBoundaryConditions();

    volSymmTensorField D(symm(gradU()));
    nuFric_ = frictionalStressModel_->nu
    (
        *this,
        alphaMinFriction_,
        alphaMax_,
        Pfric_/rho_,
        D
    );

    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_);
    nuFric_ = Foam::min(nuFric_, maxNut_ - nut_);
    nut_ += nuFric_;

    Ps_ =
        granularPressureModel_->granularPressureCoeff
        (
            *this,
            gs0_,
            rho_,
            e_
        )*Theta_;
    Ps_.correctBoundaryConditions();
}


void Foam::granularPhaseModel::encode()
{
    this->max(residualAlpha_);
    this->min(alphaMax_);
    this->correctBoundaryConditions();

    alphaRho_ = (*this)*rho_;
    alphaRho_.correctBoundaryConditions();

    alphaRhoU_ = alphaRho_*U_;
    alphaRhoU_.correctBoundaryConditions();

    E_ = he_;
    alphaRhoE_ = alphaRho_*E_;
    alphaRhoE_.correctBoundaryConditions();

    Theta_.correctBoundaryConditions();
    alphaRhoPTE_ = 1.5*alphaRho_*Theta_;
    alphaRhoPTE_.correctBoundaryConditions();
}


void Foam::granularPhaseModel::correctThermo()
{
    E_ = he_;
    thermoPtr_->correct();

    // Local references
    volScalarField alpha(*this);
    alpha.max(residualAlpha_);

    volScalarField ThetaSqrt(sqrt(Theta_));
    scalar sqrtPi = sqrt(Foam::constant::mathematical::pi);
    volScalarField da(dPtr_->d());

    gs0_ = radialModel_->g0(*this, alphaMinFriction_, alphaMax_);
    kappa_ = conductivityModel_->kappa(*this, Theta_, gs0_, rho_, da, e_);
    nut_ = viscosityModel_->nu(*this, Theta_, gs0_, rho_, da, e_);
    lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

    // Frictional pressure
    Pfric_ =
        frictionalStressModel_->frictionalPressure
        (
            *this,
            alphaMinFriction_,
            alphaMax_
        );
    Pfric_.correctBoundaryConditions();

    volSymmTensorField D(symm(gradU()));
    nuFric_ = frictionalStressModel_->nu
    (
        *this,
        alphaMinFriction_,
        alphaMax_,
        Pfric_/rho_,
        D
    );

    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_);
    nuFric_ = Foam::min(nuFric_, maxNut_ - nut_);
    nut_ += nuFric_;

    Ps_ =
        granularPressureModel_->granularPressureCoeff
        (
            *this,
            gs0_,
            rho_,
            e_
        )*Theta_;
    Ps_.correctBoundaryConditions();
}


bool Foam::granularPhaseModel::read(const dictionary& phaseProperties)
{
    e_.readIfPresent(phaseDict_);
    alphaMax_.readIfPresent(phaseDict_);
    alphaMinFriction_.readIfPresent(phaseDict_);

    viscosityModel_->read();
    conductivityModel_->read();
    radialModel_->read();
    granularPressureModel_->read();
    frictionalStressModel_->read();

    return true;
}


void Foam::granularPhaseModel::store()
{
    (*this).storeOldTime();
    alphaRho_.storeOldTime();
    alphaRhoU_.storeOldTime();
    alphaRhoE_.storeOldTime();
    alphaRhoPTE_.storeOldTime();
    Theta_.storeOldTime();
}

// ************************************************************************* //
