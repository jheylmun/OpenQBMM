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

#include "vdfPhaseModel.H"
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
        vdfPhaseModel,
        dictionary,
        vdf
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vdfPhaseModel::vdfPhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    phaseModel(fluid, phaseProperties, phaseName),
    populationBalanceProperties_
    (
        IOobject
        (
            "populationBalanceProperties",
            fluid.mesh().time().constant(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("0", sqr(dimVelocity), 0.0)
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
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        radialModel_->g0(*this, alphaMinFriction_, alphaMax_)
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
    c_
    (
        IOobject
        (
            IOobject::groupName("c", phaseName),
            fluid.time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimVelocity, 0.0)
    ),
    alphaRhoU_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoU", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        (*this)*rho_*U_,
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
        (*this)*rho_*E_,
        E_.boundaryField()
    ),
    phi_
    (
        IOobject
        (
            IOobject::groupName("phi", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U_)
    ),
    massFlux_
    (
        IOobject
        (
            IOobject::groupName("massFlux", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        phi_*fvc::interpolate((*this)*rho())
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
    gradP_
    (
        IOobject
        (
            IOobject::groupName("gradP", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fvc::grad(Ps_)
    ),
    populationBalance_
    (
        populationBalanceModel::New
        (
            phaseName, populationBalanceProperties_, phi_
        )
    ),
    quadrature_
    (
        fluid.mesh().lookupObjectRef<velocityQuadratureApproximation>
        (
            IOobject::groupName
            (
                "quadratureProperties",
                phaseName
            )
        )
    )
{
    setTurbulenceModel();

    alphaRhoUs_.setSize(quadrature_.nodes().size());
    forAll(alphaRhoUs_, nodei)
    {
        alphaRhoUs_.set
        (
            nodei,
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("alphaRhoU" + Foam::name(nodei), name()),
                    fluid.mesh().time().timeName(),
                    fluid.mesh()
                ),
                quadrature_.nodes()[nodei].primaryWeight()*rho_
               *quadrature_.nodes()[nodei].primaryAbscissa(),
                alphaRhoU_.boundaryField()
            )
        );
    }

    // Kinetic energy is not included in E
    E_ = he_;

    decode();
    updateFluxes();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vdfPhaseModel::~vdfPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vdfPhaseModel::setNSteps
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

            moments_.append
            (
                new PtrList<volScalarField>(quadrature_.nMoments())
            );
            forAll(moments_[currFieldIndex], mi)
            {
                moments_[currFieldIndex].set
                (
                    mi,
                    new volScalarField
                    (
                        IOobject
                        (
                            quadrature_.moments()[mi].name() + Foam::name(stepi),
                            fluid_.mesh().time().timeName(),
                            fluid_.mesh()
                        ),
                        fluid_.mesh(),
                        dimensionedScalar
                        (
                            "0",
                            quadrature_.moments()[mi].dimensions(),
                            0.0
                        )
                    )
                );
            }
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

            currFieldIndex++;
        }

        if (storeDeltas[stepi])
        {
            this->storedDeltaIndexes_[stepi] = currDeltaIndex;

            deltaMoments_.append
            (
                new PtrList<volScalarField>(quadrature_.nMoments())
            );
            forAll(deltaMoments_[currDeltaIndex], mi)
            {
                deltaMoments_[currDeltaIndex].set
                (
                    mi,
                    new volScalarField
                    (
                        IOobject
                        (
                            "delta" + quadrature_.moments()[mi].name()
                          + Foam::name(stepi),
                            fluid_.mesh().time().timeName(),
                            fluid_.mesh()
                        ),
                        fluid_.mesh(),
                        dimensionedScalar
                        (
                            "0",
                            quadrature_.moments()[mi].dimensions()/dimTime,
                            0.0
                        )
                    )
                );
            }
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

            currDeltaIndex++;
        }
    }
}

Foam::tmp<Foam::volScalarField> Foam::vdfPhaseModel::c() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("c", name()),
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
                )*Theta_/sqr(Foam::max((*this)*rho_, residualAlphaRho()))
            )
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::vdfPhaseModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField> Foam::vdfPhaseModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField> Foam::vdfPhaseModel::R() const
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


Foam::tmp<Foam::volScalarField> Foam::vdfPhaseModel::pPrime() const
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


Foam::tmp<Foam::surfaceScalarField> Foam::vdfPhaseModel::pPrimef() const
{
    return fvc::interpolate(pPrime());
}


Foam::tmp<Foam::volSymmTensorField>
Foam::vdfPhaseModel::devRhoReff() const
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


Foam::tmp<Foam::fvVectorMatrix> Foam::vdfPhaseModel::divDevRhoReff
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
Foam::vdfPhaseModel::dissipationCoeff() const
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
Foam::vdfPhaseModel::productionCoeff() const
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


void Foam::vdfPhaseModel::advect
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
        forAll(moments_[fi], mi)
        {
            moments_[fi][mi] = quadrature_.moments()[mi];
        }
        alphaRhoEs_[fi] = alphaRhoE_;
    }

    PtrList<volScalarField> deltaMoments(quadrature_.nMoments());
    forAll(deltaMoments, mi)
    {
        deltaMoments.set
        (
            mi,
            new volScalarField
            (
                -fvc::div(populationBalance_->momentFluxes()[mi])
            )
        );

        // Add gravitational source
        labelList order = quadrature_.moments()[mi].cmptOrders();
        forAll(quadrature_.momentOrders()[mi], cmpt)
        {
            if (order[cmpt] > 0)
            {
                labelList srcOrder = order;
                srcOrder[cmpt] = srcOrder[cmpt] - 1;
                deltaMoments[mi] +=
                    quadrature_.moments()(srcOrder)*g.component(cmpt);
            }
        }
    }
    volScalarField deltaAlphaRhoE(-fvc::div(energyFlux_));

    if (di != -1)
    {
        forAll(deltaMoments, mi)
        {
            deltaMoments_[di][mi] = deltaMoments[mi];
        }
        deltaAlphaRhoEs_[di] = deltaAlphaRhoE;
    }

    PtrList<volScalarField> moments(quadrature_.nMoments());
    forAll(moments, mi)
    {
        moments.set
        (
            mi,
            new volScalarField(quadrature_.moments()[mi]*coeffs[stepi])
        );
        deltaMoments[mi] *= Fcoeffs[stepi];
    }
    volScalarField alphaRhoE(alphaRhoE_*coeffs[stepi]);
    deltaAlphaRhoE *= Fcoeffs[stepi];

    label fieldi = 0;
    label deltai = 0;
    for (label i = 0; i < stepi; i++)
    {
        if (this->storedFieldIndexes_[i] != -1)
        {
            forAll(moments, mi)
            {
                moments[mi] += moments_[fieldi][mi]*coeffs[i];
            }
            alphaRhoE += alphaRhoEs_[fieldi]*coeffs[i];
            fieldi++;
        }

        if (this->storedDeltaIndexes_[i] != -1)
        {
            forAll(deltaMoments, mi)
            {
                deltaMoments[mi] += deltaMoments_[deltai][mi]*Fcoeffs[i];
            }
            deltaAlphaRhoE += deltaAlphaRhoEs_[deltai]*Fcoeffs[i];
            deltai++;
        }
    }

    forAll(moments, mi)
    {
        quadrature_.moments()[mi] == moments[mi] + deltaT*deltaMoments[mi];
        quadrature_.moments()[mi].correctBoundaryConditions();
    }

    alphaRhoE_ = alphaRhoE + deltaT*deltaAlphaRhoE;
    alphaRhoE_.correctBoundaryConditions();
}


void Foam::vdfPhaseModel::updateFluxes()
{
    // calculate fluxes with
    c_ = this->c();
    populationBalance_->updateAdvection();
}

void Foam::vdfPhaseModel::updateFluxes(const surfaceScalarField& alphaf)
{
    NotImplemented;
}


void Foam::vdfPhaseModel::decode()
{
    quadrature_.updateQuadrature();
    quadrature_.moments()[0].max(residualAlpha_);
    quadrature_.moments()[0].correctBoundaryConditions();

    volScalarField& alpha = *this;
    alpha = quadrature_.moments()[0];
    volScalarField alphaRho(alpha*rho_);

    const labelListList& momentOrders = quadrature_.momentOrders();
    label nDims = momentOrders[0].size();
    Theta_ = dimensionedScalar("0", sqr(dimVelocity), 0.0);
    forAll(momentOrders[0], cmpti)
    {
        labelList order1(nDims, 0);
        labelList order2(nDims, 0);
        order1[cmpti] = 1;
        order2[cmpti] = 2;

        alphaRhoU_.replace(cmpti, quadrature_.moments()(order1)*rho_);
        volScalarField meanU("meanU", quadrature_.moments()(order1)/alpha);
        U_.replace(cmpti, meanU);

        Theta_ +=
            Foam::max
            (
                quadrature_.moments()(order2)/alpha - sqr(meanU),
                dimensionedScalar("zero", sqr(dimVelocity), 0.0)
            );
    }
    Theta_ /= 3.0;
    phi_ = fvc::flux(U_);

    forAll(alphaRhoUs_, nodei)
    {
        alphaRhoUs_[nodei] =
            quadrature_.nodes()[nodei].primaryWeight()*rho_
           *quadrature_.nodes()[nodei].primaryAbscissa();
    }

    E_ = alphaRhoE_/alphaRho;
    he_ = E_;
    he_.correctBoundaryConditions();
    thermoPtr_->correct();

    gs0_ = radialModel_->g0(*this, alphaMinFriction_, alphaMax_);

    // Granular pressure
    Ps_ =
        frictionalStressModel_->frictionalPressure
        (
            *this,
            alphaMinFriction_,
            alphaMax_
        );
//       + granularPressureModel_->granularPressureCoeff
//         (
//             *this,
//             gs0_,
//             rho_,
//             e_
//         )*Theta_;
}


void Foam::vdfPhaseModel::encode()
{
    quadrature_.moments()[0].max(residualAlpha_);
    quadrature_.moments()[0].correctBoundaryConditions();
    volScalarField& alpha = *this;
    alpha = quadrature_.moments()[0];

    const labelListList& momentOrders = quadrature_.momentOrders();
    forAll(momentOrders[0], cmpti)
    {
        labelList order1(momentOrders[0].size(), 0);
        order1[cmpti] = 1;

        alphaRhoU_.replace(cmpti, quadrature_.moments()(order1)*rho_);
    }

    forAll(alphaRhoUs_, nodei)
    {
        alphaRhoUs_[nodei] =
            quadrature_.nodes()[nodei].primaryWeight()*rho_
           *quadrature_.nodes()[nodei].primaryAbscissa();
    }

    E_ = he_;
    alphaRhoE_ = alpha*rho_*E_;
    alphaRhoE_.correctBoundaryConditions();
}


void Foam::vdfPhaseModel::correctThermo()
{
    // Local references
    forAll(alphaRhoUs_, nodei)
    {
        quadrature_.nodes()[nodei].primaryAbscissa() =
            alphaRhoUs_[nodei]
           /Foam::max
            (
                quadrature_.nodes()[nodei].primaryWeight()*rho_,
                residualAlphaRho()
            );
    }
    quadrature_.updateMoments();
    const volScalarField& alpha = *this;

    // Granular pressure
    gs0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

    Ps_ =
        frictionalStressModel_->frictionalPressure
        (
            *this,
            alphaMinFriction_,
            alphaMax_
        );
//       + granularPressureModel_->granularPressureCoeff
//         (
//             alpha,
//             gs0_,
//             rho_,
//             e_
//         )*Theta_;

    //- Solve collisions
    populationBalance_->solveSources();

    Theta_ = dimensionedScalar("0", sqr(dimVelocity), 0.0);
    const labelListList& momentOrders = quadrature_.momentOrders();
    label nDims = momentOrders[0].size();
    forAll(momentOrders[0], cmpti)
    {
        labelList order1(nDims, 0);
        labelList order2(nDims, 0);
        order1[cmpti] = 1;
        order2[cmpti] = 2;

        volScalarField meanU("meanU", quadrature_.moments()(order1)/alpha);
        Theta_ +=
            Foam::max
            (
                quadrature_.moments()(order2)/alpha - sqr(meanU),
                dimensionedScalar("zero", sqr(dimVelocity), 0.0)
            );
    }
    Theta_ /= 3.0;
}


bool Foam::vdfPhaseModel::read(const dictionary& phaseProperties)
{
    e_.readIfPresent(phaseDict_);
    alphaMax_.readIfPresent(phaseDict_);
    alphaMinFriction_.readIfPresent(phaseDict_);

//     granularPressureModel_->read();
    frictionalStressModel_->read();

    return true;
}

// ************************************************************************* //
