/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "twoPhaseSystem.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "BlendedInterfacialModel.H"
#include "virtualMassModel.H"
#include "heatTransferModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"
#include "fvMatrix.H"
#include "surfaceInterpolate.H"
#include "MULES.H"
#include "subCycle.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcCurl.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fixedValueFvsPatchFields.H"
#include "blendingMethod.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseSystem::relaxPressure()
{
    if
    (
        phase1_->granular() || phase1_->slavePressure()
     || phase2_->granular() || phase2_->slavePressure()
    )
    {
        return;
    }

    dimensionedScalar localT("localT", dimTime, 0);
    bool timeComplete = false;

    volScalarField& alpha1(phase1_());
    volScalarField& alpha2(phase2_());
    volScalarField& alphaRhoE1(phase1_->alphaRhoE());
    volScalarField& alphaRhoE2(phase2_->alphaRhoE());

    volScalarField& p1(phase1_->p());
    volScalarField& p2(phase2_->p());

    while (!timeComplete)
    {
        encode();

        // Store old fields
        volScalarField alpha1Old(alpha1);
        volScalarField alphaRhoE1Old(alphaRhoE1);
        volScalarField alphaRhoE2Old(alphaRhoE2);

        // 1st Predictor
        volScalarField theta(alpha1*alpha2/(tauP_*(p1 + p2)));
        volScalarField alphaK1(pDt_*theta*(p2 - p1));
        volScalarField alphaRhoEK1(pDt_*p_*theta*(p2 - p1));

        alpha1 -= 0.5*alphaK1;
        alpha2 = 1.0 - alpha1;

        alphaRhoE1 += 0.5*alphaRhoEK1;
        alphaRhoE2 -= 0.5*alphaRhoEK1;

        decode();

        // 2nd Predictor
        encode();

        theta = alpha1*alpha2/(tauP_*(p1 + p2));
        volScalarField alphaK2(pDt_*theta*(p2 - p1));
        volScalarField alphaRhoEK2(pDt_*p_*theta*(p2 - p1));

        alpha1 -= 0.5*alphaK2;
        alpha2 = 1.0 - alpha1;

        alphaRhoE1 += 0.5*alphaRhoEK2;
        alphaRhoE2 -= 0.5*alphaRhoEK2;

        decode();

        // 3rd Predictor
        encode();

        theta = alpha1*alpha2/(tauP_*(p1 + p2));
        volScalarField alphaK3(pDt_*theta*(p2 - p1));
        volScalarField alphaRhoEK3(pDt_*p_*theta*(p2 - p1));

        alpha1 = alpha1Old - alphaK3;
        alpha2 = 1.0 - alpha1;

        alphaRhoE1 = alphaRhoE1Old + alphaRhoEK3;
        alphaRhoE2 = alphaRhoE2Old - alphaRhoEK3;

        decode();

        volScalarField alphaNew(alpha1);
        volScalarField pNew(p_);

        // Corrector
        encode();

        theta = alpha1*alpha2/(tauP_*(p1 + p2));
        volScalarField alphaK4
        (
            (
                alphaK1 + 2.0*alphaK2 + 2.0*alphaK3
            + pDt_*theta*(p2 - p1)
            )/6.0
        );
        volScalarField alphaRhoEK4
        (
            (
                alphaRhoEK1 + 2.0*alphaRhoEK2 + 2.0*alphaRhoEK3
            + pDt_*p_*theta*(p2 - p1)
            )/6.0
        );

        alpha1 = alpha1Old - alphaK4;
        alpha2 = 1.0 - alpha1;

        alphaRhoE1 = alphaRhoE1Old + alphaRhoEK4;
        alphaRhoE2 = alphaRhoE2Old - alphaRhoEK4;

        decode();

        scalar error =
            sqr
            (
                max(mag(pNew - p_)/(pTol_ + relTol_*max(p_, pNew)))
            ).value();
        error +=
            sqr
            (
                max
                (
                    mag(alphaNew - alpha1)
                /(alphaTol_ + relTol_*max(alpha1, alphaNew))
                )
            ).value();
        error = sqrt(error/2.0) + SMALL;

        if (error < 1)
        {
            pDt_ *= min(facMax_, max(facMin_, fac_/Foam::pow(error, 1.0/3.0)));
            dimensionedScalar tmpDt(pDt_);

            dimensionedScalar maxLocalDt =
                max
                (
                    mesh_.time().deltaT() - localT,
                    dimensionedScalar("0", dimTime, 0.0)
                );
            pDt_ = min(maxLocalDt, pDt_);
            if (pDt_.value() == 0.0)
            {
                timeComplete = true;
                pDt_ = tmpDt;
            }
            localT += pDt_;
        }
        else
        {
            pDt_ *= min(1, max(facMin_, fac_/Foam::pow(error, 1.0/3.0)));

            alpha1 = alpha1Old;
            alphaRhoE1 = alphaRhoE1Old;
            alphaRhoE2 = alphaRhoE2Old;
            decode();
        }
    }
    alpha1.correctBoundaryConditions();
    alpha2.correctBoundaryConditions();
    alphaRhoE1.correctBoundaryConditions();
    alphaRhoE2.correctBoundaryConditions();
}


void Foam::twoPhaseSystem::relaxVelocity()
{

    if (instantRelaxation_  || phase1_->granular() || phase2_->granular())
    {
        return;
    }

    dimensionedScalar localT("localT", dimTime, 0);
    label nItt = 0;
    bool timeComplete = false;

    volVectorField& alphaRhoU1(phase1_->alphaRhoU());
    volVectorField& alphaRhoU2(phase2_->alphaRhoU());
    volScalarField& alphaRhoE1(phase1_->alphaRhoE());
    volScalarField& alphaRhoE2(phase2_->alphaRhoE());

    volVectorField& U1(phase1_->U());
    volVectorField& U2(phase2_->U());

    while (!timeComplete)
    {
        encode();
        volVectorField alphaRhoU1Old(alphaRhoU1);
        volVectorField alphaRhoU2Old(alphaRhoU2);
        volScalarField alphaRhoE1Old(alphaRhoE1);
        volScalarField alphaRhoE2Old(alphaRhoE2);

        //- 1st predictor
        volScalarField K(Kd());
        volVectorField UK1(uDt_*K*(U2 - U1));
        volScalarField EK1(uDt_*K*(U_ & (U2 - U1)));

        alphaRhoU1 += 0.5*UK1;
        alphaRhoU2 -= 0.5*UK1;
        alphaRhoE1 += 0.5*EK1;
        alphaRhoE2 -= 0.5*EK1;

        decode();
        U_ = mixtureU();

        //- 2nd predictor
        encode();

        volVectorField UK2(uDt_*K*(U2 - U1));
        volScalarField EK2(uDt_*K*(U_ & (U2 - U1)));

        alphaRhoU1 += 0.5*UK2;
        alphaRhoU2 -= 0.5*UK2;
        alphaRhoE1 += 0.5*EK2;
        alphaRhoE2 -= 0.5*EK2;

        decode();
        U_ = mixtureU();

        //- 3rd predictor
        encode();
        volVectorField UK3(uDt_*K*(U2 - U1));
        volScalarField EK3(uDt_*K*(U_ & (U2 - U1)));

        alphaRhoU1 = alphaRhoU1Old + UK3;
        alphaRhoU2 = alphaRhoU2Old - UK3;
        alphaRhoE1 = alphaRhoE1Old + EK3;
        alphaRhoE2 = alphaRhoE2Old - EK3;

        decode();
        volVectorField U1New(U1);
        volVectorField U2New(U2);
        U_ = mixtureU();

        //- Corrector
        encode();
        volVectorField UK4
        (
            (UK1 + 2.0*UK2 + 2.0*UK3 + uDt_*K*(U2 - U1))/6.0
        );
        volScalarField EK4
        (
            (EK1 + 2.0*EK2 + 2.0*EK3 + uDt_*K*(U_ & (U2 - U1)))/6.0
        );

        alphaRhoU1 = alphaRhoU1Old + UK4;
        alphaRhoU2 = alphaRhoU2Old - UK4;
        alphaRhoE1 = alphaRhoE1Old + EK4;
        alphaRhoE2 = alphaRhoE2Old - EK4;

        decode();
        U_ = mixtureU();

        scalar error = 0.0;
        error +=
            sqr
            (
                max(mag(U1New - U1)/(uTol_ + relTol_*mag(U1)))
            ).value();
        error +=
            sqr
            (
                max(mag(U2New - U2)/(uTol_ + relTol_*mag(U2)))
            ).value();
        error = Foam::sqrt(error/2.0) + SMALL;

        if (error < 1)
        {
            uDt_ *= min(facMax_, max(facMin_, fac_/Foam::pow(error, 1.0/3.0)));
            dimensionedScalar tmpDt(uDt_);
            dimensionedScalar maxLocalDt =
                max
                (
                    mesh_.time().deltaT() - localT,
                    dimensionedScalar("0", dimTime, 0.0)
                );
            uDt_ = min(maxLocalDt, uDt_);
            if (uDt_.value() == 0.0)
            {
                timeComplete = true;
                uDt_ = tmpDt;
            }
            localT += uDt_;
        }
        else
        {
            uDt_ *= min(1, max(facMin_, fac_/Foam::pow(error, 1.0/3.0)));

            alphaRhoU1 = alphaRhoU1Old;
            alphaRhoU2 = alphaRhoU2Old;
            alphaRhoE1 = alphaRhoE1Old;
            alphaRhoE2 = alphaRhoE2Old;
            decode();
        }
        nItt++;
    }
    alphaRhoU1.correctBoundaryConditions();
    alphaRhoU2.correctBoundaryConditions();
    alphaRhoE1.correctBoundaryConditions();
    alphaRhoE2.correctBoundaryConditions();

    Info<< "number of velocity relaxation iterations: " << nItt << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::twoPhaseSystem
(
    const fvMesh& mesh,
    const dimensionedVector& g
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    phase1_
    (
        phaseModel::New
        (
            *this,
            *this,
            wordList(lookup("phases"))[0]
        )
    ),

    phase2_
    (
        phaseModel::New
        (
            *this,
            *this,
            wordList(lookup("phases"))[1]
        )
    ),

    U_
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (
            phase1_->alphaRho()*phase1_->U()
          + phase2_->alphaRho()*phase2_->U()
        )/(phase1_->alphaRho() + phase2_->alphaRho())
    ),
    p_
    (
        IOobject
        (
            "pMixture",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase1_->p()
    ),
    g_(g),
    fluxIntegrator_
    (
        phaseFluxIntegrator::New(phase1_(), phase2_())
    ),
    relaxationDict_(subDict("relaxation")),
    instantRelaxation_(relaxationDict_.lookup("instantRelaxation")),
    uDt_(mesh_.time().deltaT()/100),
    pDt_(uDt_),
    maxItt_(5000),
    fac_((relaxationDict_.lookupOrDefault("fac", 1.0))),
    facMin_((relaxationDict_.lookupOrDefault("facMin", 0.5))),
    facMax_((relaxationDict_.lookupOrDefault("facMax", 2.0))),
    relTol_((relaxationDict_.lookupOrDefault("relTol", 1e-6))),
    pTol_
    (
        dimensionedScalar::lookupOrDefault
        (
            "pTol",
            relaxationDict_,
            dimPressure,
            1.0
        )
    ),
    uTol_
    (
        dimensionedScalar::lookupOrDefault
        (
            "uTol",
            relaxationDict_,
            dimVelocity,
            1e-4
        )
    ),
    alphaTol_
    (
        dimensionedScalar::lookupOrDefault("alphaTol", relaxationDict_, 1e-8)
    ),
    tauP_("pressureRelaxationTime", dimTime, relaxationDict_)
{
    volScalarField& alpha2 = phase2_();
    alpha2 = (scalar(1) - phase1_());
    encode();

    // Blending
    forAllConstIter(dictionary, subDict("blending"), iter)
    {
        blendingMethods_.insert
        (
            iter().dict().dictName(),
            blendingMethod::New
            (
                iter().dict(),
                wordList(lookup("phases"))
            )
        );
    }

    // Pairs
    phasePair::scalarTable sigmaTable(lookup("sigma"));
    phasePair::dictTable aspectRatioTable(lookup("aspectRatio"));

    pair_.set
    (
        new phasePair
        (
            phase1_(),
            phase2_(),
            g,
            sigmaTable
        )
    );

    pair1In2_.set
    (
        new orderedPhasePair
        (
            phase1_(),
            phase2_(),
            g,
            sigmaTable,
            aspectRatioTable
        )
    );

    pair2In1_.set
    (
        new orderedPhasePair
        (
            phase2_(),
            phase1_(),
            g,
            sigmaTable,
            aspectRatioTable
        )
    );


    // Models

    drag_.set
    (
        new BlendedInterfacialModel<dragModel>
        (
            lookup("drag"),
            (
                blendingMethods_.found("drag")
              ? blendingMethods_["drag"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_,
            false // Do not zero drag coefficent at fixed-flux BCs
        )
    );

    virtualMass_.set
    (
        new BlendedInterfacialModel<virtualMassModel>
        (
            lookup("virtualMass"),
            (
                blendingMethods_.found("virtualMass")
              ? blendingMethods_["virtualMass"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    heatTransfer_.set
    (
        new BlendedInterfacialModel<heatTransferModel>
        (
            lookup("heatTransfer"),
            (
                blendingMethods_.found("heatTransfer")
              ? blendingMethods_["heatTransfer"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    lift_.set
    (
        new BlendedInterfacialModel<liftModel>
        (
            lookup("lift"),
            (
                blendingMethods_.found("lift")
              ? blendingMethods_["lift"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    wallLubrication_.set
    (
        new BlendedInterfacialModel<wallLubricationModel>
        (
            lookup("wallLubrication"),
            (
                blendingMethods_.found("wallLubrication")
              ? blendingMethods_["wallLubrication"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    turbulentDispersion_.set
    (
        new BlendedInterfacialModel<turbulentDispersionModel>
        (
            lookup("turbulentDispersion"),
            (
                blendingMethods_.found("turbulentDispersion")
              ? blendingMethods_["turbulentDispersion"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    if (phase1_->granular() || phase1_->slavePressure())
    {
        p_ = phase2_->p();
    }
    else if (phase2_->granular() || phase2_->slavePressure())
    {
        p_ = phase1_->p();
    }
    else
    {
        p_ =
        (
            phase1_->alphaRho()*phase1_->p()
          + phase2_->alphaRho()*phase2_->p()
        )/(phase1_->alphaRho() + phase2_->alphaRho());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::~twoPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::rho() const
{
    return phase1_().alphaRho() + phase2_().alphaRho();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Kd() const
{
    return drag_->K();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::Kdf() const
{
    return drag_->Kf();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Vm() const
{
    return virtualMass_->K();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::Vmf() const
{
    return virtualMass_->Kf();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Kh() const
{
    return heatTransfer_->K();
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseSystem::F() const
{
    return lift_->F<vector>() + wallLubrication_->F<vector>();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::Ff() const
{
    return lift_->Ff() + wallLubrication_->Ff();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::D() const
{
    return turbulentDispersion_->D();
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseSystem::mixtureU() const
{
    return
    (
        phase1_->alphaRho()*phase1_->U() + phase2_->alphaRho()*phase2_->U()
    )/(phase1_->alphaRho() + phase2_->alphaRho());
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::mixturep() const
{
    if (phase1_->granular() || phase1_->slavePressure())
    {
        return phase2_->p();
    }
    else if (phase2_->granular() || phase2_->slavePressure())
    {
        return phase1_->p();
    }

    return
    (
        phase1_->alphaRho()*phase1_->p() + phase2_->alphaRho()*phase2_->p()
    )/(phase1_->alphaRho() + phase2_->alphaRho());
}


void Foam::twoPhaseSystem::advect()
{
    fluxIntegrator_->integrateFluxes(g_, U_, p_);
}


void Foam::twoPhaseSystem::solve()
{
}


void Foam::twoPhaseSystem::correct()
{
    phase1_->correct();
    phase2_->correct();
}


void Foam::twoPhaseSystem::correctTurbulence()
{
    phase1_->turbulence().correct();
    phase2_->turbulence().correct();
}


void Foam::twoPhaseSystem::decode()
{
    phase1_->decode();
    phase2_->decode();

    p_ = mixturep();
    U_ = mixtureU();
}


void Foam::twoPhaseSystem::encode()
{
    phase1_->encode();
    phase2_->encode();
}


void Foam::twoPhaseSystem::relax()
{
//     relaxPressure();
//     relaxVelocity();

    dimensionedScalar deltaT(mesh_.time().deltaT());

    // Momentum and heat transfer
    volScalarField XiD(1.0/phase1_->alphaRho() + 1.0/phase2_->alphaRho());
    volScalarField XiE
    (
        1.0/(phase1_->alphaRho()*phase1_->Cv())
      + 1.0/(phase2_->alphaRho()*phase2_->Cv())
    );

    volScalarField dragCoeff(Kd());
    volVectorField deltaM
    (
        (phase1_->U() - phase2_->U())/XiD
       *(1.0/(dragCoeff*XiD*deltaT + 1.0) - 1.0)
    );
    volScalarField deltaE
    (
        (phase1_->thermo().T() - phase2_->thermo().T())/XiE
       *(exp(-max(Kh()*XiE*deltaT, vSmall)) - 1.0)
    );
    phase1_->store();
    phase2_->store();

    phase1_->alphaRhoU() += deltaM;
    phase2_->alphaRhoU() -= deltaM;

    phase1_->alphaRhoE() += deltaE;
    phase2_->alphaRhoE() -= deltaE;

    if (phase1_->granular() || phase2_->granular())
    {
        phaseModel* particles;
        phaseModel* gas;
        if (phase1_->granular())
        {
            particles = &phase1_();
            gas = &phase2_();
        }
        else
        {
            particles = &phase2_();
            gas = &phase1_();
        }
        const volScalarField& alphaRhop = particles->alphaRho();
        volScalarField& Theta = particles->Theta();

        volScalarField magSqrUp(magSqr(particles->alphaRhoU()/alphaRhop));

        gas->alphaRhoE() -=
            particles->alphaRho()*0.5
           *(
                magSqrUp - magSqr(particles->U())
            );
        volScalarField ThetaStar
        (
            Theta*exp(-max(2.0*dragCoeff*deltaT/alphaRhop, vSmall))
        );
        ThetaStar.max(vSmall);

        volScalarField XiSlip
        (
            particles->productionCoeff()
           /(dragCoeff*XiD*deltaT + 1.0)
        );

        volScalarField ThetaStarStar
        (
            pow
            (
                XiSlip/alphaRhop*deltaT
              + pow(ThetaStar, 1.5), 2.0/3.0
            )
        );
        gas->alphaRhoE() -= 1.5*alphaRhop*(ThetaStarStar - Theta);

        Theta =
            ThetaStarStar*9.0*sqr(alphaRhop)
           /sqr
            (
                3.0*alphaRhop
              + deltaT*particles->dissipationCoeff()*sqrt(ThetaStarStar)
            );

        particles->alphaRhoPTE() = 1.5*alphaRhop*Theta;
        particles->alphaRhoPTE().correctBoundaryConditions();

        particles->alphaRhoE() -= 1.5*alphaRhop*(Theta - ThetaStarStar);
    }
    else
    {
        volScalarField magSqrU1
        (
            magSqr(phase1_->alphaRhoU()/phase1_->alphaRho())
        );
        volScalarField magSqrU2
        (
            magSqr(phase2_->alphaRhoU()/phase2_->alphaRho())
        );
        phase1_->alphaRhoE() -=
            phase1_->alphaRho()*0.5
           *(
                magSqrU1 - magSqr(phase1_->U())
            );
        phase2_->alphaRhoE() -=
            phase2_->alphaRho()*0.5
           *(
                magSqrU2 - magSqr(phase2_->U())
            );
    }

    phase1_->alphaRhoU().correctBoundaryConditions();
    phase2_->alphaRhoU().correctBoundaryConditions();
    phase1_->alphaRhoE().correctBoundaryConditions();
    phase2_->alphaRhoE().correctBoundaryConditions();

    decode();
    U_ = mixtureU();
    p_ = mixturep();

    relaxPressure();

}


bool Foam::twoPhaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        readOK &= phase1_->read(*this);
        readOK &= phase2_->read(*this);

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


const Foam::dimensionedScalar& Foam::twoPhaseSystem::sigma() const
{
    return pair_->sigma();
}


// ************************************************************************* //
