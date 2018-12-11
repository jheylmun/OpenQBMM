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
//     if
//     (
//         phase1_->granular() || phase1_->slavePressure()
//      || phase2_->granular() || phase2_->slavePressure()
//     )
//     {
//         return;
//     }
//
//     dimensionedScalar localT("localT", dimTime, 0);
//     bool timeComplete = false;
//
//     volScalarField& alpha1(phase1_());
//     volScalarField& alpha2(phase2_());
//     volScalarField& alphaRhoE1(phase1_->alphaRhoERef());
//     volScalarField& alphaRhoE2(phase2_->alphaRhoERef());
//
//     volScalarField& p1(phase1_->pRef());
//     volScalarField& p2(phase2_->pRef());
//
//     while (!timeComplete)
//     {
//         encode();
//
//         // Store old fields
//         volScalarField alpha1Old(alpha1);
//         volScalarField alphaRhoE1Old(alphaRhoE1);
//         volScalarField alphaRhoE2Old(alphaRhoE2);
//
//         // 1st Predictor
//         volScalarField theta(alpha1*alpha2/(tauP_*(p1 + p2)));
//         volScalarField alphaK1(pDt_*theta*(p2 - p1));
//         volScalarField alphaRhoEK1(pDt_*p_*theta*(p2 - p1));
//
//         alpha1 -= 0.5*alphaK1;
//         alpha2 = 1.0 - alpha1;
//
//         alphaRhoE1 += 0.5*alphaRhoEK1;
//         alphaRhoE2 -= 0.5*alphaRhoEK1;
//
//         decode();
//
//         // 2nd Predictor
//         encode();
//
//         theta = alpha1*alpha2/(tauP_*(p1 + p2));
//         volScalarField alphaK2(pDt_*theta*(p2 - p1));
//         volScalarField alphaRhoEK2(pDt_*p_*theta*(p2 - p1));
//
//         alpha1 -= 0.5*alphaK2;
//         alpha2 = 1.0 - alpha1;
//
//         alphaRhoE1 += 0.5*alphaRhoEK2;
//         alphaRhoE2 -= 0.5*alphaRhoEK2;
//
//         decode();
//
//         // 3rd Predictor
//         encode();
//
//         theta = alpha1*alpha2/(tauP_*(p1 + p2));
//         volScalarField alphaK3(pDt_*theta*(p2 - p1));
//         volScalarField alphaRhoEK3(pDt_*p_*theta*(p2 - p1));
//
//         alpha1 = alpha1Old - alphaK3;
//         alpha2 = 1.0 - alpha1;
//
//         alphaRhoE1 = alphaRhoE1Old + alphaRhoEK3;
//         alphaRhoE2 = alphaRhoE2Old - alphaRhoEK3;
//
//         decode();
//
//         volScalarField alphaNew(alpha1);
//         volScalarField pNew(p_);
//
//         // Corrector
//         encode();
//
//         theta = alpha1*alpha2/(tauP_*(p1 + p2));
//         volScalarField alphaK4
//         (
//             (
//                 alphaK1 + 2.0*alphaK2 + 2.0*alphaK3
//             + pDt_*theta*(p2 - p1)
//             )/6.0
//         );
//         volScalarField alphaRhoEK4
//         (
//             (
//                 alphaRhoEK1 + 2.0*alphaRhoEK2 + 2.0*alphaRhoEK3
//             + pDt_*p_*theta*(p2 - p1)
//             )/6.0
//         );
//
//         alpha1 = alpha1Old - alphaK4;
//         alpha2 = 1.0 - alpha1;
//
//         alphaRhoE1 = alphaRhoE1Old + alphaRhoEK4;
//         alphaRhoE2 = alphaRhoE2Old - alphaRhoEK4;
//
//         decode();
//
//         scalar error =
//             sqr
//             (
//                 max(mag(pNew - p_)/(pTol_ + relTol_*max(p_, pNew)))
//             ).value();
//         error +=
//             sqr
//             (
//                 max
//                 (
//                     mag(alphaNew - alpha1)
//                 /(alphaTol_ + relTol_*max(alpha1, alphaNew))
//                 )
//             ).value();
//         error = sqrt(error/2.0) + SMALL;
//
//         if (error < 1)
//         {
//             pDt_ *= min(facMax_, max(facMin_, fac_/Foam::pow(error, 1.0/3.0)));
//             dimensionedScalar tmpDt(pDt_);
//
//             dimensionedScalar maxLocalDt =
//                 max
//                 (
//                     mesh_.time().deltaT() - localT,
//                     dimensionedScalar("0", dimTime, 0.0)
//                 );
//             pDt_ = min(maxLocalDt.value(), pDt_);
//             if (pDt_.value() == 0.0)
//             {
//                 timeComplete = true;
//                 pDt_ = tmpDt;
//             }
//             localT += pDt_;
//         }
//         else
//         {
//             pDt_ *= min(1, max(facMin_, fac_/Foam::pow(error, 1.0/3.0)));
//
//             alpha1 = alpha1Old;
//             alphaRhoE1 = alphaRhoE1Old;
//             alphaRhoE2 = alphaRhoE2Old;
//             decode();
//         }
//     }
//     alpha1.correctBoundaryConditions();
//     alpha2.correctBoundaryConditions();
//     alphaRhoE1.correctBoundaryConditions();
//     alphaRhoE2.correctBoundaryConditions();
}


void Foam::twoPhaseSystem::relaxVelocity()
{

    if (instantRelaxation_  || phase1_->granular() || phase2_->granular())
    {
        return;
    }

    volVectorField& U1(phase1_->URef());
    volVectorField& U2(phase2_->URef());
    volScalarField& E1(phase1_->ERef());
    volScalarField& E2(phase2_->ERef());

    const volScalarField& alphaRho1(phase1_->alphaRho());
    const volScalarField& alphaRho2(phase2_->alphaRho());
    volVectorField& alphaRhoU1(phase1_->alphaRhoURef());
    volVectorField& alphaRhoU2(phase2_->alphaRhoURef());
    volScalarField& alphaRhoE1(phase1_->alphaRhoERef());
    volScalarField& alphaRhoE2(phase2_->alphaRhoERef());


    volScalarField K(Kd());

    forAll(alphaRhoU1, celli)
    {
        scalar localT = 0.0;
        bool timeComplete = false;

        while (!timeComplete)
        {
            vector alphaRhoU1Old(alphaRhoU1[celli]);
            vector alphaRhoU2Old(alphaRhoU2[celli]);
            scalar alphaRhoE1Old(alphaRhoE1[celli]);
            scalar alphaRhoE2Old(alphaRhoE2[celli]);

            //- 1st predictor
            vector UK1 =
                uDt_[celli]*K[celli]*(U2[celli] - U1[celli]);
            scalar EK1 =
                uDt_[celli]*K[celli]*(U_[celli] & (U2[celli] - U1[celli]));

            alphaRhoU1[celli] += 0.5*UK1;
            alphaRhoU2[celli] -= 0.5*UK1;
            alphaRhoE1[celli] += 0.5*EK1;
            alphaRhoE2[celli] -= 0.5*EK1;

            U1[celli] = alphaRhoU1[celli]/alphaRho1[celli];
            U2[celli] = alphaRhoU2[celli]/alphaRho2[celli];
            E1[celli] = alphaRhoE1[celli]/alphaRho1[celli];
            E2[celli] = alphaRhoE2[celli]/alphaRho2[celli];
            U_ =
                (alphaRhoU1[celli] + alphaRhoU2[celli])
               /(alphaRho1[celli] + alphaRho2[celli]);

            //- 2nd predictor
            vector UK2 =
                uDt_[celli]*K[celli]*(U2[celli] - U1[celli]);
            scalar EK2 =
                uDt_[celli]*K[celli]*(U_[celli] & (U2[celli] - U1[celli]));

            alphaRhoU1[celli] += 0.5*UK2;
            alphaRhoU2[celli] -= 0.5*UK2;
            alphaRhoE1[celli] += 0.5*EK2;
            alphaRhoE2[celli] -= 0.5*EK2;

            U1[celli] = alphaRhoU1[celli]/alphaRho1[celli];
            U2[celli] = alphaRhoU2[celli]/alphaRho2[celli];
            E1[celli] = alphaRhoE1[celli]/alphaRho1[celli];
            E2[celli] = alphaRhoE2[celli]/alphaRho2[celli];
            U_ =
                (alphaRhoU1[celli] + alphaRhoU2[celli])
               /(alphaRho1[celli] + alphaRho2[celli]);

            //- 3rd predictor
            vector UK3 =
                uDt_[celli]*K[celli]*(U2[celli] - U1[celli]);
            scalar EK3 =
                uDt_[celli]*K[celli]*(U_[celli] & (U2[celli] - U1[celli]));

            alphaRhoU1[celli] += 0.5*UK3;
            alphaRhoU2[celli] -= 0.5*UK3;
            alphaRhoE1[celli] += 0.5*EK3;
            alphaRhoE2[celli] -= 0.5*EK3;

            U1[celli] = alphaRhoU1[celli]/alphaRho1[celli];
            U2[celli] = alphaRhoU2[celli]/alphaRho2[celli];
            E1[celli] = alphaRhoE1[celli]/alphaRho1[celli];
            E2[celli] = alphaRhoE2[celli]/alphaRho2[celli];
            U_ =
                (alphaRhoU1[celli] + alphaRhoU2[celli])
               /(alphaRho1[celli] + alphaRho2[celli]);

            vector U1New = U1[celli];
            vector U2New = U2[celli];

            //- Corrector
            vector UK4 =
            (
                UK1 + 2.0*UK2 + 2.0*UK3
              + uDt_[celli]*K[celli]*(U2[celli] - U1[celli])
            )/6.0;
            scalar EK4 =
            (
                EK1 + 2.0*EK2 + 2.0*EK3
              + uDt_[celli]*K[celli]*(U_[celli] & (U2[celli] - U1[celli]))
            )/6.0;

            alphaRhoU1[celli] += 0.5*UK4;
            alphaRhoU2[celli] -= 0.5*UK4;
            alphaRhoE1[celli] += 0.5*EK4;
            alphaRhoE2[celli] -= 0.5*EK4;

            U1[celli] = alphaRhoU1[celli]/alphaRho1[celli];
            U2[celli] = alphaRhoU2[celli]/alphaRho2[celli];
            E1[celli] = alphaRhoE1[celli]/alphaRho1[celli];
            E2[celli] = alphaRhoE2[celli]/alphaRho2[celli];
            U_ =
                (alphaRhoU1[celli] + alphaRhoU2[celli])
               /(alphaRho1[celli] + alphaRho2[celli]);

            scalar error = 0.0;
            error +=
                sqr
                (
                    mag(U1New - U1[celli])/(uTol_ + relTol_*mag(U1[celli]))
                );
            error +=
                sqr
                (
                    mag(U2New - U2[celli])/(uTol_ + relTol_*mag(U2[celli]))
                );
            error = Foam::sqrt(error/2.0) + SMALL;

            if (error < 1)
            {
                uDt_[celli] *=
                    min(facMax_, max(facMin_, fac_/Foam::pow(error, 1.0/3.0)));
                scalar maxLocalDt =
                    max
                    (
                        mesh_.time().deltaTValue() - localT,
                        0.0
                    );
                if (maxLocalDt == 0.0)
                {
                    timeComplete = true;
                    localT = 0.0;
                }
                else
                {
                    uDt_[celli] = min(maxLocalDt, uDt_[celli]);
                    localT += uDt_[celli];
                }
            }
            else
            {
                uDt_[celli] *=
                    min(1, max(facMin_, fac_/Foam::pow(error, 1.0/3.0)));

                alphaRhoU1[celli] = alphaRhoU1Old;
                alphaRhoU2[celli] = alphaRhoU2Old;
                alphaRhoE1[celli] = alphaRhoE1Old;
                alphaRhoE2[celli] = alphaRhoE2Old;

                U1[celli] = alphaRhoU1[celli]/alphaRho1[celli];
                U2[celli] = alphaRhoU2[celli]/alphaRho2[celli];
                E1[celli] = alphaRhoE1[celli]/alphaRho1[celli];
                E2[celli] = alphaRhoE2[celli]/alphaRho2[celli];
                U_ =
                    (alphaRhoU1[celli] + alphaRhoU2[celli])
                   /(alphaRho1[celli] + alphaRho2[celli]);

            }
        }
    }
    alphaRhoU1.correctBoundaryConditions();
    alphaRhoU2.correctBoundaryConditions();
    alphaRhoE1.correctBoundaryConditions();
    alphaRhoE2.correctBoundaryConditions();
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
    uDt_(mesh.nCells(), mesh_.time().deltaTValue()),
    pDt_(mesh_.nCells(), mesh_.time().deltaTValue()),
    maxItt_(5000),
    fac_((relaxationDict_.lookupOrDefault("fac", 1.0))),
    facMin_((relaxationDict_.lookupOrDefault("facMin", 0.5))),
    facMax_((relaxationDict_.lookupOrDefault("facMax", 2.0))),
    relTol_((relaxationDict_.lookupOrDefault("relTol", 1e-6))),
    pTol_(relaxationDict_.lookupOrDefault("pTol", 1.0)),
    uTol_(relaxationDict_.lookupOrDefault("uTol", 1e-4)),
    alphaTol_(relaxationDict_.lookupOrDefault("alphaTol", 1e-8)),
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


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Kd
(
    const label nodei,
    const label nodej
) const
{
    return drag_->K(nodei, nodej);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Kdf(const label nodei, const label nodej) const
{
    return drag_->Kf(nodei, nodej);
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::Vm(const label nodei, const label nodej) const
{
    return virtualMass_->K(nodei, nodej);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Vmf(const label nodei, const label nodej) const
{
    return virtualMass_->Kf(nodei, nodej);
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::Kh(const label nodei, const label nodej) const
{
    return heatTransfer_->K(nodei, nodej);
}


Foam::tmp<Foam::volVectorField>
Foam::twoPhaseSystem::F(const label nodei, const label nodej) const
{
    return
        lift_->F<vector>(nodei, nodej)
      + wallLubrication_->F<vector>(nodei, nodej);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Ff(const label nodei, const label nodej) const
{
    return lift_->Ff(nodei, nodej) + wallLubrication_->Ff(nodei, nodej);
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::D(const label nodei, const label nodej) const
{
    return turbulentDispersion_->D(nodei, nodej);
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
    volScalarField XiDMean
    (
        IOobject
        (
            "XiDmean",
            mesh_.time().timeName(),
            mesh_
        ),
        (1.0/phase1_->alphaRho() + 1.0/phase2_->alphaRho())
    );
    volScalarField dragCoeffMean
    (
        IOobject
        (
            "dragCoeffMean",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("0", dimDensity/dimTime, 0.0)
    );

    for (label nodei = 0; nodei < phase1_->nNodes(); nodei++)
    {
        for (label nodej = 0; nodej < phase2_->nNodes(); nodej++)
        {
            volScalarField alphaRho1
            (
                max(phase1_->alphaRho(nodei), phase1_->residualAlphaRho())
            );
            volScalarField alphaRho2
            (
                max(phase2_->alphaRho(nodej), phase2_->residualAlphaRho())
            );
            // Momentum and heat transfer
            volScalarField XiD
            (
                1.0/alphaRho1 + 1.0/alphaRho2
            );
            volScalarField XiE
            (
                1.0/(alphaRho1*phase1_->Cv())
              + 1.0/(alphaRho2*phase2_->Cv())
            );
            volScalarField dragCoeff(Kd(nodei, nodej));
            dragCoeffMean += dragCoeff;
            volVectorField deltaM
            (
                (phase1_->U(nodei) - phase2_->U(nodej))/XiD
               *(1.0/(dragCoeff*XiD*deltaT + 1.0) - 1.0)
            );
            volScalarField deltaE
            (
                (phase1_->thermo().T() - phase2_->thermo().T())/XiE
               *(exp(-max(Kh(nodei, nodej)*XiE*deltaT, vSmall)) - 1.0)
            );

            phase1_->alphaRhoURef(nodei) += deltaM;
            phase2_->alphaRhoURef(nodei) -= deltaM;

            phase1_->alphaRhoERef() += deltaE;
            phase2_->alphaRhoERef() -= deltaE;

            if (!phase1_->granular())
            {
                phase1_->alphaRhoERef() -=
                    alphaRho2*0.5
                   *(
                        magSqr(phase2_->alphaRhoU(nodej)/alphaRho2)
                      - magSqr(phase2_->U(nodej))
                    );
            }
            else
            {
                phase1_->alphaRhoURef(nodei) -=
                    deltaT*phase1_->gradp()*phase1_->alphas(nodei)
                   /max(phase1_(), phase1_->residualAlpha());
            }
            if (!phase2_->granular())
            {
                phase2_->alphaRhoERef() -=
                    alphaRho1*0.5
                   *(
                        magSqr(phase1_->alphaRhoU(nodei)/alphaRho1)
                      - magSqr(phase1_->U(nodei))
                    );
            }
            else
            {
                phase2_->alphaRhoURef(nodei) -=
                    deltaT*phase2_->gradp()*phase2_->alphas(nodei)
                   /max(phase2_(), phase2_->residualAlpha());
            }
        }
    }

    if
    (
        (phase1_->granular() && phase1_->nNodes() > 1)
     || (phase2_->granular() && phase2_->nNodes() > 1)
    )
    {
        phaseModel* particles;
        phaseModel* gas;
        if (phase1_->granular() && phase1_->nNodes() > 1)
        {
            particles = &phase1_();
            gas = &phase2_();
        }
        else
        {
            particles = &phase2_();
            gas = &phase1_();
        }

        volScalarField ThetaOld = particles->Theta();
        particles->correctThermo();

        volScalarField deltaAlphaRhoPTEp
        (
            3.0/2.0*(*particles)*particles->rho()*(particles->Theta() - ThetaOld)
        );

        particles->alphaRhoERef() -= deltaAlphaRhoPTEp;
    }
    else if (phase1_->granular() || phase2_->granular())
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
        const volScalarField& alphaRhop = particles->alphaRho()();
        volScalarField& Theta = particles->ThetaRef();
        volScalarField ThetaStar
        (
            Theta*exp(-max(2.0*dragCoeffMean*deltaT/alphaRhop, vSmall))
        );
        ThetaStar.max(vSmall);

        volScalarField XiSlip
        (
            particles->productionCoeff()
           /(dragCoeffMean*XiDMean*deltaT + 1.0)
        );

        volScalarField ThetaStarStar
        (
            pow
            (
                XiSlip/alphaRhop*deltaT
              + pow(ThetaStar, 1.5), 2.0/3.0
            )
        );
        gas->alphaRhoERef() -= 1.5*alphaRhop*(ThetaStarStar - Theta);

        Theta =
            ThetaStarStar*9.0*sqr(alphaRhop)
           /sqr
            (
                3.0*alphaRhop
              + deltaT*particles->dissipationCoeff()*sqrt(ThetaStarStar)
            );

        particles->alphaRhoPTERef() = 1.5*alphaRhop*Theta;
        particles->alphaRhoPTERef().correctBoundaryConditions();

        particles->alphaRhoERef() -= 1.5*alphaRhop*(Theta - ThetaStarStar);
    }

    phase1_->alphaRhoURef().correctBoundaryConditions();
    phase2_->alphaRhoURef().correctBoundaryConditions();
    phase1_->alphaRhoERef().correctBoundaryConditions();
    phase2_->alphaRhoERef().correctBoundaryConditions();

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
