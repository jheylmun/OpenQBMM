/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Alberto Passalacqua
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

#include "BoltzmannCollision.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace collisionKernels
{
    defineTypeNameAndDebug(BoltzmannCollision, 0);

    addToRunTimeSelectionTable
    (
        collisionKernel,
        BoltzmannCollision,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::updateCells1D(const label celli)
{
}

void Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::updateCells2D(const label celli)
{
}

void Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::updateCells3D(const label celli)
{
}


void Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::updateFields1D()
{
}

void Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::updateFields2D()
{
}

void Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::updateFields3D()
{
    const volVectorField& U1 = quadrature_.nodes()[node1].primaryAbscissa();
    const volVectorField& U2 = quadrature_.nodes()[node2].primaryAbscissa();
    volScalarField v1(U1.component(0));
    volScalarField v2(U1.component(1));
    volScalarField v3(U1.component(2));
    volVectorField g = U1 - U2;
    volScalarField g1(g.component(0));
    volScalarField g1Sqr(sqr(g1));
    volScalarField g2(g.component(1));
    volScalarField g2Sqr(sqr(g2));
    volScalarField g3(g.component(2));
    volScalarField g3Sqr(sqr(g3));
    volScalarField gMag = mag(g);
    volScalarField gMagSqr = sqr(gMag);
    scalar omega = (1 + e)/2.0;
    scalar omegaSqr = sqr(omega);

    Isf_(0,0,0) = 0;
    Isf_(1,0,0) = -(omega/2.0)*g1;
    Isf_(0,1,0) = -(omega/2.0)*g2;
    Isf_(0,0,1) = -(omega/2.0)*g3;
    Isf_(2,0,0) = (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g1Sqr - omega*g1*v1;
    Isf_(1,1,0) = (omegaSqr/4.0)*g1*g2 + (omega/4.0)*g1*g2 - 0.5*g1*v2
    Isf_(1,0,1) = (omegaSqr/4.0)*g1*g3 + (omega/4.0)*g1*g3 - 0.5*g1*v3
    Isf_(0,2,0) = (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g2Sqr - omega*g2*v2;
    Isf_(0,1,1) = (omegaSqr/4.0)*g2*g3 + (omega/4.0)*g2*g3 - 0.5*g2*v3
    Isf_(0,0,2) = (omegaSqr/12.0)*gMagSqr + (omegaSqr/4.0)*g3Sqr - omega*g3*v3;
    Isf_(3,0,0) =
        (omega*omegaSqr/8.0)*(qSqr + g1Sqr)
      + (omegaSqr/4.0)*(gSqr + 3.0*g1Sqr)*v1
      - (1.5*omega)*g1*sqr(v1);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision::BoltzmannCollision
(
    const dictionary& dict,
    const fvMesh& mesh,
    const velocityQuadratureApproximation& quadrature,
    const bool ode
)
:
    collisionKernel(dict, mesh, quadrature, ode),
    tauCollisional_(dict.lookup("tau")),
    e_(dict.lookupType<scalar>("e")),
    Meqf_(quadrature.moments().size(), momentOrders_),
    Meq_(quadrature.moments().size(), momentOrders_),
    Isf_(quadrature.moments().size(), momentOrders_),
    Is_(quadrature.moments().size(), momentOrders_)
{
    if (!ode)
    {
//         Meqf_.setSize(quadrature.moments().size());

        forAll(Meqf_, mi)
        {
            const labelList& momentOrder = momentOrders_[mi];
            Meqf_.set
            (
                momentOrder,
                new volScalarField
                (
                    IOobject
                    (
                        "Meq" + mappedList<scalar>::listToWord(momentOrder),
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar
                    (
                        "zero",
                        quadrature.moments()[mi].dimensions(),
                        0.0
                    )
                )
            );

            Isf.set
            (
                momentOrder,
                new volScalarField
                (
                    IOobject
                    (
                        "I" + mappedList<scalar>::listToWord(momentOrder),
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar
                    (
                        "zero",
                        quadrature.moments()[mi].dimensions(),
                        0.0
                    )
                )
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::~BoltzmannCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::explicitCollisionSource(const label mi, const label celli) const
{
    return (quadrature_.moments()[mi][celli] - Meq_[mi])/tauCollisional_.value();
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::populationBalanceSubModels::collisionKernels::BoltzmannCollision
::implicitCollisionSource(const volVectorMoment& m) const
{
    return
    (
        Meqf_(m.cmptOrders())/tauCollisional_
      - fvm::Sp(1/tauCollisional_, m)
    );
}
// ************************************************************************* //
