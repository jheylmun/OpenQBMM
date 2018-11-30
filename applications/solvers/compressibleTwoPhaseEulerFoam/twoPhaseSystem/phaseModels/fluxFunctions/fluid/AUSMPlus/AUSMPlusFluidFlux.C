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

#include "AUSMPlusFluidFlux.H"
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidFluxFunctions
{
    defineTypeNameAndDebug(AUSMPlusFlux, 0);
    addToRunTimeSelectionTable(fluidFluxFunction, AUSMPlusFlux, dictionary);
}
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::fluidFluxFunctions::AUSMPlusFlux::f
(
    const surfaceScalarField& M,
    const label sign
)
{
    surfaceScalarField magM(mag(M));

    return
        0.5*(M + sign*mag(M))*pos0(magM)
      + (sign*0.25*sqr(M + sign) + sign*0.125*sqr(sqr(M) - 1))*neg(magM);
}

Foam::tmp<Foam::surfaceScalarField>
Foam::fluidFluxFunctions::AUSMPlusFlux::beta
(
    const surfaceScalarField& M,
    const label sign,
    const surfaceScalarField& fa
)
{
    surfaceScalarField A(3.0/16.0*(5.0*sqr(fa) - 4.0));
    surfaceScalarField magM(mag(M));

    return
        0.5*(1.0 + Foam::sign(sign*M))*pos(magM)
      + (
            0.25*(2.0 - sign*M)*sqr(M + sign*1.0)
          + sign*A*M*sqr(sqr(M) - 1)
        )*neg(magM);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidFluxFunctions::AUSMPlusFlux::AUSMPlusFlux
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    fluidFluxFunction(mesh, phaseName),
    ku_(dict_.lookupOrDefault("ku", 1.0)),
    kp_(dict_.lookupOrDefault("kp", 1.0)),
    cutOffMa_("small", dimless, epsilon_),
    residualRho_("small", dimDensity, epsilon_),
    residualU_("small", dimVelocity, epsilon_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidFluxFunctions::AUSMPlusFlux::~AUSMPlusFlux()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluidFluxFunctions::AUSMPlusFlux::updateFluxes
(
    surfaceScalarField& massFlux,
    surfaceVectorField& momentumFlux,
    surfaceScalarField& energyFlux,
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& E,
    const volScalarField& p,
    const volScalarField& a,
    const volVectorField& Ui,
    const volScalarField& pi
)
{
    surfaceVectorField normal(mesh_.Sf()/mesh_.magSf());
    surfaceScalarField alphaOwn
    (
        fvc::interpolate(alpha, own_, interpScheme(alpha.name()))
    );
    surfaceScalarField alphaNei
    (
        fvc::interpolate(alpha, nei_, interpScheme(alpha.name()))
    );

    surfaceScalarField rhoOwn
    (
        fvc::interpolate(rho, own_, interpScheme(rho.name()))
    );
    surfaceScalarField rhoNei
    (
        fvc::interpolate(rho, nei_, interpScheme(rho.name()))
    );

    surfaceVectorField UOwn(fvc::interpolate(U, own_, interpScheme(U.name())));
    surfaceVectorField UNei(fvc::interpolate(U, nei_, interpScheme(U.name())));

    surfaceScalarField EOwn(fvc::interpolate(E, own_, interpScheme(E.name())));
    surfaceScalarField ENei(fvc::interpolate(E, nei_, interpScheme(E.name())));

    surfaceScalarField pOwn(fvc::interpolate(p, own_, interpScheme(p.name())));
    surfaceScalarField pNei(fvc::interpolate(p, nei_, interpScheme(p.name())));

    surfaceScalarField HOwn(EOwn + pOwn/rhoOwn);
    surfaceScalarField HNei(ENei + pNei/rhoNei);

    surfaceScalarField aOwn(fvc::interpolate(a, own_, interpScheme(a.name())));
    surfaceScalarField aNei(fvc::interpolate(a, nei_, interpScheme(a.name())));

    surfaceScalarField UvOwn(UOwn & normal);
    surfaceScalarField UvNei(UNei & normal);

    surfaceScalarField aStar(sqrt(aOwn*aNei));

    // Compute slpit Mach numbers
    surfaceScalarField MaOwn("MaOwn", UvOwn/aStar);
    surfaceScalarField MaNei("MaNei", UvNei/aStar);
    surfaceScalarField magMaOwn(mag(MaOwn));
    surfaceScalarField magMaNei(mag(MaNei));

    surfaceScalarField MaStar(f(MaOwn, 1.0) + f(MaNei, -1.0));
    surfaceScalarField M0
    (
        sqrt(min(1.0, max(sqr((MaOwn + MaNei)*0.5), cutOffMa_)))
    );
    surfaceScalarField fa(M0*(2.0 - M0));

    surfaceScalarField BetaOwn(beta(MaOwn, 1.0, fa));
    surfaceScalarField BetaNei(beta(MaNei, -1.0, fa));

    surfaceScalarField deltaMa
    (
        f(MaOwn, 1.0) - pos0(MaOwn) - f(MaNei, -1.0) + neg(MaNei)
    );
    surfaceScalarField Du
    (
        -ku_*BetaOwn*BetaNei
       *0.5*(alphaOwn*rhoOwn + alphaNei*rhoNei)
       *fa*aStar*(UvNei - UvOwn)
    );
    surfaceScalarField Dp
    (
        -kp_/fa*deltaMa*max(1.0 - sqr(0.5*(MaOwn - MaNei)), 0.0)
       *(alphaNei*pNei - alphaOwn*pOwn)/aStar
    );

    surfaceScalarField mDot
    (
        "mDot",
        0.5*aStar
       *(
            alphaOwn*rhoOwn*max(MaStar, 0.0)
          + alphaNei*rhoNei*min(MaStar, 0.0)
        ) + Dp
    );
    surfaceScalarField alphaP
    (
        BetaOwn*alphaOwn*pOwn + BetaNei*alphaNei*pNei + Du
    );

    alphaf_ = pos0(mDot)*alphaOwn + neg(mDot)*alphaNei;
    pf_ = alphaP/max((alphaOwn + alphaNei)*0.5, residualAlpha_);

    massFlux = mesh_.magSf()*mDot;
    Uf_ = pos0(mDot)*UOwn + neg(mDot)*UNei;
    phi_ = Uf_ & mesh_.Sf();

    momentumFlux =
        mesh_.magSf()*0.5
       *(
            mDot*(UOwn + UNei)
          + mag(mDot)*(UOwn - UNei)
        )
      + alphaP*mesh_.Sf();

    energyFlux =
        mesh_.magSf()*0.5
       *(
            mDot*(HOwn + HNei)
          + mag(mDot)*(HOwn - HNei)
        );
}


void Foam::fluidFluxFunctions::AUSMPlusFlux::updateFluxes
(
    surfaceScalarField& massFlux,
    surfaceVectorField& momentumFlux,
    surfaceScalarField& energyFlux,
    const surfaceScalarField& alphaf,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& E,
    const volScalarField& p,
    const volScalarField& a,
    const volVectorField& Ui,
    const volScalarField& pi
)
{
    alphaf_ = alphaf;
    surfaceVectorField normal(mesh_.Sf()/mesh_.magSf());

    surfaceScalarField rhoOwn
    (
        fvc::interpolate(rho, own_, interpScheme(rho.name()))
    );
    surfaceScalarField rhoNei
    (
        fvc::interpolate(rho, nei_, interpScheme(rho.name()))
    );

    surfaceVectorField UOwn(fvc::interpolate(U, own_, interpScheme(U.name())));
    surfaceVectorField UNei(fvc::interpolate(U, nei_, interpScheme(U.name())));

    surfaceScalarField EOwn(fvc::interpolate(E, own_, interpScheme(E.name())));
    surfaceScalarField ENei(fvc::interpolate(E, nei_, interpScheme(E.name())));

    surfaceScalarField pOwn(fvc::interpolate(p, own_, interpScheme(p.name())));
    surfaceScalarField pNei(fvc::interpolate(p, nei_, interpScheme(p.name())));

    surfaceScalarField HOwn(EOwn + pOwn/rhoOwn);
    surfaceScalarField HNei(ENei + pNei/rhoNei);

    surfaceScalarField aOwn(fvc::interpolate(a, own_, interpScheme(a.name())));
    surfaceScalarField aNei(fvc::interpolate(a, nei_, interpScheme(a.name())));

    surfaceScalarField UvOwn(UOwn & normal);
    surfaceScalarField UvNei(UNei & normal);

    surfaceScalarField aStar(sqrt(aOwn*aNei));

    // Compute slpit Mach numbers
    surfaceScalarField MaOwn("MaOwn", UvOwn/aStar);
    surfaceScalarField MaNei("MaNei", UvNei/aStar);
    surfaceScalarField magMaOwn(mag(MaOwn));
    surfaceScalarField magMaNei(mag(MaNei));

    surfaceScalarField MaStar(f(MaOwn, 1.0) + f(MaNei, -1.0));
    surfaceScalarField M0
    (
        sqrt(min(1.0, max(sqr((MaOwn + MaNei)*0.5), cutOffMa_)))
    );
    surfaceScalarField fa(M0*(2.0 - M0));

    surfaceScalarField BetaOwn(beta(MaOwn, 1.0, fa));
    surfaceScalarField BetaNei(beta(MaNei, -1.0, fa));

    surfaceScalarField deltaMa
    (
        f(MaOwn, 1.0) - pos0(MaOwn) - f(MaNei, -1.0) + neg(MaNei)
    );
    surfaceScalarField Du
    (
        -ku_*BetaOwn*BetaNei
       *0.5*alphaf*(rhoOwn + rhoNei)
       *fa*aStar*(UvNei - UvOwn)
    );
    surfaceScalarField Dp
    (
        -kp_/fa*deltaMa*max(1.0 - sqr(0.5*(MaOwn - MaNei)), 0.0)
       *alphaf*(pNei - pOwn)/aStar
    );

    surfaceScalarField mDot
    (
        "mDot",
        0.5*alphaf*aStar*(rhoOwn + rhoNei)*max(MaStar, 0.0) + Dp
    );
    surfaceScalarField alphaP
    (
        alphaf*(BetaOwn*pOwn + BetaNei*pNei) + Du
    );

    pf_ = alphaP/max(alphaf, residualAlpha_);

    massFlux = mesh_.magSf()*mDot;
    Uf_ = pos0(mDot)*UOwn + neg(mDot)*UNei;
    phi_ = Uf_ & mesh_.Sf();

    momentumFlux =
        mesh_.magSf()*0.5
       *(
            mDot*(UOwn + UNei)
          + mag(mDot)*(UOwn - UNei)
        )
      + alphaP*mesh_.Sf();

    energyFlux =
        mesh_.magSf()*0.5
       *(
            mDot*(HOwn + HNei)
          + mag(mDot)*(HOwn - HNei)
        );
}
