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

#include "HLLPhaseFlux.H"
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxFunctions
{
    defineTypeNameAndDebug(HLLPhase, 0);
    addToRunTimeSelectionTable(phaseFluxFunction, HLLPhase, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxFunctions::HLLPhase::HLLPhase
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    phaseFluxFunction(mesh, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxFunctions::HLLPhase::~HLLPhase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseFluxFunctions::HLLPhase::updateFluxes
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

    // Averages
    surfaceScalarField aBar("aBar", 0.5*(aOwn + aNei));
    surfaceScalarField rhoBar("rhoBar", 0.5*(rhoOwn + rhoNei));

    // Fluxes
    surfaceScalarField massFluxOwn(alphaOwn*rhoOwn*UvOwn);
    surfaceScalarField massFluxNei(alphaNei*rhoNei*UvNei);

    surfaceVectorField momentumFluxOwn
    (
        UOwn*massFluxOwn
      + alphaOwn*pOwn*normal
    );
    surfaceVectorField momentumFluxNei
    (
        UNei*massFluxNei
      + alphaNei*pNei*normal
    );

    surfaceScalarField energyFluxOwn(HOwn*massFluxOwn);
    surfaceScalarField energyFluxNei(HNei*massFluxNei);

    surfaceScalarField SOwn("SOwn", min(UvOwn - aOwn, UvNei - aNei));
    surfaceScalarField SNei("SNei", max(UvNei + aNei, UvOwn + aNei));
    surfaceScalarField rDeltaS("rDeltaS", 1.0/(SNei - SOwn));

    pf_ =
        pos0(SOwn)*pOwn
      + pos0(SNei)*neg(SOwn)*(SNei*pOwn - SOwn*pNei)*rDeltaS
      + neg(SNei)*pNei;

    // Compute fluxes
    alphaf_ =
        pos0(SOwn)*alphaOwn
      + pos0(SNei)*neg(SOwn)*(SNei*alphaNei - SOwn*alphaOwn)*rDeltaS
      + neg(SNei)*alphaNei;

    massFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*massFluxOwn
          + pos0(SNei)*neg(SOwn)
           *(
                SNei*massFluxOwn - SOwn*massFluxNei
              + SOwn*SNei*(alphaNei*rhoNei - alphaOwn*rhoOwn)
            )*rDeltaS
          + neg(SNei)*massFluxNei
        );
    Uf_ =
        pos0(SOwn)*UOwn
      + pos0(SNei)*neg(SOwn)
       *(
            SNei*UOwn - SOwn*UNei + SOwn*SNei*normal
        )*rDeltaS
      + neg(SNei)*UNei;
    phi_ = Uf_ & mesh_.Sf();

    momentumFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*momentumFluxOwn
          + pos0(SNei)*neg(SOwn)
           *(
                SNei*momentumFluxOwn - SOwn*momentumFluxNei
              + SOwn*SNei*(alphaNei*rhoNei*UNei - alphaOwn*rhoOwn*UOwn)
            )*rDeltaS
          + neg(SNei)*momentumFluxNei
        );

    energyFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*energyFluxOwn
          + pos0(SNei)*neg(SOwn)
           *(
                SNei*energyFluxOwn - SOwn*energyFluxNei
              + SOwn*SNei*(alphaNei*rhoNei*ENei - alphaOwn*rhoOwn*EOwn)
            )*rDeltaS
          + neg(SNei)*energyFluxNei
        );
}


void Foam::phaseFluxFunctions::HLLPhase::updateFluxes
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
    surfaceVectorField normal(mesh_.Sf()/mesh_.magSf());
    alphaf_ = alphaf;

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

    // Averages
    surfaceScalarField aBar("aBar", 0.5*(aOwn + aNei));
    surfaceScalarField rhoBar("rhoBar", 0.5*(rhoOwn + rhoNei));

    // Fluxes
    surfaceScalarField massFluxOwn(alphaf*rhoOwn*UvOwn);
    surfaceScalarField massFluxNei(alphaf*rhoNei*UvNei);

    surfaceVectorField momentumFluxOwn
    (
        UOwn*massFluxOwn
      + alphaf*pOwn*normal
    );
    surfaceVectorField momentumFluxNei
    (
        UNei*massFluxNei
      + alphaf*pNei*normal
    );

    surfaceScalarField energyFluxOwn(HOwn*massFluxOwn);
    surfaceScalarField energyFluxNei(HNei*massFluxNei);

    surfaceScalarField SOwn("SOwn", min(UvOwn - aOwn, UvNei - aNei));
    surfaceScalarField SNei("SNei", max(UvNei + aNei, UvOwn + aNei));
    surfaceScalarField rDeltaS("rDeltaS", 1.0/(SNei - SOwn));

    pf_ =
        pos0(SOwn)*pOwn
      + pos0(SNei)*neg(SOwn)*(SNei*pOwn - SOwn*pNei)*rDeltaS
      + neg(SNei)*pNei;

    // Compute fluxes
    massFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*massFluxOwn
          + pos0(SNei)*neg(SOwn)
           *(
                SNei*massFluxOwn - SOwn*massFluxNei
              + SOwn*SNei*alphaf*(rhoNei - rhoOwn)
            )*rDeltaS
          + neg(SNei)*massFluxNei
        );
    Uf_ =
        pos0(SOwn)*UOwn
      + pos0(SNei)*neg(SOwn)
       *(
            SNei*UOwn - SOwn*UNei + SOwn*SNei*normal
        )*rDeltaS
      + neg(SNei)*UNei;
    phi_ = Uf_ & mesh_.Sf();

    momentumFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*momentumFluxOwn
          + pos0(SNei)*neg(SOwn)
           *(
                SNei*momentumFluxOwn - SOwn*momentumFluxNei
              + SOwn*SNei*alphaf*(rhoNei*UNei - rhoOwn*UOwn)
            )*rDeltaS
          + neg(SNei)*momentumFluxNei
        );

    energyFlux =
        mesh_.magSf()
       *(
            pos0(SOwn)*energyFluxOwn
          + pos0(SNei)*neg(SOwn)
           *(
                SNei*energyFluxOwn - SOwn*energyFluxNei
              + SOwn*SNei*alphaf*(rhoNei*ENei - rhoOwn*EOwn)
            )*rDeltaS
          + neg(SNei)*energyFluxNei
        );
}
