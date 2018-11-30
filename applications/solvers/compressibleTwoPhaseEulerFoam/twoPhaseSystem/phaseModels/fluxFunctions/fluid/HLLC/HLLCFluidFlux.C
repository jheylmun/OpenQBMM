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

#include "HLLCFluidFlux.H"
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidFluxFunctions
{
    defineTypeNameAndDebug(HLLCFlux, 0);
    addToRunTimeSelectionTable(fluidFluxFunction, HLLCFlux, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidFluxFunctions::HLLCFlux::HLLCFlux
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    fluidFluxFunction(mesh, phaseName),
    residualU_("small", dimVelocity, epsilon_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidFluxFunctions::HLLCFlux::~HLLCFlux()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluidFluxFunctions::HLLCFlux::updateFluxes
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

    surfaceScalarField UvOwn(UOwn & normal);
    surfaceScalarField UvNei(UNei & normal);

    surfaceScalarField aOwn
    (
        fvc::interpolate(a, own_, interpScheme(a.name()))
    );
    surfaceScalarField aNei
    (
        fvc::interpolate(a, nei_, interpScheme(a.name()))
    );

    // Averages
    surfaceScalarField aTilde
    (
        "aTilde",
        (sqrt(rhoOwn)*aOwn + sqrt(rhoNei)*aNei)
       /(sqrt(rhoOwn) + sqrt(rhoNei))
    );
    surfaceVectorField uTilde
    (
        "uTilde",
        (sqrt(rhoOwn)*UOwn + sqrt(rhoNei)*UNei)
       /(sqrt(rhoOwn) + sqrt(rhoNei))
    );

    // Compute wave speeds
    surfaceScalarField SOwn
    (
        "SOwn", min(UvOwn - aOwn, (uTilde & normal) - aTilde)
    );
    surfaceScalarField SNei
    (
        "SNei", max(UvNei + aNei, (uTilde & normal) + aTilde)
    );

    //- Star quantities
    surfaceVectorField SStar
    (
        "SStar",
        (
            (pNei - pOwn)*normal
          + rhoOwn*UOwn*(SOwn - UvOwn)
          - rhoNei*UNei*(SNei - UvNei)
        )/(rhoOwn*(SOwn - UvOwn) - rhoNei*(SNei - UvNei))
    );
    surfaceScalarField sstar(SStar & normal);
    surfaceScalarField pStar
    (
        "pOwnNei",
        (
            pOwn + rhoOwn*(SOwn - UvOwn)*(sstar - UvOwn)
          + pNei + rhoNei*(SNei - UvNei)*(sstar - UvNei)
        )*0.5
    );
    surfaceScalarField rhoOwnStar
    (
        "rhoOwnStar",
        rhoOwn*(SOwn - UvOwn)/(SOwn - sstar)
    );
    surfaceScalarField rhoNeiStar
    (
        "rhoNeiStar",
        rhoNei*(SNei - UvNei)/(SNei - sstar)
    );
    surfaceScalarField EOwnStar
    (
        "EOwnStar",
        EOwn + (pStar*sstar - pOwn*UvOwn)/(rhoOwn*(SOwn - UvOwn))
    );
    surfaceScalarField ENeiStar
    (
        "ENeiStar",
        ENei + (pStar*sstar - pNei*UvNei)/(rhoNei*(SNei - UvNei))
    );

    // Reimann primitive quantities
    alphaf_ = alphaOwn;

    // Reimann primitive quantities
    surfaceScalarField rhoR("rhoR", rhoOwn);
    Uf_ = UOwn;
    surfaceScalarField ER("ER", EOwn);
    pf_ = pOwn;

    forAll(rhoR, facei)
    {
        if (SOwn[facei] < 0 && sstar[facei] >= 0)
        {
            rhoR[facei] = rhoOwnStar[facei];
            Uf_[facei] = SStar[facei];
            pf_[facei] = pStar[facei];
            ER[facei] = EOwnStar[facei];
        }
        else if (sstar[facei] < 0 && SNei[facei] >= 0)
        {
            alphaf_[facei] = alphaNei[facei];
            rhoR[facei] = rhoNeiStar[facei];
            Uf_[facei] = SStar[facei];
            pf_[facei] = pStar[facei];
            ER[facei] = ENeiStar[facei];
        }
        else if (SNei[facei] < 0)
        {
            alphaf_[facei] = alphaNei[facei];
            rhoR[facei] = rhoNei[facei];
            Uf_[facei] = UNei[facei];
            pf_[facei] = pNei[facei];
            ER[facei] = ENei[facei];
        }
    }
    phi_ = Uf_ & mesh_.Sf();

    // Set total fluxes
    massFlux = alphaf_*rhoR*phi();
    momentumFlux = massFlux*Uf_ + alphaf_*pf_*mesh_.Sf();
    energyFlux = alphaf_*phi()*(rhoR*ER + pf_);
}


void Foam::fluidFluxFunctions::HLLCFlux::updateFluxes
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

    surfaceScalarField UvOwn(UOwn & normal);
    surfaceScalarField UvNei(UNei & normal);

    surfaceScalarField aOwn(fvc::interpolate(a, own_, interpScheme(a.name())));
    surfaceScalarField aNei(fvc::interpolate(a, nei_, interpScheme(a.name())));

    // Averages
    surfaceScalarField aTilde
    (
        "aTilde",
        (sqrt(rhoOwn)*aOwn + sqrt(rhoNei)*aNei)
       /(sqrt(rhoOwn) + sqrt(rhoNei))
    );
    surfaceVectorField uTilde
    (
        "uTilde",
        (sqrt(rhoOwn)*UOwn + sqrt(rhoNei)*UNei)
       /(sqrt(rhoOwn) + sqrt(rhoNei))
    );

    // Compute wave speeds
    surfaceScalarField SOwn
    (
        "SOwn", min(UvOwn - aOwn, (uTilde & normal) - aTilde)
    );
    surfaceScalarField SNei
    (
        "SNei", max(UvNei + aNei, (uTilde & normal) + aTilde)
    );

    //- Star quantities
    surfaceVectorField SStar
    (
        "SStar",
        (
            (pNei - pOwn)*normal
          + rhoOwn*UOwn*(SOwn - UvOwn)
          - rhoNei*UNei*(SNei - UvNei)
        )/(rhoOwn*(SOwn - UvOwn) - rhoNei*(SNei - UvNei))
    );
    surfaceScalarField sstar(SStar & normal);
    surfaceScalarField pStar
    (
        "pOwnNei",
        (
            pOwn + rhoOwn*(SOwn - UvOwn)*(sstar - UvOwn)
          + pNei + rhoNei*(SNei - UvNei)*(sstar - UvNei)
        )*0.5
    );
    surfaceScalarField rhoOwnStar
    (
        "rhoOwnStar",
        rhoOwn*(SOwn - UvOwn)/(SOwn - sstar)
    );
    surfaceScalarField rhoNeiStar
    (
        "rhoNeiStar",
        rhoNei*(SNei - UvNei)/(SNei - sstar)
    );
    surfaceScalarField EOwnStar
    (
        "EOwnStar",
        EOwn + (pStar*sstar - pOwn*UvOwn)/(rhoOwn*(SOwn - UvOwn))
    );
    surfaceScalarField ENeiStar
    (
        "ENeiStar",
        ENei + (pStar*sstar - pNei*UvNei)/(rhoNei*(SNei - UvNei))
    );

    // Reimann primitive quantities
//     surfaceScalarField posSOwn(pos0(SOwn));
//     surfaceScalarField SOwnSStar(neg(SOwn)*pos0(sstar));
//     surfaceScalarField SStarSNei(neg(sstar)*pos0(SNei));
//     surfaceScalarField negSNei(neg(SNei));
//
//     surfaceScalarField rhoR
//     (
//         "rhoR",
//         posSOwn*rhoOwn
//       + SOwnSStar*rhoOwnStar
//       + SStarSNei*rhoNeiStar
//       + negSNei*rhoNei
//     );
//     Uf_ =
//         posSOwn*UOwn
//       + min(neg(SOwn) + pos0(SNei), 1.0)*SStar
//       + negSNei*UNei;
//     phi_ = Uf_ & mesh_.Sf();
//
//     surfaceScalarField ER
//     (
//         "ER",
//         posSOwn*EOwn
//       + SOwnSStar*EOwnStar
//       + SStarSNei*ENeiStar
//       + negSNei*ENei
//     );
//
//     pf_ =
//         posSOwn*pOwn
//       + min(neg(SOwn) + pos0(SNei), 1.0)*pStar
//       + negSNei*pNei;
    surfaceScalarField rhoR
    (
        IOobject
        (
            "rhoR",
            rho.time().timeName(),
            rho.mesh()
        ),
        rho.mesh(),
        dimensionedScalar("0", dimDensity, 0.0)
    );
    Uf_ = dimensionedVector("0", dimVelocity, Zero);
    pf_ = dimensionedScalar("0", dimPressure, 0.0);
    surfaceScalarField ER
    (
        IOobject
        (
            "ER",
            E.time().timeName(),
            E.mesh()
        ),
        E.mesh(),
        dimensionedScalar("0", E.dimensions(), 0.0)
    );
    forAll(rhoR, facei)
    {
        if (SOwn[facei] >= 0)
        {
            rhoR[facei] = rhoOwn[facei];
            Uf_[facei] = UOwn[facei];
            pf_[facei] = pOwn[facei];
            ER[facei] = EOwn[facei];
        }
        else if (SNei[facei] <= 0)
        {
            rhoR[facei] = rhoNei[facei];
            Uf_[facei] = UNei[facei];
            pf_[facei] = pNei[facei];
            ER[facei] = ENei[facei];
        }
        else if (SOwn[facei] <= 0 && 0 <= sstar[facei])
        {
            rhoR[facei] = rhoOwnStar[facei];
            Uf_[facei] = SStar[facei];
            pf_[facei] = pStar[facei];
            ER[facei] = EOwnStar[facei];
        }
        else// if (sstar[facei] <= 0 && 0 <= SNei[facei])
        {
            rhoR[facei] = rhoNeiStar[facei];
            Uf_[facei] = SStar[facei];
            pf_[facei] = pStar[facei];
            ER[facei] = ENeiStar[facei];
        }
    }
    forAll(rhoR.boundaryField(), patchi)
    {
        forAll(rhoR.boundaryField()[patchi], facei)
        {
            if (SOwn[facei] >= 0)
            {
                rhoR.boundaryFieldRef()[patchi][facei] =
                    rhoOwn.boundaryField()[patchi][facei];
                Uf_.boundaryFieldRef()[patchi][facei] =
                    UOwn.boundaryField()[patchi][facei];
                pf_.boundaryFieldRef()[patchi][facei] =
                    pOwn.boundaryField()[patchi][facei];
                ER.boundaryFieldRef()[patchi][facei] =
                    EOwn.boundaryField()[patchi][facei];
            }
            else if (SNei[facei] <= 0)
            {
                rhoR.boundaryFieldRef()[patchi][facei] =
                    rhoNei.boundaryField()[patchi][facei];
                Uf_.boundaryFieldRef()[patchi][facei] =
                    UNei.boundaryField()[patchi][facei];
                pf_.boundaryFieldRef()[patchi][facei] =
                    pNei.boundaryField()[patchi][facei];
                ER.boundaryFieldRef()[patchi][facei] =
                    ENei.boundaryField()[patchi][facei];
            }
            else if (SOwn[facei] <= 0 && 0 <= sstar[facei])
            {
                rhoR.boundaryFieldRef()[patchi][facei] =
                    rhoOwnStar.boundaryField()[patchi][facei];
                Uf_.boundaryFieldRef()[patchi][facei] =
                    SStar.boundaryField()[patchi][facei];
                pf_.boundaryFieldRef()[patchi][facei] =
                    pStar.boundaryField()[patchi][facei];
                ER.boundaryFieldRef()[patchi][facei] =
                    EOwnStar.boundaryField()[patchi][facei];
            }
            else if (sstar[facei] <= 0 && 0 <= SNei[facei])
            {
                rhoR.boundaryFieldRef()[patchi][facei] =
                    rhoNeiStar.boundaryField()[patchi][facei];
                Uf_.boundaryFieldRef()[patchi][facei] =
                    SStar.boundaryField()[patchi][facei];
                pf_.boundaryFieldRef()[patchi][facei] =
                    pStar.boundaryField()[patchi][facei];
                ER.boundaryFieldRef()[patchi][facei] =
                    ENeiStar.boundaryField()[patchi][facei];
            }
        }
    }
    phi_ = Uf_ & mesh_.Sf();


    // Set total fluxes

    massFlux = alphaf_*rhoR*phi();
    momentumFlux = massFlux*Uf_ + alphaf_*pf_*mesh_.Sf();
    energyFlux = alphaf_*phi()*(rhoR*ER + pf_);
}