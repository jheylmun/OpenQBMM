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

Class
    Foam::granularFluxFunction

Description
    Interpolates conservative and primative variables and updates fluxes

SourceFiles
    granularFluxFunction.C
    newGranularFluxFunction.C

\*---------------------------------------------------------------------------*/

#ifndef granularFluxFunction_H
#define granularFluxFunction_H

#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class granularFluxFunction Declaration
\*---------------------------------------------------------------------------*/

class granularFluxFunction
{
protected:
    // Protected Data

        //- Refrence to mesh
        const fvMesh& mesh_;

        const dictionary& dict_;

        //- Owner interpolation
        surfaceScalarField own_;

        //- Neighbour interpolation
        surfaceScalarField nei_;

        //- Residual alpha
        dimensionedScalar residualAlpha_;

        //- Cutoff Mach number
        scalar epsilon_;

        //- Face interpolated volume fraction
        surfaceScalarField alphaf_;

        //- Volumetric flux
        surfaceScalarField phi_;

        //- Interpolated velocity
        surfaceVectorField Uf_;

        //- Face interpolated pressure
        surfaceScalarField pf_;


    //- Protected functions

        //- Return interpolation scheme name
        word interpScheme(const word& fieldName) const
        {
            return "reconstruct(" + fieldName + ")";
        }

public:

    //- Runtime type information
    TypeName("granularFluxFunction");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            granularFluxFunction,
            dictionary,
            (
                const fvMesh& mesh,
                const word& phaseName
            ),
            (mesh, phaseName)
        );


    // Constructor
    granularFluxFunction(const fvMesh& mesh, const word& phaseName);

    // Selector
    static autoPtr<granularFluxFunction>
    New(const fvMesh& mesh, const word& phaseName);

    //- Destructor
    virtual ~granularFluxFunction();


    // Public Functions

        //- Return the initerpolated volume fraction
        const surfaceScalarField& alphaf() const
        {
            return alphaf_;
        }

        //- Return volumetric flux
        const surfaceScalarField& phi() const;

        //- Return the volume flux of volume fraction
        tmp<surfaceScalarField> alphaPhi() const;

        //- Return the volume fraction gradient
        tmp<volVectorField> gradAlpha() const;

        //- Return the velocity gradient
        tmp<volTensorField> gradU() const;

        //- Return the pressure gradient
        tmp<volVectorField> gradp() const;

        virtual void updateFluxes
        (
            surfaceScalarField& massFlux,
            surfaceVectorField& momentumFlux,
            surfaceScalarField& energyFlux,
            surfaceScalarField& PTEFlux,
            const volScalarField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& E,
            const volScalarField& Theta,
            const volScalarField& p,
            const volScalarField& a,
            const volScalarField& pi
        ) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif