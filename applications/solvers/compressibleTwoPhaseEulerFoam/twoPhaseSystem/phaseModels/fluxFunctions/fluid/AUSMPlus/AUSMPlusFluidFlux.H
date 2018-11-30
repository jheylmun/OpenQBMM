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
    Foam::fluidFluxFunctions::AUSMPlusFlux

Description
    Interpolates conservative and primative variables and updates fluxes

SourceFiles
    AUSMPlusFluidFlux.C

\*---------------------------------------------------------------------------*/

#ifndef AUSMPlusFluidFlux_H
#define AUSMPlusFluidFlux_H

#include "fluidFluxFunction.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidFluxFunctions
{

/*---------------------------------------------------------------------------*\
                  Class fluxFunction Declaration
\*---------------------------------------------------------------------------*/

class AUSMPlusFlux
:
    public fluidFluxFunction
{
    // Private Data

        //- Coefficients
        scalar ku_;
        scalar kp_;

        //- Cut off mach Number
        dimensionedScalar cutOffMa_;

        //- Residual density
        dimensionedScalar residualRho_;

        //- Minimum velocity
        dimensionedScalar residualU_;


        // Private functions

            //- First polynomial
            tmp<surfaceScalarField> f
            (
                const surfaceScalarField& M,
                const label sign
            );

            //- Second polynomial
            tmp<surfaceScalarField> beta
            (
                const surfaceScalarField& M,
                const label sign,
                const surfaceScalarField& fa
            );


public:

    //- Runtime type information
    TypeName("AUSMPlus");

    // Costructor
    AUSMPlusFlux(const fvMesh& mesh, const word& phaseName);


    //- Destructor
    virtual ~AUSMPlusFlux();


    // Public Functions

        //- Update fluxes
        virtual void updateFluxes
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
        );

        //- Update phase fluxes with a predefined alpha
        virtual void updateFluxes
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
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidFluxFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif