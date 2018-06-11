/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2018-05-22 Jeff Heylmun:    Added the return of one cell
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

Class
    Foam::kineticTheoryModels::radialModels::LunSavage

Description

SourceFiles
    LunSavage.C

\*---------------------------------------------------------------------------*/

#ifndef LunSavage_H
#define LunSavage_H

#include "radialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{

/*---------------------------------------------------------------------------*\
                           Class LunSavage Declaration
\*---------------------------------------------------------------------------*/

class LunSavage
:
    public radialModel
{

public:

    //- Runtime type information
    TypeName("LunSavage");


    // Constructors

        //- Construct from components
        LunSavage(const dictionary& dict);


    //- Destructor
    virtual ~LunSavage();


    // Member Functions

        //- Radial distribution function
        virtual tmp<volScalarField> g0
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax
        ) const;

        //- Derivative of the radial distribution function
        virtual tmp<volScalarField> g0prime
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax
        ) const;

        //- Radial distribution function in celli
        virtual scalar g0
        (
            const scalar& alpha,
            const scalar& alphaMinFriction,
            const scalar& alphaMax
        ) const;

        //- Derivative of the radial distribution function in celli
        virtual scalar g0prime
        (
            const scalar& alpha,
            const scalar& alphaMinFriction,
            const scalar& alphaMax
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radialModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
