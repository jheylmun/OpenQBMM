/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2018 Alberto Passalacqua
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

Class
            Foam::momentAdvection

Description
    Abstract class for the advection of a moment set.

SourceFiles
    momentAdvection.C
    momentAdvectionI.H

\*---------------------------------------------------------------------------*/

#ifndef momentAdvection_H
#define momentAdvection_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "mappedPtrList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "quadratureNodes.H"
#include "momentFieldSets.H"
#include "basicFieldMomentInversion.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class momentAdvection Declaration
\*---------------------------------------------------------------------------*/

class momentAdvection
{
    // Private member functions

        //- Disallow default bitwise copy construct
        momentAdvection(const momentAdvection&);

        //- Disallow default bitwise assignment
        void operator=(const momentAdvection&);


protected:

    // Protected data

        //- Name of the distribution associated to the quadrature approximation.
        const word name_;

        //- Reference to the moments to advect
        const volMomentFieldSet& moments_;

        //- Number of moments
        const label nMoments_;

        //- Advection terms for moment transport equations
        mappedPtrList<volScalarField> divMoments_;

        //- Field used to find upwind values of the owner
        surfaceScalarField own_;

        //- Field used to find upwind values of the neighbour
        surfaceScalarField nei_;

        //- Conservative advection flux
        const surfaceScalarField& phi_;

        //- Support of the distribution function
        word support_;

         //- Moment indicies ues to map moments
        labelListList momentOrders_;

        //- Node indicies used to map nodes
        labelListList nodeIndexes_;

        //- Dimensions of the underlying distribution
        label nDimensions_;


public:

    //- Runtime type information
    TypeName("momentAdvection");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            momentAdvection,
            dictionary,
            (
                const dictionary& dict,
                const quadratureApproximation& quadrature,
                const surfaceScalarField& phi,
                const word& support
            ),
            (dict, quadrature, phi, support)
        );


    // Constructors

        //- Construct from univariateMomentSet
        momentAdvection
        (
            const dictionary& dict,
            const quadratureApproximation& quadrature,
            const surfaceScalarField& phi,
            const word& support
        );


    //- Destructor
    virtual ~momentAdvection();


    // Selectors

        static autoPtr<momentAdvection> New
        (
            const dictionary& dict,
            const quadratureApproximation& quadrature,
            const surfaceScalarField& phi,
            const word& support
        );


    // Member Functions

        //- Return the advection term for the moment transport equation
        inline const mappedPtrList<volScalarField>& divMoments() const;

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const = 0;

        //- Return the maximum Courant number based on velocities
        virtual scalar CoNum() const = 0;

        //- Update moment advection
        virtual void update
        (
            const bool localPhi = true,
            const bool wallCollisions = false
        ) = 0;

        //- Update moment advection using a specified list of velocity fields
        virtual void update
        (
            const mappedPtrList<volVectorField>& Us,
            const bool wallCollisions = false
        ) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "momentAdvectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
