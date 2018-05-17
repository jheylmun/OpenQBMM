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
    Foam::firstOrderKineticUnivariate

Description
    First-order univariate moment advection with kinetic scheme.

SourceFiles
    firstOrderKineticUnivariateAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef firstOrderKineticUnivariateAdvection_H
#define firstOrderKineticUnivariateAdvection_H

#include "momentAdvection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace momentAdvectionSchemes
{

/*---------------------------------------------------------------------------*\
                    Class firstOrderKinetic Declaration
\*---------------------------------------------------------------------------*/

class firstOrderKineticUnivariate
:
    public momentAdvection
{
    // Private data

        //- List of quadrature nodes used for moment advection
        autoPtr<mappedPtrList<volNode>> nodes_;

        //- List of interpolated nodes (neighbour)
        autoPtr<mappedPtrList<surfaceNode>> nodesNei_;

        //- List of interpolated nodes (owner)
        autoPtr<mappedPtrList<surfaceNode>> nodesOwn_;

        //- List of interpolated moments (neighbour)
        surfaceMomentFieldSet momentsNei_;

        //- List of interpolated moments (owner)
        surfaceMomentFieldSet momentsOwn_;

        //- Flux moment inversion method
        autoPtr<basicFieldMomentInversion> momentFieldInverter_;

        //- Number of nodes used for advection
        label nNodes_;


protected:

        //- Protected member functions
        void interpolateNodes();


public:

    //- Runtime type information
    TypeName("firstOrderKineticUnivariate");


    // Constructors

        //- Construct from univariateMomentSet
        firstOrderKineticUnivariate
        (
            const dictionary& dict,
            const quadratureApproximation& quadrature,
            const surfaceScalarField& phi,
            const word& support
        );


    //- Destructor
    virtual ~firstOrderKineticUnivariate();


    // Public member functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const;

        //- Return the maximum Courant number based on velocities
        virtual scalar CoNum() const;

        //- Update moment advection
        virtual void update
        (
            const bool localPhi = true,
            const bool wallCollisions = false
        );

        //- Update moment advection using a specified list of velocity fields
        virtual void update
        (
            const mappedPtrList<volVectorField>& Us,
            const bool wallCollisions = false
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace momentAdvectionSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
