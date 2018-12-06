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

Class
    Foam::firstOrderKinetic

Description
    First-order velocity moment advection with kinetic scheme.

SourceFiles
    firstOrderKineticVelocityAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef AUSMPlusKineticVelocityAdvection_H
#define AUSMPlusKineticVelocityAdvection_H

#include "velocityMomentAdvection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityAdvection
{

/*---------------------------------------------------------------------------*\
                    Class AUSMPlusKinetic Declaration
\*---------------------------------------------------------------------------*/

class AUSMPlusKinetic
:
    public velocityMomentAdvection
{
    // Private data

        //- List of quadrature nodes used for moment advection
        const PtrList<volVectorNode>& nodes_;

        //- List of interpolated nodes (neighbour)
        autoPtr<PtrList<surfaceVectorNode>> nodesNei_;

        //- List of interpolated nodes (owner)
        autoPtr<PtrList<surfaceVectorNode>> nodesOwn_;

        PtrList<surfaceScalarField> momentsOwn_;

        PtrList<surfaceScalarField> momentsNei_;

        PtrList<surfaceScalarField> momentsf_;

        //- Reference to density field
        const volScalarField& rho_;

        //- Reference to frictional pressure field
        const volScalarField& pFric_;

        //- Maximum volume fraction
        dimensionedScalar alphaMax_;

        //- Critical volume fraction
        dimensionedScalar alphaMinFriction_;

        //- Cut off mach Number
        dimensionedScalar cutOffMa_;

        //- Residual density
        dimensionedScalar residualRho_;

        //- Minimum velocity
        dimensionedScalar residualU_;

        //- Model Coefficient
        scalar beta_ = 0.125;
        scalar fa_;
        scalar D_;
        scalar xi_;


        //- Flux moment inversion method
//         autoPtr<hyperbolicFieldMomentInversion> momentFieldInverter_;

        // Private functions

            //- Return interpolation scheme name
            word interpScheme(const word& fieldName) const
            {
                return "reconstruct(" + fieldName + ")";
            }

            //- First polynomial
            tmp<surfaceScalarField> M1
            (
                const surfaceScalarField& M,
                const label sign
            );

            //- Second polynomial
            tmp<surfaceScalarField> M2
            (
                const surfaceScalarField& M,
                const label sign
            );

protected:

    // Protected member functions

        // Interpolate weights and abscissae
        void interpolateNodes();

        //- Compute collisions with walls
        void updateWallCollisions();


public:

    //- Runtime type information
    TypeName("AUSMPlusKinetic");


    // Constructors

        //- Construct from velocityMomentSet
        AUSMPlusKinetic
        (
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            const word& support
        );


    //- Destructor
    virtual ~AUSMPlusKinetic();


    // Public member functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const;

        //- Return the maximum Courant number based on velocities
        virtual scalar CoNum() const;

        //- Update moment advection
        virtual void update();

        //- Update moment advection using a specified velocity field
        virtual void update
        (
            const surfaceScalarField& U,
            const bool wallCollisions
        );

        //- Update moment advection using a specified list of velocity fields
        virtual void update
        (
            const mappedPtrList<volVectorField>& Us,
            const bool wallCollisions
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace velocityAdvection
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
