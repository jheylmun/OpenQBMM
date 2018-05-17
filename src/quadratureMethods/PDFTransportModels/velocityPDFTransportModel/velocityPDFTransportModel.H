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
    Foam::velocityPDFTransportModel

Description
    Solve a velocity PDF transport equation with the hyperbolic conditional
    quadrature method of moments

SourceFiles
    velocityPDFTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef velocityPDFTransportModel_H
#define velocityPDFTransportModel_H

#include "PDFTransportModel.H"
#include "momentAdvection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{

/*---------------------------------------------------------------------------*\
                Class velocityPDFTransportModel Declaration
\*---------------------------------------------------------------------------*/

class velocityPDFTransportModel
:
    public PDFTransportModel
{
protected:

    // Protected data

        //- Moment advection
        autoPtr<momentAdvection> momentAdvection_;

    // Protected member functions

        //- Return if collisions are used
        virtual bool collision() const = 0;

        //- Update implicit collison source terms
        virtual void updateImplicitCollisionSource() = 0;

        //- Update explicit collison source terms
        virtual void updateExplicitCollisionSource(const label celli) = 0;

        //- Calculate implicit source terms
        virtual tmp<fvScalarMatrix> implicitCollisionSource
        (
            const volMoment& moment
        ) = 0;

        virtual scalar explicitCollisionSource
        (
            const label momenti,
            const label celli
        ) = 0;


public:

    // Constructors

        //- Construct from components
        velocityPDFTransportModel
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh,
            const surfaceScalarField& phi,
            const word& support
        );


    //- Destructor
    virtual ~velocityPDFTransportModel();


    // Member Functions

        //- Access

            //- Return the maximum Courant number ensuring moment realizability
            virtual inline scalar realizableCo() const;

            //- Return the maximum Courant number based on velocities
            virtual inline scalar CoNum() const;


        //-Edit

            //- Solve velocity PDF transport equation
            virtual void solve();

            //- Solve mean tranport
            virtual void meanTransport(const bool wallCollisions);

            //- Solve relative transport
            virtual void relativeTransport
            (
                const mappedPtrList<volVectorField>& Vs,
                const bool wallCollisions
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PDFTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "velocityPDFTransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
