/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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
    Foam::PDFTransportModel

Description
    Abstract class for the run-time selection of the PDF transport model.

SourceFiles
    PDFTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef PDFTransportModel_H
#define PDFTransportModel_H

#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "quadratureApproximation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PDFTransportModel Declaration
\*---------------------------------------------------------------------------*/

class PDFTransportModel
{
    // Private data

        //- Name of the PDFTransportModel
        const word name_;

        //- Absolute tolerance for realizable ODE solver
        scalar ATol_;

        //- Relative tolerance for realizable ODE solver
        scalar RTol_;

        //- Time-step change factor
        scalar fac_;

        //- Minimum time-step change factor
        scalar facMin_;

        //- Maximum time-step change factor
        scalar facMax_;

        //- Minimum local step in realizable ODE solver
        scalar minLocalDt_;
    // Private Member Functions

        //- Disallow default bitwise copy construct
        PDFTransportModel(const PDFTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const PDFTransportModel&);


protected:

    // Protected data

        //- Constant reference to the mesh
        const fvMesh& mesh_;

        //- Quadrature approximation used to solve PDF transport
        quadratureApproximation quadrature_;

        //- Are source Terms solved with the ode solver
        Switch solveODESource_;


    // Protected member functions

        //- Calculate explicit source terms with realizable ODE solver
        virtual void explicitMomentSource(quadratureApproximation&);

        //- Calculate source term for moment equation in the specified cell
        virtual scalar cellMomentSource
        (
            const label momentOrder,
            const label celli
        ) = 0;


public:

    // Constructors

        //- Construct from components
        PDFTransportModel
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh,
            const word& support
        );

    //- Destructor
    virtual ~PDFTransportModel();


    // Member Functions

        //- Returns a const reference to the quadrature approximation
        inline const quadratureApproximation& quadrature() const;

        //- Returns a non-const reference to the quadrature approximation
        inline quadratureApproximation& quadrature();

        //- Solve PDF transport equation
        virtual void solve() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PDFTransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
