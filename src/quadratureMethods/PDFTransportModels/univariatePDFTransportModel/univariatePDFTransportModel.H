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
    Foam::univariatePDFTransportModel

Description
    Solve a univariate PDF transport equation with the extended quadrature
    method of moments.

SourceFiles
    univariatePDFTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef univariatePDFTransportModel_H
#define univariatePDFTransportModel_H

#include "PDFTransportModel.H"
#include "momentAdvection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{

/*---------------------------------------------------------------------------*\
                Class univariatePDFTransportModel Declaration
\*---------------------------------------------------------------------------*/

class univariatePDFTransportModel
:
    public PDFTransportModel
{
protected:

    // Protected data

        //- Moment advection
        autoPtr<momentAdvection> momentAdvection_;

    // Protected member functions

        //- Diffusion term
        virtual tmp<fvScalarMatrix> momentDiffusion
        (
            const volMoment& moment
        ) = 0;

        //- Calculate implicit source terms
        virtual tmp<fvScalarMatrix> implicitMomentSource
        (
            const volMoment& moment
        ) = 0;


public:

    // Constructors

        //- Construct from components
        univariatePDFTransportModel
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh,
            const surfaceScalarField& phi,
            const word& support
        );


    //- Destructor
    virtual ~univariatePDFTransportModel();


    // Member Functions

        //- Access

            //- Return the maximum Courant number ensuring moment realizability
            inline scalar realizableCo() const;

        //-Edit

            //- Solve univariate PDF transport equation
            virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PDFTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariatePDFTransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
