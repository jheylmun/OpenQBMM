/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 Alberto Passalacqua
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
    Foam::momentGenerationModels::noneMomentGenerationModel

Description
    Moments are read directly from the dictionary to create the moment fields

SourceFiles
    noneMomentGenerationModel.C

\*---------------------------------------------------------------------------*/

#ifndef noneMomentGenerationModel_H
#define noneMomentGenerationModel_H

#include "momentGenerationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationModels
{

/*---------------------------------------------------------------------------*\
                        Class noneMomentGenerationModel Declaration
\*---------------------------------------------------------------------------*/

class noneMomentGenerationModel
:
    public momentGenerationModel
{
public:

        //- Runtime type information
        TypeName("moments");


    // Constructors

        //- Construct from components
        noneMomentGenerationModel
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const label nNodes
        );


    //- Destructor
    virtual ~noneMomentGenerationModel();


    // Public member functions

        //- Update quadrature given a dictionary
        virtual void updateQuadrature(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace momentGenerationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
