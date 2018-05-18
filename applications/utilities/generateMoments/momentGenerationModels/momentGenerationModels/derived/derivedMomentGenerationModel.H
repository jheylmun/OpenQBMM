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
    Foam::momentGenerationModel

Description
    Base class for runtime selection for moment generation. One object is used
    for all node. Moments are constructed outside of class.

SourceFiles
    derivedMomentGenerationModel.C

\*---------------------------------------------------------------------------*/

#ifndef derivedMomentGenerationModel_H
#define derivedMomentGenerationModel_H

#include "dictionary.H"
#include "mappedLists.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "momentGenerationModel.H"
#include "nodeGenerationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationModels
{

class derivedMomentGenerationModel
:
    public momentGenerationModel
{
    // Private data

        PtrList<nodeGenerationModel> nodeGenerators_;


public:

    //- Runtime information
    TypeName("derived");


    // Constructor

        //- Construct from dictionary and number of nodes
        derivedMomentGenerationModel
        (
            const dictionary& dict,
            const labelListList& momentsOrders,
            const labelListList& nodeIndexes,
            const label nNodes
        );


    //- Destructor
    virtual ~derivedMomentGenerationModel();


    // Public Functions

        // Member Functions

            //- Update weights and abscissae given a dictionary
            virtual void updateQuadrature(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace momentGenerationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
