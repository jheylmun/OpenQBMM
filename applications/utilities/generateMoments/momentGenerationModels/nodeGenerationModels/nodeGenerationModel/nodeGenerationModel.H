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
    Foam::nodeGenerationModel

Description
    Base class for runtime selection for moment generation. One object is used
    for all node. Moments are constructed outside of class.

SourceFiles
    nodeGenerationModel.C
    newNodeGenerationModel.C

\*---------------------------------------------------------------------------*/

#ifndef nodeGenerationModel_H
#define nodeGenerationModel_H

#include "dictionary.H"
#include "mappedLists.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nodeGenerationModel
{

    // Private member functions

        //- Disallow default bitwise copy construct
        nodeGenerationModel(const nodeGenerationModel&);

        //- Disallow default bitwise assignment
        void operator=(const nodeGenerationModel&);


protected:

    // Protected data

        //- Node indexes
        const labelListList& nodeIndexes_;

        //- List of weights
        mappedList<scalar>& weights_;

        //- List of abscissae
        mappedList<scalarList>& abscissae_;

        //- Index of abscissa component
        label cmpt_;


public:

    //- Runtime information
    TypeName("nodeGenerationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        nodeGenerationModel,
        dictionary,
        (
            const word& nodeGenerationType,
            const labelListList& nodeIndexes,
            mappedList<scalar>& weights,
            mappedList<scalarList>& abscissae,
            const label cmpt
        ),
        (nodeGenerationType, nodeIndexes, weights, abscissae, cmpt)
    );


    // Constructor

        //- Construct from dictionary and number of nodes
        nodeGenerationModel
        (
            const word& nodeGenerationType,
            const labelListList& nodeIndexes,
            mappedList<scalar>& weights,
            mappedList<scalarList>& abscissae,
            const label cmpt
        );


    // Selector

        static autoPtr<nodeGenerationModel> New
        (
            const word& nodeGenerationType,
            const labelListList& nodeIndexes,
            mappedList<scalar>& weights,
            mappedList<scalarList>& abscissae,
            const label cmpt
        );


    //- Destructor
    virtual ~nodeGenerationModel();


    // Public Functions

        //- Update weights and abscissae given a node index and dictionary
        virtual void updateNodes
        (
            const dictionary& dict,
            const label vCmpt
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
