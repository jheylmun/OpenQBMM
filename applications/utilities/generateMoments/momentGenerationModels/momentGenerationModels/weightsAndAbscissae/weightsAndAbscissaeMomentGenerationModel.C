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

\*---------------------------------------------------------------------------*/

#include "weightsAndAbscissaeMomentGenerationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationModels
{
    defineTypeNameAndDebug(weightsAndAbscissaeMomentGenerationModel, 0);

    addToRunTimeSelectionTable
    (
        momentGenerationModel,
        weightsAndAbscissaeMomentGenerationModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationModels::weightsAndAbscissaeMomentGenerationModel
::weightsAndAbscissaeMomentGenerationModel
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const label nNodes
)
:
    momentGenerationModel(dict, momentOrders, nodeIndexes, nNodes)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationModels::weightsAndAbscissaeMomentGenerationModel
::~weightsAndAbscissaeMomentGenerationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationModels::weightsAndAbscissaeMomentGenerationModel::
updateQuadrature
(
    const dictionary& dict
)
{
    reset();
    forAll(nodeIndexes_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        word nodeName = "node" + mappedList<scalar>::listToWord(nodeIndex);
        Info<<nodeName<<endl;
        if(dict.found(nodeName))
        {
            dictionary nodeDict(dict.subDict(nodeName));
            abscissae_(nodeIndex) = nodeDict.lookupType<scalarList>("abscissa");
            weights_(nodeIndex) = nodeDict.lookupType<scalar>("weight");
        }
    }

    updateMoments();
}


// ************************************************************************* //
