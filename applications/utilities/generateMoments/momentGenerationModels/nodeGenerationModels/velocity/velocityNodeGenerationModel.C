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

#include "velocityNodeGenerationModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace nodeGenerationModels
{
    defineTypeNameAndDebug(velocityNodeGenerationModel, 0);

    addToRunTimeSelectionTable
    (
        nodeGenerationModel,
        velocityNodeGenerationModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nodeGenerationModels::velocityNodeGenerationModel::
velocityNodeGenerationModel
(
    const word& nodeGenerationType,
    const labelListList& nodeIndexes,
    mappedList<scalar>& weights,
    mappedList<scalarList>& abscissae,
    const label cmpt
)
:
    nodeGenerationModel
    (
        nodeGenerationType,
        nodeIndexes,
        weights,
        abscissae,
        cmpt
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nodeGenerationModels::velocityNodeGenerationModel::
~velocityNodeGenerationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nodeGenerationModels::velocityNodeGenerationModel::updateNodes
(
    const dictionary& dict,
    const label vCmpt
)
{
    forAll(weights_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        word nodeName = "node" + mappedList<scalar>::listToWord(nodeIndex);
        if(dict.found(nodeName))
        {
            dictionary nodeDict(dict.subDict(nodeName));
            vector U(nodeDict.lookup("U"));

            abscissae_(nodeIndex)[cmpt_] = U.component(vCmpt);
        }
    }
}


// ************************************************************************* //
