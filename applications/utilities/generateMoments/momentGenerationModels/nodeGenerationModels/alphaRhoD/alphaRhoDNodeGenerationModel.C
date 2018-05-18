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

#include "alphaRhoDNodeGenerationModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace nodeGenerationModels
{
    defineTypeNameAndDebug(alphaRhoDNodeGenerationModel, 0);

    addToRunTimeSelectionTable
    (
        nodeGenerationModel,
        alphaRhoDNodeGenerationModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nodeGenerationModels::alphaRhoDNodeGenerationModel::
alphaRhoDNodeGenerationModel
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

Foam::nodeGenerationModels::alphaRhoDNodeGenerationModel::
~alphaRhoDNodeGenerationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nodeGenerationModels::alphaRhoDNodeGenerationModel::updateNodes
(
    const dictionary& dict,
    const label vCmpt
)
{
    forAll(nodeIndexes_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        word nodeName = "node" + mappedList<scalar>::listToWord(nodeIndex);

        if(dict.found(nodeName))
        {
            dictionary nodeDict(dict.subDict(nodeName));
            scalar dia(readScalar(nodeDict.lookup("d")));
            scalar alpha(readScalar(nodeDict.lookup("alpha")));
            scalar rho(readScalar(nodeDict.lookup("rho")));

            abscissae_(nodeIndex)[cmpt_]
                = (4.0/3.0)*Foam::constant::mathematical::pi*rho*pow3(dia/2.0);

            weights_(nodeIndex) *= rho*alpha/abscissae_(nodeIndex)[cmpt_];
        }
    }
}


// ************************************************************************* //
