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

#include "derivedMomentGenerationModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationModels
{
    defineTypeNameAndDebug(derivedMomentGenerationModel, 0);

    addToRunTimeSelectionTable
    (
        momentGenerationModel,
        derivedMomentGenerationModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationModels::derivedMomentGenerationModel::
derivedMomentGenerationModel
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const label nNodes
)
:
    momentGenerationModel(dict, momentOrders, nodeIndexes, nNodes),
    nodeGenerators_(momentOrders[0].size())
{
    wordList abscissaeTypes = dict.lookup("abscissaeTypes");
    forAll(momentOrders_[0], cmpt)
    {
        nodeGenerators_.set
        (
            cmpt,
            nodeGenerationModel::New
            (
                abscissaeTypes[cmpt],
                nodeIndexes_,
                weights_,
                abscissae_,
                cmpt
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationModels::derivedMomentGenerationModel::
~derivedMomentGenerationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::momentGenerationModels::derivedMomentGenerationModel::updateQuadrature
(
    const dictionary& dict
)
{
    reset();
    forAll(nodeIndexes_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        word nodeName = "node" + mappedList<scalar>::listToWord(nodeIndex);
        if(dict.found(nodeName))
        {
            dictionary nodeDict(dict.subDict(nodeName));
            weights_(nodeIndex) =
                nodeDict.lookupOrDefault<scalar>("weight", 1.0);
        }
    }

    label vCmpt = 0;
    forAll(abscissae_[0], cmpt)
    {
        Info<<dict.parent()<<endl;
        word absName = "abscissa" + Foam::name(cmpt) + "Dimension";
        if (dict.parent().found(absName))
        {
            dimensionSet absDims = dict.parent().lookup(absName);
            word prevAbsName =
                "abscissa" + Foam::name(cmpt - 1) + "Dimension";
            if
            (
                dict.parent().found(prevAbsName)
             && (
                    absDims
                 == dict.parent().lookupType<dimensionSet>(prevAbsName)
                )
             && cmpt > 0
            )
            {
                vCmpt++;
            }
            else
            {
                vCmpt = 0;
            }
        }
        else if (cmpt > 0)
        {
            vCmpt++;
        }
        else
        {
            vCmpt = 0;
        }
        Info<<vCmpt<<endl;
        nodeGenerators_[cmpt].updateNodes(dict, vCmpt);
    }

    updateMoments();
}


// ************************************************************************* //
