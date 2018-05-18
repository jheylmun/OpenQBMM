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

#include "momentGenerationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(momentGenerationModel, 0);
    defineRunTimeSelectionTable(momentGenerationModel, dictionary);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationModel::updateMoments()
{
    forAll(moments_, mi)
    {
        const labelList& momentOrder = momentOrders_[mi];
        moments_[mi].value() = 0.0;

        forAll(abscissae_, nodei)
        {

            scalar absCmpt = 1.0;
            forAll(abscissae_[0], cmpti)
            {
                absCmpt *=
                    pow
                    (
                        abscissae_[nodei][cmpti],
                        momentOrder[cmpti]
                    );
            }
            moments_[mi].value() += weights_[nodei]*absCmpt;
        }
    }
}

void Foam::momentGenerationModel::reset()
{
    forAll(abscissae_, nodei)
    {
        abscissae_[nodei] = scalarList(momentOrders_[0].size(), 0.0);
        weights_[nodei] = 0.0;
    }
    forAll(moments_, mi)
    {
        moments_[mi].value() = 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationModel::momentGenerationModel
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const label nNodes
)
:
    dict_(dict),
    momentOrders_(momentOrders),
    nodeIndexes_(nodeIndexes),
    weights_(nodeIndexes_.size(), nodeIndexes_, 0.0),
    abscissae_
    (
        nodeIndexes_.size(),
        nodeIndexes_,
        scalarList(momentOrders_[0].size(), 0.0)
    ),
    moments_(momentOrders_.size(), momentOrders_, 0.0)
{
    dimensionSet wDims(dict.lookup("weightDimension"));
    forAll(moments_, mi)
    {
        const labelList& momentOrder = momentOrders_[mi];
        dimensionSet absCmptDims(dimless);
        forAll(momentOrders_[mi], cmpti)
        {
            word absName = "abscissa" + Foam::name(cmpti) + "Dimension";
            dimensionSet absDim
            (
                (
                    dict.found(absName)
                  ? dict.lookup(absName)
                  : dict.lookup("abscissaDimension")
                )
            );
            absCmptDims *= pow(absDim, momentOrder[cmpti]);
        }
        moments_[mi].dimensions().reset(wDims*absCmptDims);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationModel::~momentGenerationModel()
{}


// ************************************************************************* //
