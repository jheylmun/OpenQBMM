/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
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

#include "sizeHyCQMOM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multivariateMomentInversions
{
    defineTypeNameAndDebug(sizeHyCQMOM, 0);

    addToRunTimeSelectionTable
    (
        multivariateMomentInversion,
        sizeHyCQMOM,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::sizeHyCQMOM::sizeHyCQMOM
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes
)
:
    multivariateMomentInversion(dict, momentOrders, nodeIndexes),
    nGeometricDimensions_(nDimensions_ - 1),
    nSizeMoments_(calcNSizeMoments(momentOrders)),
    velocityMomentOrders_
    (
        hyperbolicConditionalMomentInversion::hyperbolicMomentOrders
        (
            nGeometricDimensions_
        )
    ),
    velocityNodeIndexes_
    (
        hyperbolicConditionalMomentInversion::hyperbolicNodeIndexes
        (
            nGeometricDimensions_
        )
    ),
    supports_({"RPlus", "R", "R", "R"}),
    sizeInverter_
    (
        univariateMomentInversion::New(dict.subDict("basicQuadrature"))
    ),
    velocityInverter_
    (
        dict,
        velocityMomentOrders_,
        velocityNodeIndexes_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::sizeHyCQMOM::~sizeHyCQMOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label
Foam::multivariateMomentInversions::sizeHyCQMOM::calcNSizeMoments
(
    const labelListList& momentOrders
)
{
    label maxOrder = 0;
    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];
        if (momentOrder[0] > maxOrder)
        {
            maxOrder = momentOrder[0];
        }
    }
    return maxOrder + 1;
}


void
Foam::multivariateMomentInversions::sizeHyCQMOM::invert
(
    const multivariateMomentSet& moments
)
{
    reset();

    //- Invert size moments and build VR matrix
    univariateMomentSet sizeMoments(nSizeMoments_, supports_[0], 0.0);
    labelList order(nDimensions_, 0);
    forAll(sizeMoments, mi)
    {
        order[0] = mi;
        sizeMoments[mi] = moments(order);
    }
    sizeInverter_->invert(sizeMoments);
    const scalarList& sizeWeights(sizeInverter_->weights());
    const scalarList& sizeAbscissae(sizeInverter_->abscissae());

    forAll(nodeIndexes_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        label sizeNode = nodeIndex[0];
        if (sizeNode <= sizeInverter_->nNodes())
        {
            weights_(nodeIndex) = sizeWeights[sizeNode - 1];
            abscissae_(nodeIndex)[0] = sizeAbscissae[sizeNode - 1];
        }
    }

    label nNonZeroNodes = 0;
    boolList nonZeroNodes(sizeInverter_->nNodes(), false);
    forAll(sizeWeights, nodei)
    {
        // Check if bubble moments are large enough.
        //  If yes make matricies 1 component larger,
        //  if no the rest of the nodes are assumed to
        //  be too small as well.
        //  This is done to avoid a divide by 0 error,
        //  and to reduce unneeded computation time
        if
        (
            sizeWeights[nodei] > 0
         && sizeAbscissae[nodei] > 0
        )
        {
            nonZeroNodes[nodei] = true;
            nNonZeroNodes++;
        }
    }

    if (nNonZeroNodes)
    {
        scalarDiagonalMatrix x(nNonZeroNodes, 0.0);
        scalarSquareMatrix R(nNonZeroNodes, 0.0);
        scalarSquareMatrix invR(nNonZeroNodes, 0.0);

        label i = 0;
        forAll(sizeWeights, nodei)
        {
            if (nonZeroNodes[nodei])
            {
                x[i] = sizeAbscissae[nodei];
                invR[i][i] = 1.0/sizeWeights[nodei];
                i++;
            }
        }
        Vandermonde V(x);
        scalarSquareMatrix invVR = invR*V.inv();

        // Compute conditional velocity moments and invert
        PtrList<mappedList<scalar>> conditionalMoments(sizeInverter_->nNodes());
        forAll(conditionalMoments, sNodei)
        {
            conditionalMoments.set
            (
                sNodei,
                new mappedList<scalar>
                (
                    velocityMomentOrders_.size(),
                    velocityMomentOrders_,
                    0.0
                )
            );
        }

        forAll(velocityMomentOrders_, mi)
        {
            const labelList& velocityMomentOrder = velocityMomentOrders_[mi];
            labelList pureMomentOrder(nDimensions_, 0);
            for (label dimi = 1; dimi < nDimensions_; dimi++)
            {
                pureMomentOrder[dimi] = velocityMomentOrder[dimi - 1];
            }

            scalarRectangularMatrix M(nNonZeroNodes, 1, 0);
            for (label sNodei = 0; sNodei < nNonZeroNodes; sNodei++)
            {
                pureMomentOrder[0] = sNodei;
                M(sNodei, 0) = moments(pureMomentOrder);
            }
            scalarRectangularMatrix nu = invVR*M;

            forAll(conditionalMoments, sNodei)
            {
                conditionalMoments[sNodei](velocityMomentOrder) = nu(sNodei, 0);
            }
        }

        forAll(conditionalMoments, sNodei)
        {
            multivariateMomentSet momentsToInvert
            (
                velocityMomentOrders_.size(),
                velocityMomentOrders_,
                "R"
            );
            forAll(momentsToInvert, mi)
            {
                momentsToInvert(velocityMomentOrders_[mi]) =
                    conditionalMoments[sNodei](velocityMomentOrders_[mi]);
            }
            velocityInverter_.invert(momentsToInvert);

            forAll(velocityNodeIndexes_, nodei)
            {
                const labelList& velocityNodeIndex = velocityNodeIndexes_[nodei];
                labelList nodeIndex(nDimensions_, 0);
                nodeIndex[0] = sNodei + 1;
                for (label dimi = 1; dimi < nDimensions_; dimi++)
                {
                    nodeIndex[dimi] = velocityNodeIndex[dimi - 1];
                }

                weights_(nodeIndex) *=
                    velocityInverter_.weights()(velocityNodeIndex);
                for (label dimi = 0; dimi < nGeometricDimensions_; dimi++)
                {
                    abscissae_(nodeIndex)[dimi + 1] =
                        velocityInverter_.abscissae()
                        (
                            velocityNodeIndex
                        )[dimi];
                }
            }
        }
    }
}


// ************************************************************************* //
