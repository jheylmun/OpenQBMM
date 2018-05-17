/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
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
Application
    Test-ExtendedMomentInversion.C
Description
    Test the extendedMomentInversion class and its subclasses.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "scalarMatrices.H"
#include "mappedList.H"
#include "hyperbolicConditionalMomentInversion.H"
#include "Random.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "createFields.H"

    Random rand(clock::getTime());

    mappedList<scalar> w(nNodes, nodeIndexes, 0.0);
    mappedList<vector> u(nNodes, nodeIndexes, Zero);

    forAll(nodeIndexes, nodei)
    {
        const labelList& nodeIndex = nodeIndexes[nodei];

        w(nodeIndex) = rand.scalar01();
        if (nDims == 1)
        {
            u(nodeIndex) =
                vector(rand.scalar01() - 0.5, 0.0, 0.0)*10.0;
        }
        else if (nDims == 2)
        {
            u(nodeIndex) =
                vector
                (
                    rand.scalar01() - 0.5,
                    rand.scalar01() - 0.5,
                    0.0
                )*10;
        }
        else
        {
            u(nodeIndex) =
                vector
                (
                    rand.scalar01() - 0.5,
                    rand.scalar01() - 0.5,
                    rand.scalar01() - 0.5
                )*10;
        }
    }

    multivariateMomentSet moments(nMoments, momentOrders, "R");
    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];
        moments(momentOrder) = 0.0;

        forAll(nodeIndexes, nodei)
        {
            const labelList& nodeIndex = nodeIndexes[nodei];

            scalar cmpt = w(nodeIndex);
            forAll(nodeIndex, dimi)
            {
                cmpt *= pow(u(nodeIndex)[dimi], momentOrder[dimi]);
            }
            moments(momentOrder) += cmpt;
        }
    }
//     moments(0,0) = 5.1441;
//     moments(1,0) = -2.3417;
//     moments(0,1) = 4.0371;
//     moments(2,0) = 35.8806;
//     moments(1,1) = -25.7800;
//     moments(0,2) = 43.6027;
//     moments(3,0) = 11.0056;
//     moments(0,3) = 111.2670;
//     moments(4,0) = 368.2310;
//     moments(0,4) = 764.4060;

    hyperbolicConditionalMomentInversion momentInverter
    (
        quadratureProperties, nDims
    );

    Info<< "\nInverting moments" << endl;

    momentInverter.invert(moments);

    Info<< "\nReconstructed moments:" << endl;

    const mappedList<scalar>& weights = momentInverter.weights();
    const mappedList<vector>& abscissae = momentInverter.abscissae();

    mappedList<scalar> newMoments(nMoments, momentOrders);
    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];
        newMoments(momentOrder) = 0.0;

        forAll(nodeIndexes, nodei)
        {
            const labelList& nodeIndex = nodeIndexes[nodei];

            scalar cmpt = weights(nodeIndex);
            forAll(momentOrder, dimi)
            {
                cmpt *= pow(abscissae(nodeIndex)[dimi], momentOrder[dimi]);
            }
            newMoments(momentOrder) += cmpt;
        }

        Info<< "moment.";
        forAll(momentOrder, dimi)
        {
            Info<< momentOrder[dimi];
        }
        Info<< ":\told: " << moments(momentOrder)
            << "\tnew: " << newMoments(momentOrder)
            << ",\terror: "
            << (mag(moments(momentOrder) - newMoments(momentOrder))/moments(momentOrder))<< endl;
    }
    Info << nl << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //