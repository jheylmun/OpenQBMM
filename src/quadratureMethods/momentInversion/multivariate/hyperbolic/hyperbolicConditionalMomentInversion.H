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

Class
    Foam::multivariateMomentInversions::hyperbolicConditionalMomentInversion

Description
    Implementation of the hyperbolic conditional quadrature method of moment
    for velocity. Three nodes are used in each direction. Moment realizability
    is verified analytically.

SourceFiles
    hyperbolicConditionalMomentInversion.C
    hyperbolicConditionalMomentInversionI.C

\*---------------------------------------------------------------------------*/

#ifndef hyperbolicConditionalMomentInversion_H
#define hyperbolicConditionalMomentInversion_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "vector2D.H"
#include "scalarMatrices.H"
#include "dictionary.H"
#include "univariateMomentSet.H"
#include "multivariateMomentInversion.H"
#include "hyperbolicMomentInversion.H"


namespace Foam
{
namespace multivariateMomentInversions
{

/*---------------------------------------------------------------------------*\
            Class hyperbolicConditionalMomentInversion Declaration
\*---------------------------------------------------------------------------*/
class hyperbolicConditionalMomentInversion
:
    public multivariateMomentInversion
{
    // Private member data

        //- Univariate moment inversion method
        autoPtr<hyperbolicMomentInversion> univariateInverter_;

        //- Minimum value of the variance to compute eta
        scalar etaMin_;

        //- Maximum normalized skewness
        scalar qMax_;

        //- Small negative value for realizability
        scalar smallNegRealizability_;

        //- Minimum value of variance
        scalar varMin_;

        //- Minimum value of correlation
        scalar minCorrelation_;


    // Private member functions

        // Calculate q from q and eta
        scalar calcQ(scalar q, scalar eta);

        // Check realizability of univariate moments
        void realizabilityUnivariateMoments
        (
            scalar& c2,
            scalar& c3,
            scalar& c4
        );

        //- One-dimensional HyQMOM inversion
        void invert1D
        (
            const multivariateMomentSet& moments,
            scalarList& weights1D,
            scalarList& abscissae1D
        );

        //- Two-dimensional HyCQMOM inversion
        void invert2D
        (
            const multivariateMomentSet& moments,
            mappedList<scalar>& weights2D,
            mappedList<vector2D>& abscissae2D
        );

        //- Three-dimensional HyCQMOM inversion
        void invert3D
        (
            const multivariateMomentSet& moments
        );


public:

    TypeName("HyCQMOM");

    // Constructors

        // Construct from nMoments, momentMap, nNodes, and support type
        hyperbolicConditionalMomentInversion
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes
        );


    //- Destructor
    virtual ~hyperbolicConditionalMomentInversion();

    // Static data

        //- Return the number of nodes needed
        static label nHyperbolicNodes(const label nDims);

        //- Return the node indexes
        static labelListList hyperbolicNodeIndexes(const label nDims);

        //- Return the number of moments needed
        static label nHyperbolicMoments(const label nDims);

        //- Return the momentOrders
        static labelListList hyperbolicMomentOrders(const label nDims);


    // Member Functions

        //- Invert moments to find weights and  abscissae
        virtual void invert(const multivariateMomentSet& moments);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace multivariateMomentInversions
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //
