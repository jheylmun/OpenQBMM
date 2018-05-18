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
    Foam::multivariateMomentInversion

Description
    Implementation of the hyperbolic conditional quadrature method of moment
    for velocity. Three nodes are used in each direction. Moment realizability
    is verified analytically.

SourceFiles
    multivariateMomentInversion.C
    multivariateMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateMomentInversion_H
#define multivariateMomentInversion_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "mappedLists.H"
#include "multivariateMomentSet.H"
#include "runTimeSelectionTables.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class multivariateMomentInversion Declaration
\*---------------------------------------------------------------------------*/
class multivariateMomentInversion
{
protected:

    // Protected member data

        //- Number of dimensions in velocity space
        const label nDimensions_;

        //- Number of moments
        const label nMoments_;

        //- List of the number of nodes
        const label nNodes_;

        //- Support of the NDF
        const word support_;

        //- Moment orders
        const labelListList& momentOrders_;

        //- Node indexes
        const labelListList& nodeIndexes_;

        //- Moments
        mappedList<scalar> moments_;

        //- Quadratrure weights
        mappedList<scalar> weights_;

        //- Quadratrure abscissae
        mappedList<scalarList> abscissae_;


    // Protected member functions

        //- Reset inverter
        void reset();

// private:
//
//     // Private member functions
//
//         //- Disallow default bitwise copy construct
//         multivariateMomentInversion(const multivariateMomentInversion&);
//
//         //- Disallow default bitwise assignment
//         void operator=(const multivariateMomentInversion&);


public:

    //- Runtime type information
    TypeName("multivariateMomentInversion");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            multivariateMomentInversion,
            dictionary,
            (
                const dictionary& dict,
                const labelListList& momentOrders,
                const labelListList& nodeIndexes
            ),
            (dict, momentOrders, nodeIndexes)
        );


    // Constructors

        //- Construct from univariateMomentSet
        multivariateMomentInversion
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes
        );


    //- Destructor
    virtual ~multivariateMomentInversion();


    // Selectors

        static autoPtr<multivariateMomentInversion> New
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes
        );


    // Member Functions

        //- Invert moments to find weights and  abscissae
        virtual void invert(const multivariateMomentSet& moments) = 0;

        //- Return number of moments
        inline label nMoments();

        //- Return number of number of nodes in each dimension
        inline label nNodes();

        //- Return mapped quadrature weights for each dimension
        inline const mappedList<scalar>& weights() const;

        //- Return mapped quadrature abscissae for each dimension
        inline const mappedList<scalarList>& abscissae() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multivariateMomentInversionI.H"

#endif


// ************************************************************************* //
