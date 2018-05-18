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
    Foam::momentGenerationSubModels::alphaAndDiameter

Description
    Moment generation done using volume fraction, density, and bubble diameter.

SourceFiles
    alphaRhoDNodeGenerationModel.C

\*---------------------------------------------------------------------------*/

#ifndef alphaRhoDNodeGenerationModel_H
#define alphaRhoDNodeGenerationModel_H

#include "nodeGenerationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nodeGenerationModels
{

/*---------------------------------------------------------------------------*\
                 Class alphaRhoDNodeGenerationModel Declaration
\*---------------------------------------------------------------------------*/

class alphaRhoDNodeGenerationModel
:
    public nodeGenerationModel
{
public:

        //- Runtime type information
        TypeName("alphaRhoD");


    // Constructors

        //- Construct from components
        alphaRhoDNodeGenerationModel
        (
            const word& nodeGenerationType,
            const labelListList& nodeIndexes,
            mappedList<scalar>& weights,
            mappedList<scalarList>& abscissae,
            const label cmpt
        );


    //- Destructor
    virtual ~alphaRhoDNodeGenerationModel();


    // Public member functions

        //- Update weights and abscissae given a node index and dictionary
        virtual void updateNodes
        (
            const dictionary& dict,
            const label vCmpt
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nodeGenerationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
