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

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::momentGenerationModel> Foam::momentGenerationModel::New
(
    const dictionary& dict,
    const label nNodes,
    const bool extended,
    const bool Radau
)
{
    word momentGenerationModelType(dict.lookup("type"));

    Info<< "Selecting momentGenerationModel "
        << momentGenerationModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(momentGenerationModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "momentGenerationModel::New" << nl
            << "(" << nl
            << "    const dictionary&," << nl
            << "    const label" << nl
            << ") : " << endl
            << "    unknown momentGenerationModel type "
            << momentGenerationModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid momentGenerationModel types are :" << endl;

        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<momentGenerationModel>
    (
        cstrIter()
        (
            dict,
            nNodes,
            extended,
            Radau
        )
    );
}


// ************************************************************************* //
