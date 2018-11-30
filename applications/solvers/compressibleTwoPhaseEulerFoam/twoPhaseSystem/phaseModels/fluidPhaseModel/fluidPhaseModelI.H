/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

inline const Foam::volScalarField& Foam::phaseModel::alphaRho() const
{
    return alphaRho_;
}

inline Foam::volScalarField& Foam::phaseModel::alphaRho()
{
    return alphaRho_;
}

inline const Foam::volVectorField& Foam::phaseModel::alphaRhoU() const
{
    return alphaRhoU_;
}

inline Foam::volVectorField& Foam::phaseModel::alphaRhoU()
{
    return alphaRhoU_;
}

inline const Foam::volScalarField& Foam::phaseModel::alphaRhoE() const
{
    return alphaRhoE_;
}

inline Foam::volScalarField& Foam::phaseModel::alphaRhoE()
{
    return alphaRhoE_;
}

inline const Foam::surfaceScalarField& Foam::phaseModel::massFlux() const
{
    return massFlux_;
}

inline Foam::surfaceScalarField& Foam::phaseModel::massFlux()
{
    return massFlux_;
}

inline const Foam::surfaceVectorField& Foam::phaseModel::momentumFlux() const
{
    return momentumFlux_;
}

inline Foam::surfaceVectorField& Foam::phaseModel::momentumFlux()
{
    return momentumFlux_;
}

inline const Foam::surfaceScalarField& Foam::phaseModel::energyFlux() const
{
    return energyFlux_;
}

inline Foam::surfaceScalarField& Foam::phaseModel::energyFlux()
{
    return energyFlux_;
}

inline const Foam::volVectorField& Foam::phaseModel::gradAlpha() const
{
    return gradAlpha_;
}

inline Foam::volVectorField& Foam::phaseModel::gradAlpha()
{
    return gradAlpha_;
}

// ************************************************************************* //
