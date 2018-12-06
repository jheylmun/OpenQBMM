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

inline const Foam::word& Foam::phaseModel::name() const
{
    return name_;
}

inline const Foam::twoPhaseSystem& Foam::phaseModel::fluid() const
{
    return fluid_;
}

inline const Foam::dimensionedScalar& Foam::phaseModel::residualAlpha() const
{
    return residualAlpha_;
}

inline const Foam::dimensionedScalar& Foam::phaseModel::residualRho() const
{
    return residualRho_;
}

inline Foam::dimensionedScalar Foam::phaseModel::residualAlphaRho() const
{
    return residualAlpha_*residualRho_;
}

inline const Foam::rhoThermo& Foam::phaseModel::thermo() const
{
    return thermoPtr_();
}

inline Foam::rhoThermo& Foam::phaseModel::thermoRef()
{
    return thermoPtr_();
}

inline const Foam::volScalarField& Foam::phaseModel::rho() const
{
    return rho_;
}

inline Foam::volScalarField& Foam::phaseModel::rhoRef()
{
    return rho_;
}

inline const Foam::volScalarField& Foam::phaseModel::E() const
{
    return E_;
}

inline Foam::volScalarField& Foam::phaseModel::ERef()
{
    return E_;
}

inline Foam::tmp<Foam::volScalarField> Foam::phaseModel::nu() const
{
    return thermoPtr_->nu();
}

inline Foam::tmp<Foam::scalarField> Foam::phaseModel::nu
(
    const label patchi
) const
{
    return thermoPtr_->nu(patchi);
}

inline Foam::tmp<Foam::volScalarField> Foam::phaseModel::mu() const
{
    return thermoPtr_->mu();
}

inline Foam::tmp<Foam::scalarField> Foam::phaseModel::mu
(
    const label patchi
) const
{
    return thermoPtr_->mu(patchi);
}

inline Foam::tmp<Foam::scalarField> Foam::phaseModel::kappa
(
    const label patchi
) const
{
    return thermoPtr_->kappa(patchi);
}

inline Foam::tmp<Foam::volScalarField> Foam::phaseModel::kappa() const
{
    return thermoPtr_->kappa();
}

inline Foam::tmp<Foam::volScalarField> Foam::phaseModel::alphahe() const
{
    return thermoPtr_->alphahe();
}


inline Foam::tmp<Foam::scalarField> Foam::phaseModel::alphahe
(
    const label patchi
) const
{
    return thermoPtr_->alphahe(patchi);
}


inline Foam::tmp<Foam::volScalarField> Foam::phaseModel::kappaEff
(
    const volScalarField& alphat
) const
{
    return thermoPtr_->kappaEff(alphat);
}

inline Foam::tmp<Foam::scalarField> Foam::phaseModel::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermoPtr_->kappaEff(alphat, patchi);
}

inline Foam::tmp<Foam::volScalarField> Foam::phaseModel::alpha() const
{
    return thermoPtr_->alpha();
}

inline Foam::tmp<Foam::scalarField> Foam::phaseModel::alpha
(
    const label patchi
) const
{
    return thermoPtr_->alpha(patchi);
}

inline Foam::tmp<Foam::volScalarField> Foam::phaseModel::alphaEff
(
    const volScalarField& alphat
) const
{
    return thermoPtr_->alphaEff(alphat);
}

inline Foam::tmp<Foam::scalarField> Foam::phaseModel::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermoPtr_->alphaEff(alphat, patchi);
}

inline Foam::tmp<Foam::volScalarField> Foam::phaseModel::Cp() const
{
    return thermoPtr_->Cp();
}

inline Foam::tmp<Foam::volScalarField> Foam::phaseModel::Cv() const
{
    return thermoPtr_->Cv();
}

inline Foam::tmp<Foam::volScalarField> Foam::phaseModel::gamma() const
{
    return thermoPtr_->gamma();
}


// ************************************************************************* //
