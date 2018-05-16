/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class weightType, class abscissaType, class sigmaType, class vType>
const Foam::word&
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::name() const
{
    return name_;
}

template <class weightType, class abscissaType, class sigmaType, class vType>
inline bool Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
extended() const
{
    return extended_;
}

template <class weightType, class abscissaType, class sigmaType, class vType>
Foam::label
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
nSecondaryNodes() const
{
    return nSecondaryNodes_;
}

template <class weightType, class abscissaType, class sigmaType, class vType>
const weightType&
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
primaryWeight() const
{
    return weight_;
}

template <class weightType, class abscissaType, class sigmaType, class vType>
weightType&
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
primaryWeight()
{
    return weight_;
}

template <class weightType, class abscissaType, class sigmaType, class vType>
const abscissaType&
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
primaryAbscissa(const label cmpt) const
{
    return abscissa_[cmpt];
}

template <class weightType, class abscissaType, class sigmaType, class vType>
abscissaType&
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
primaryAbscissa(const label cmpt)
{
    return abscissa_[cmpt];
}

template <class weightType, class abscissaType, class sigmaType, class vType>
const Foam::PtrList<weightType>&
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
secondaryWeights() const
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access secondary weights of a quadrature node" << nl
            << "    not of extended type. "
            << abort(FatalError);
    }

    return secondaryWeights_;
}

template <class weightType, class abscissaType, class sigmaType, class vType>
Foam::PtrList<weightType>&
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
secondaryWeights()
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access secondary weights of a quadrature node" << nl
            << "    not of extended type. "
            << abort(FatalError);
    }
    return secondaryWeights_;
}

template <class weightType, class abscissaType, class sigmaType, class vType>
const Foam::PtrList<abscissaType>&
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>
::secondaryAbscissae(const label cmpt) const
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access secondary abscissae of a quadrature" << nl
            << "   node not of extended type. "
            << abort(FatalError);
    }

    return secondaryAbscissae_[cmpt];
}

template <class weightType, class abscissaType, class sigmaType, class vType>
Foam::PtrList<abscissaType>&
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>
::secondaryAbscissae(const label cmpt)
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access secondary abscissae of a quadrature" << nl
            << "   node not of extended type. "
            << abort(FatalError);
    }

    return secondaryAbscissae_[cmpt];
}

template <class weightType, class abscissaType, class sigmaType, class vType>
const sigmaType&
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>
::sigma() const
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access sigma field of a quadrature node not" << nl
            << "   of extended type. "
            << abort(FatalError);
    }

    return sigma_();
}

template <class weightType, class abscissaType, class sigmaType, class vType>
sigmaType&
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::sigma()
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access sigma field of a quadrature node not" << nl
            << "   of extended type. "
            << abort(FatalError);
    }

    return sigma_();
}


// ************************************************************************* //
