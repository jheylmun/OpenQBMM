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

#include "quadratureNode.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template <class weightType, class abscissaType, class sigmaType, class vType>
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
quadratureNode
(
    const word& name,
    const word& distributionName,
    const fvMesh& mesh,
    const dimensionSet& weightDimensions,
    const PtrList<dimensionSet>& abscissaDimensions,
    const bool extended,
    const label nSecondaryNodes
)
:
    name_(IOobject::groupName(name, distributionName)),
    weight_
    (
        IOobject
        (
            IOobject::groupName("weight", name_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<typename weightType::value_type>
        (
            "zeroWeight",
            weightDimensions,
            pTraits<typename weightType::value_type>::zero
        )
    ),
    abscissa_(abscissaDimensions.size()),
    secondaryWeights_(),
    secondaryAbscissae_(),
    sigma_(),
    nSecondaryNodes_(nSecondaryNodes),
    extended_(extended),
    velocityIndexes_(3, -1)
{
    label dimi = 0;
    forAll(abscissa_, cmpt)
    {
        if (abscissaDimensions[cmpt] == dimVelocity)
        {
            velocityIndexes_[dimi] = cmpt;
            dimi++;
        }

        abscissa_.set
        (
            cmpt,
            new abscissaType
            (
                IOobject
                (
                    IOobject::groupName("abscissa" + Foam::name(cmpt), name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<typename abscissaType::value_type>
                (
                    "zeroAbscissa",
                    abscissaDimensions[cmpt],
                    pTraits<typename abscissaType::value_type>::zero
                )
            )
        );
    }

    if (extended_)
    {
        // Allocating secondary abscissae
        secondaryAbscissae_.setSize(abscissa_.size());
        forAll (secondaryAbscissae_, cmpt)
        {
            secondaryAbscissae_.set
            (
                cmpt,
                new PtrList<abscissaType>(nSecondaryNodes_)
            );

            forAll(secondaryAbscissae_[cmpt], nodei)
            {
                secondaryAbscissae_[cmpt].set
                (
                    nodei,
                    new abscissaType
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                IOobject::groupName
                                (
                                    "secondaryAbscissae_" + Foam::name(cmpt),
                                    Foam::name(nodei)
                                ),
                                name_
                            ),
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimensioned<typename abscissaType::value_type>
                        (
                            "zeroAbscissa",
                            abscissaDimensions[cmpt],
                            pTraits<typename abscissaType::value_type>::zero
                        )
                    )
                );
            }
        }

        // Allocating secondary weights
        secondaryWeights_.setSize(nSecondaryNodes_);
        forAll(secondaryWeights_, nodei)
        {
            secondaryWeights_.set
            (
                nodei,
                new weightType
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "secondaryWeight." + Foam::name(nodei),
                            name_
                        ),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensioned<typename weightType::value_type>
                    (
                        "zeroWeight",
                        dimless,
                        pTraits<typename weightType::value_type>::zero
                    )
                )
            );
        }

        // Allocating sigma
        sigma_ = autoPtr<sigmaType>
        (
            new sigmaType
            (
                IOobject
                (
                    IOobject::groupName("sigma", name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<typename sigmaType::value_type>
                (
                    "zeroSigma",
                    dimless,
                    pTraits<typename sigmaType::value_type>::zero
                )
            )
        );
    }
}


template <class weightType, class abscissaType, class sigmaType, class vType>
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
quadratureNode
(
    const word& name,
    const word& distributionName,
    const fvMesh& mesh,
    const dimensionSet& weightDimensions,
    const PtrList<dimensionSet>& abscissaDimensions,
    const wordList& boundaryTypes,
    const bool extended,
    const label nSecondaryNodes
)
:
    name_(IOobject::groupName(name, distributionName)),
    weight_
    (
        IOobject
        (
            IOobject::groupName("weight", name_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<typename weightType::value_type>
        (
            "zeroWeight",
            weightDimensions,
            pTraits<typename weightType::value_type>::zero
        )
    ),
    abscissa_(abscissaDimensions.size()),
    secondaryWeights_(),
    secondaryAbscissae_(),
    sigma_(),
    nSecondaryNodes_(nSecondaryNodes),
    extended_(extended),
    velocityIndexes_(3, -1)
{
    label dimi = 0;
    forAll(abscissa_, cmpt)
    {
        if (abscissaDimensions[cmpt] == dimVelocity)
        {
            velocityIndexes_[dimi] = cmpt;
            dimi++;
        }

        abscissa_.set
        (
            cmpt,
            new abscissaType
            (
                IOobject
                (
                    IOobject::groupName("abscissa" + Foam::name(cmpt), name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<typename abscissaType::value_type>
                (
                    "zeroAbscissa",
                    abscissaDimensions[cmpt],
                    pTraits<typename abscissaType::value_type>::zero
                )
            )
        );
    }

    if (extended_)
    {
        // Allocating secondary abscissae
        secondaryAbscissae_.setSize(abscissa_.size());
        forAll (secondaryAbscissae_, cmpt)
        {
            secondaryAbscissae_.set
            (
                cmpt,
                new PtrList<abscissaType>(nSecondaryNodes_)
            );

            forAll(secondaryAbscissae_[cmpt], nodei)
            {
                secondaryAbscissae_[cmpt].set
                (
                    nodei,
                    new abscissaType
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                IOobject::groupName
                                (
                                    "secondaryAbscissae_" + Foam::name(cmpt),
                                    Foam::name(nodei)
                                ),
                                name_
                            ),
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimensioned<typename abscissaType::value_type>
                        (
                            "zeroAbscissa",
                            abscissaDimensions[cmpt],
                            pTraits<typename abscissaType::value_type>::zero
                        )
                    )
                );
            }
        }

        // Allocating secondary weights
        secondaryWeights_.setSize(nSecondaryNodes_);
        forAll(secondaryWeights_, nodei)
        {
            secondaryWeights_.set
            (
                nodei,
                new weightType
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "secondaryWeight." + Foam::name(nodei),
                            name_
                        ),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensioned<typename weightType::value_type>
                    (
                        "zeroWeight",
                        dimless,
                        pTraits<typename weightType::value_type>::zero
                    ),
                    boundaryTypes
                )
            );
        }

        sigma_ = autoPtr<sigmaType>
        (
            new sigmaType
            (
                IOobject
                (
                    IOobject::groupName("sigma", name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<typename sigmaType::value_type>
                (
                    "zeroSigma",
                    dimless,
                    pTraits<typename sigmaType::value_type>::zero
                ),
                boundaryTypes
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class weightType, class abscissaType, class sigmaType, class vType>
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
~quadratureNode()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class weightType, class abscissaType, class sigmaType, class vType>
Foam::autoPtr
<
    Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>
>
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::clone() const
{
    notImplemented("quadratureNode::clone() const");
    return autoPtr
    <
        quadratureNode<weightType, abscissaType, sigmaType, vType>
    >(NULL);
}


template <class weightType, class abscissaType, class sigmaType, class vType>
Foam::tmp<vType>
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
velocityAbscissae() const
{
    if (velocityIndexes_[0] == -1)
    {
        FatalErrorInFunction
            << "Attempt to return velocity abscissa of a quadrature node" << nl
            << "    but no abscissa have dimensions of velocity. "
            << abort(FatalError);
    }

    tmp<vType> v
    (
        new vType
        (
            IOobject
            (
                IOobject::groupName("velocityAbscissa", name_),
                weight_.mesh().time().timeName(),
                weight_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            weight_.mesh(),
            dimensionedVector("zeroVelocity", dimVelocity, Zero)
        )
    );

    forAll(velocityIndexes_, cmpt)
    {
        if (velocityIndexes_[cmpt] != -1)
        {
            v.ref().replace(cmpt, abscissa_[velocityIndexes_[cmpt]]);
        }
    }
    return v;
}


template <class weightType, class abscissaType, class sigmaType, class vType>
Foam::vector
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
velocityAbscissae(const label celli) const
{
    if (velocityIndexes_[0] == -1)
    {
        FatalErrorInFunction
            << "Attempt to return velocity abscissa of a quadrature node" << nl
            << "    but no abscissa have dimensions of velocity. "
            << abort(FatalError);
    }

    vector v(Zero);

    forAll(velocityIndexes_, cmpt)
    {
        if (velocityIndexes_[cmpt] != -1)
        {
            v[cmpt] = abscissa_[velocityIndexes_[cmpt]][celli];
        }
    }
    return v;
}


template <class weightType, class abscissaType, class sigmaType, class vType>
void
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
setVelocityAbscissae(const vType& U)
{
    if (velocityIndexes_[0] == -1)
    {
        FatalErrorInFunction
            << "Attempt to set velocity abscissa of a quadrature node" << nl
            << "    but no abscissa have dimensions of velocity. "
            << abort(FatalError);
    }

    forAll(velocityIndexes_, cmpt)
    {
        if (velocityIndexes_[cmpt] != -1)
        {
            abscissa_[velocityIndexes_[cmpt]] = U.component(cmpt);
        }
    }
}


template <class weightType, class abscissaType, class sigmaType, class vType>
void
Foam::quadratureNode<weightType, abscissaType, sigmaType, vType>::
setVelocityAbscissae(const label celli, const vector& U)
{
    if (velocityIndexes_[0] == -1)
    {
        FatalErrorInFunction
            << "Attempt to set velocity abscissa of a quadrature node" << nl
            << "    but no abscissa have dimensions of velocity. "
            << abort(FatalError);
    }

    forAll(velocityIndexes_, cmpt)
    {
        if (velocityIndexes_[cmpt] != -1)
        {
            abscissa_[velocityIndexes_[cmpt]][celli] = component(U, cmpt);
        }
    }
}
// ************************************************************************* //
