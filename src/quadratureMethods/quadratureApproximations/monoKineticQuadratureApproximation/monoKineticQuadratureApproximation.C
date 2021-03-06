/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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

#include "monoKineticQuadratureApproximation.H"
#include "fixedValueFvPatchFields.H"
#include "cyclicFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "emptyFvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "Vandermonde.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::monoKineticQuadratureApproximation::monoKineticQuadratureApproximation
(
    const word& name,
    const fvMesh& mesh,
    const word& support
)
:
    quadratureApproximation(name, mesh, support, 1),
    U_
    (
        mesh_.lookupObject<volVectorField>
        (
            IOobject::groupName("U", name_)
        )
    ),
    nNodes_(nodes_().size()),
    velocityMoments_(max(nNodes_, 2)),
    velocityAbscissae_(nNodes_),
    nodesNei_(),
    velocitiesNei_(nNodes_),
    nodesOwn_(),
    velocitiesOwn_(nNodes_),
    minM0_(readScalar((*this).subDict("residuals").lookup("minM0"))),
    minM1_(readScalar((*this).subDict("residuals").lookup("minM1")))
{
    //  Set boundary cconditions for velocity abscissae based on
    //  mean velocity
    wordList UTypes(U_.boundaryField().types());

    forAll(U_.boundaryField(), i)
    {
        if
        (
            isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
         || isA<slipFvPatchVectorField>(U_.boundaryField()[i])
         || isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
        )
        {
            UTypes[i] = fixedValueFvPatchVectorField::typeName;
        }
        else if
        (
            isA<directionMixedFvPatchVectorField>(U_.boundaryField()[i])
        )
        {
            UTypes[i] = zeroGradientFvPatchVectorField::typeName;
        }
    }

    forAll(velocityMoments_, mi)
    {
        velocityMoments_.set
        (
            mi,
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "Up",
                        IOobject::groupName
                        (
                            Foam::name(mi),
                            "air"
                        )
                    ),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                moments_[mi]*U_,
                UTypes
            )
        );
    }
    nodesNei_ = autoPtr<PtrList<surfaceScalarNode>>
    (
        new PtrList<surfaceScalarNode>(nNodes_)
    );
    nodesOwn_ = autoPtr<PtrList<surfaceScalarNode>>
    (
        new PtrList<surfaceScalarNode>(nNodes_)
    );

    PtrList<surfaceScalarNode>& nodesNei = nodesNei_();
    PtrList<surfaceScalarNode>& nodesOwn = nodesOwn_();


    // Populating interpolated nodes
    forAll(nodes_(), nodei)
    {
        velocityAbscissae_.set
        (
            nodei,
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "U",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(nodei)
                        )
                    ),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                U_,
                U_.boundaryField().types()
            )
        );
        nodesNei.set
        (
            nodei,
            new surfaceScalarNode
            (
                "node" + Foam::name(nodei) + "Nei",
                name_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions(),
                false,
                0
            )
        );
        velocitiesNei_.set
        (
            nodei,
            new surfaceVectorField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "UNei",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(nodei)
                        )
                    ),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::interpolate(U_)
            )
        );

        nodesOwn.set
        (
            nodei,
            new surfaceScalarNode
            (
                "node" + Foam::name(nodei) + "Own",
                name_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions(),
                false,
                0
            )
        );
        velocitiesOwn_.set
        (
            nodei,
            new surfaceVectorField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "UOwn",
                        IOobject::groupName
                        (
                            name_,
                            Foam::name(nodei)
                        )
                    ),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::interpolate(U_)
            )
        );
    }

    updateAllQuadrature();
    interpolateNodes();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::monoKineticQuadratureApproximation::~monoKineticQuadratureApproximation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::monoKineticQuadratureApproximation::interpolateNodes()
{
    surfaceScalarField nei
    (
        IOobject
        (
            "nei",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("nei", dimless, -1.0)
    );

    surfaceScalarField own
    (
        IOobject
        (
            "own",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("own", dimless, 1.0)
    );

    const PtrList<volScalarNode>& nodes = nodes_();
    PtrList<surfaceScalarNode>& nodesNei = nodesNei_();
    PtrList<surfaceScalarNode>& nodesOwn = nodesOwn_();

    forAll(nodes, nodei)
    {
        const volScalarNode& node(nodes[nodei]);
        surfaceScalarNode& nodeNei(nodesNei[nodei]);
        surfaceScalarNode& nodeOwn(nodesOwn[nodei]);

        nodeOwn.primaryWeight() =
            fvc::interpolate(node.primaryWeight(), own, "reconstruct(weight)");

        nodeOwn.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                own,
                "reconstruct(abscissa)"
            );

        velocitiesOwn_[nodei] =
            fvc::interpolate(velocityAbscissae_[nodei], own, "reconstruct(U)");


        nodeNei.primaryWeight() =
            fvc::interpolate(node.primaryWeight(), nei, "reconstruct(weight)");

        nodeNei.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                nei,
                "reconstruct(abscissa)"
            );

        velocitiesNei_[nodei] =
            fvc::interpolate(velocityAbscissae_[nodei], nei, "reconstruct(U)");
    }
}

void Foam::monoKineticQuadratureApproximation::updateBoundaryVelocities()
{
    const volScalarField& m0 = moments_[0];

    // Update boundary node velocities
    const volScalarField::Boundary& m0Bf = m0.boundaryField();
    forAll(m0Bf, patchi)
    {
        const fvPatchScalarField& m0Patch = m0Bf[patchi];

        forAll(m0Patch, facei)
        {
            //- number of nodes with a non-negligable number of bubble or mass
            label nNonZeroNodes = 0;
            boolList nonZeroNodes(nNodes_, false);

            // Check if moment.0 is large enough to be meaningfull
            if (m0Bf[patchi][facei] > minM0_)
            {
                forAll(nodes_(), nodei)
                {
                    // Check if bubble moments are large enough.
                    //  If yes make matricies 1 component larger,
                    //  if no the rest of the nodes are assumed to
                    //  be too small as well.
                    //  This is done to avoid a divide by 0 error,
                    //  and to reduce unneeded computation time
                    if
                    (
                        nodes_()[nodei].primaryWeight().boundaryField()[patchi][facei]
                      > minM0_
                     &&
                        nodes_()[nodei].primaryAbscissa().boundaryField()[patchi][facei]
                      > SMALL
                    )
                    {
                        nonZeroNodes[nodei] = true;
                        nNonZeroNodes++;
                    }
                }
            }

            if (nNonZeroNodes == 1)
            {
                label index = -1;
                forAll(nonZeroNodes, nodei)
                {
                    if (nonZeroNodes[nodei])
                    {
                        index = nodei;
                        break;
                    }
                }
                velocityAbscissae_[index].boundaryFieldRef()[patchi][facei] =
                    velocityMoments_[1].boundaryField()[patchi][facei]
                   /(
                        nodes_()[index].primaryWeight().boundaryField()[patchi][facei]
                       *nodes_()[index].primaryAbscissa().boundaryField()[patchi][facei]
                    );
            }
            else if (nNonZeroNodes > 1)
            {
                // Create invV and invR matrices outside of cmptI loop to save time
                scalarSquareMatrix invR(nNonZeroNodes, 0.0);
                scalarDiagonalMatrix x(nNonZeroNodes, 0.0);
                label nodej = 0;
                for (label nodei = 0; nodei < nNodes_; nodei++)
                {
                    if (nonZeroNodes[nodei])
                    {
                        x[nodej] =
                            nodes_()[nodei].primaryAbscissa().boundaryField()[patchi][facei];

                        invR[nodej][nodej] =
                            1.0
                           /nodes_()[nodei].primaryWeight().boundaryField()[patchi][facei];

                        nodej++;
                    }
                }

                // Invert V martix and create invVR matrix
                Vandermonde V(x);
                scalarRectangularMatrix invVR = invR*V.inv();

                // Loop over all components of U_{\alpha}
                for (label cmpti = 0; cmpti < vector::nComponents; cmpti++)
                {
                    scalarRectangularMatrix Upcmpt(nNonZeroNodes, 1, 0.0);
                    label nodej = 0;
                    for (label nodei = 0; nodei < nNodes_; nodei++)
                    {
                        if (nonZeroNodes[nodei])
                        {
                            Upcmpt[nodej][0] =
                                velocityMoments_[nodei].boundaryField()[patchi][facei].component(cmpti);

                            nodej++;
                        }
                    }

                    // Compute U_{\alpha} cmptI component using invVR matrix
                    scalarRectangularMatrix Ucmpt = invVR*Upcmpt;
                    nodej = 0;
                    for (label nodei = 0; nodei < nNodes_; nodei++)
                    {
                        if (nonZeroNodes[nodei])
                        {
                            velocityAbscissae_[nodei].boundaryFieldRef()[patchi][facei].component(cmpti) = Ucmpt[nodej][0];
                            nodej++;
                        }
                    }
                }
            }

            // Set nodes with very small bubble mass or number to zero velocity
            for (label nodei = 0; nodei < nNodes_; nodei++)
            {
                if (!nonZeroNodes[nodei])
                {
                    velocityAbscissae_[nodei].boundaryFieldRef()[patchi][facei]
                    = U_.boundaryField()[patchi][facei];
                }
            }
        }
    }
}

void Foam::monoKineticQuadratureApproximation::updateAllQuadrature()
{
    const volScalarField& m0 = moments_[0];
    const volScalarField::Boundary m0Bf = m0.boundaryField();

    // Check for small moments at cell centers
    forAll(m0, celli)
    {
        //- Make sure moments are below 0 before checking if they
        //  are small enough to be neglected
        if
        (
            m0[celli] < 0
         && mag(m0[celli]) < minM0_
        )
        {
            forAll(moments_, mi)
            {
                moments_[mi][celli] = 0.0;
            }
        }
        else if
        (
            moments_[1][celli] < 0
         && mag(moments_[1][celli]) < minM1_
        )
        {
            for (label mi = 1; mi < nMoments_; mi++)
            {
                moments_[mi][celli] = 0.0;
            }
        }
    }

    // Check for small moments on boundaries
    forAll(m0Bf, patchi)
    {
        forAll(m0Bf[patchi], facei)
        {
            if
            (
                m0Bf[patchi][facei] < 0
             && mag(m0Bf[patchi][facei]) < minM0_
            )
            {
                forAll(moments_, mi)
                {
                    moments_[mi].boundaryFieldRef()[patchi][facei] = 0.0;
                }
            }
            else if
            (
                moments_[1].boundaryField()[patchi][facei] < 0
             && mag(moments_[1].boundaryField()[patchi][facei]) < minM1_
            )
            {
                for (label mi = 1; mi < nMoments_; mi++)
                {
                    moments_[mi].boundaryFieldRef()[patchi][facei] = 0.0;
                }
            }
        }
    }

    updateQuadrature();

    updateVelocities();
    updateBoundaryVelocities();

    forAll(nodes_(), nodei)
    {
        velocityAbscissae_[nodei].correctBoundaryConditions();
    }

    updateAllMoments();
}

void Foam::monoKineticQuadratureApproximation::updateVelocities()
{
    const volScalarField& m0 = moments_[0];

    forAll(m0, celli)
    {
        //- number of nodes with a non-negligable number of bubble or mass
        label nNonZeroNodes = 0;
        boolList nonZeroNodes(nNodes_, false);

        // Check if moment.0 is large enough to be meaningful
        if (m0[celli] > minM0_)
        {
            forAll(nodes_(), nodei)
            {
                // Check if size moments are large enough.
                //  If yes make matricies 1 component larger.
                //  This is done to avoid a divide by 0 error,
                //  and to reduce unneeded computation time
                if
                (
                    nodes_()[nodei].primaryWeight()[celli] > minM0_
                 && nodes_()[nodei].primaryAbscissa()[celli] > SMALL
                )
                {
                    nonZeroNodes[nodei] = true;
                    nNonZeroNodes++;
                }
            }
        }

        if (nNonZeroNodes == 1)
        {
            label index = -1;
            forAll(nonZeroNodes, nodei)
            {
                if (nonZeroNodes[nodei])
                {
                    index = nodei;
                    break;
                }
            }
            velocityAbscissae_[index][celli] =
                velocityMoments_[1][celli]
               /(
                   nodes_()[index].primaryWeight()[celli]
                  *nodes_()[index].primaryAbscissa()[celli]
                );
        }
        else if (nNonZeroNodes > 1)
        {
            // Create invV and invR matrices outside of cmptI loop to save time
            scalarSquareMatrix invR(nNonZeroNodes, 0.0);
            scalarDiagonalMatrix x(nNonZeroNodes, 0.0);
            label nodej = 0;
            for (label nodei = 0; nodei < nNodes_; nodei++)
            {
                if (nonZeroNodes[nodei])
                {
                    x[nodej] =
                        nodes_()[nodei].primaryAbscissa()[celli];

                    invR[nodej][nodej] =
                        1.0/nodes_()[nodei].primaryWeight()[celli];

                    nodej++;
                }
            }

            // Invert V martix and create invVR matrix
            Vandermonde V(x);
            scalarRectangularMatrix invVR = invR*V.inv();

            // Loop over all components of U_{\alpha}
            for (label cmpti = 0; cmpti < vector::nComponents; cmpti++)
            {
                scalarRectangularMatrix Upcmpt(nNonZeroNodes, 1, 0.0);
                label nodej = 0;
                for (label nodei = 0; nodei < nNodes_; nodei++)
                {
                    if (nonZeroNodes[nodei])
                    {
                        Upcmpt[nodej][0] =
                            velocityMoments_[nodei][celli].component(cmpti);
                        nodej++;
                    }
                }

                // Compute U_{\alpha} cmptI component using invVR matrix
                scalarRectangularMatrix Ucmpt = invVR*Upcmpt;
                nodej = 0;
                for (label nodei = 0; nodei < nNodes_; nodei++)
                {
                    if (nonZeroNodes[nodei])
                    {
                        velocityAbscissae_[nodei][celli].component(cmpti) =
                            Ucmpt[nodej][0];
                        nodej++;
                    }
                }
            }
        }

        // Set nodes with very small bubble mass or number to mean bubble
        // velocity
        for (label nodei = 0; nodei < nNodes_; nodei++)
        {
            if (!nonZeroNodes[nodei])
            {
                velocityAbscissae_[nodei][celli] = U_[celli];
            }
        }
    }
}

void Foam::monoKineticQuadratureApproximation::updateVelocityMoments()
{
    // Update velocity moments
    forAll(velocityMoments_, mi)
    {
        velocityMoments_[mi] =
            dimensionedVector
            (
                "zero",
                velocityMoments_[mi].dimensions(),
                Zero
            );

        if (mi == 0)
        {
            forAll(nodes_(), nodei)
            {
                velocityMoments_[mi] +=
                    nodes_()[nodei].primaryWeight()*velocityAbscissae_[nodei];
            }

            velocityMoments_[mi].correctBoundaryConditions();
        }
        else
        {
            forAll(nodes_(), nodei)
            {
                velocityMoments_[mi] +=
                    nodes_()[nodei].primaryWeight()
                   *pow(nodes_()[nodei].primaryAbscissa(), mi)
                   *velocityAbscissae_[nodei];
            }

            velocityMoments_[mi].correctBoundaryConditions();
        }
    }
}



void Foam::monoKineticQuadratureApproximation::updateAllMoments()
{
    // Update size moments
    updateMoments();

    // Update velocity moments
    updateVelocityMoments();
}
// ************************************************************************* //
