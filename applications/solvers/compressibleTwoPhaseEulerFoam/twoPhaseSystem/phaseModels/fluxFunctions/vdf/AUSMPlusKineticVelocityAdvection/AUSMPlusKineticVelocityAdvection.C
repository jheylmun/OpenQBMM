/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Alberto Passalacqua
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

#include "AUSMPlusKineticVelocityAdvection.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityAdvection
{
    defineTypeNameAndDebug(AUSMPlusKinetic, 0);

    addToRunTimeSelectionTable
    (
        velocityMomentAdvection,
        AUSMPlusKinetic,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * Private Data Members * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::velocityAdvection::AUSMPlusKinetic::M1
(
    const surfaceScalarField& M,
    const label sign
)
{
    return 0.5*(M + sign*mag(M));
}

Foam::tmp<Foam::surfaceScalarField>
Foam::velocityAdvection::AUSMPlusKinetic::M2
(
    const surfaceScalarField& M,
    const label sign
)
{
    return sign*0.25*sqr(M + sign);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityAdvection::AUSMPlusKinetic::AUSMPlusKinetic
(
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    const word& support
)
:
    velocityMomentAdvection(dict, quadrature, support),
    nodes_(quadrature.nodes()),
    nodesNei_(),
    nodesOwn_(),
    momentsOwn_(moments_.size()),
    momentsNei_(moments_.size()),
    momentsf_(moments_.size()),
    rho_
    (
        moments_[0].mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("thermo:rho", moments_[0].group())
        )
    ),
    pFric_
    (
        moments_[0].mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("Pfric", moments_[0].group())
        )
    ),
    alphaMax_("alphaMax", dimless, dict),
    alphaMinFriction_("alphaMinFriction", dimless, dict),
    cutOffMa_("small", dimless, 1e-10),
    residualRho_("small", dimDensity, 1e-10),
    residualU_("small", dimVelocity, 1e-10),
    fa_(dict.lookupOrDefault("fa", 1.0)),
    D_(dict.lookupOrDefault("D", 1.0)),
    xi_(3.0/16.0*(5.0*sqr(fa_) - 4.0))
{
    nodesNei_ = autoPtr<PtrList<surfaceVectorNode> >
    (
        new PtrList<surfaceVectorNode>(nodes_.size())
    );

    nodesOwn_ = autoPtr<PtrList<surfaceVectorNode> >
    (
        new PtrList<surfaceVectorNode>(nodes_.size())
    );

    PtrList<surfaceVectorNode>& nodesNei = nodesNei_();
    PtrList<surfaceVectorNode>& nodesOwn = nodesOwn_();

    // Populating nodes and interpolated nodes
    forAll(nodes_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        nodesNei.set
        (
            nodei,
            new surfaceVectorNode
            (
                "nodeNei" + mappedList<scalar>::listToWord(nodeIndex),
                name_,
                moments_[0].mesh(),
                dimless,
                dimVelocity,
                false
            )
        );

        nodesOwn.set
        (
            nodei,
            new surfaceVectorNode
            (
                "nodeOwn" + mappedList<scalar>::listToWord(nodeIndex),
                name_,
                moments_[0].mesh(),
                dimless,
                dimVelocity,
                false
            )
        );
    }
    forAll(momentsNei_, momenti)
    {
        momentsNei_.set
        (
            momenti,
            new surfaceScalarField
            (
                "momentNei" + Foam::name(momenti),
                fvc::interpolate(moments_[momenti])
            )
        );

        momentsOwn_.set
        (
            momenti,
            new surfaceScalarField
            (
                "momentOwn" + Foam::name(momenti),
                fvc::interpolate(moments_[momenti])
            )
        );
        momentsf_.set
        (
            momenti,
            new surfaceScalarField
            (
                "momentf" + Foam::name(momenti),
                fvc::interpolate(moments_[momenti])
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityAdvection::AUSMPlusKinetic::~AUSMPlusKinetic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityAdvection::AUSMPlusKinetic::interpolateNodes()
{
    PtrList<surfaceVectorNode>& nodesNei = nodesNei_();
    PtrList<surfaceVectorNode>& nodesOwn = nodesOwn_();

    forAll(nodes_, nodei)
    {
        const volVectorNode& node(nodes_[nodei]);
        surfaceVectorNode& nodeNei(nodesNei[nodei]);
        surfaceVectorNode& nodeOwn(nodesOwn[nodei]);

        nodeOwn.primaryWeight() =
            fvc::interpolate
            (
                node.primaryWeight(),
                own_,
                "reconstruct(weight)"
            );

        nodeOwn.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                own_,
                "reconstruct(U)"
            );

        nodeNei.primaryWeight() =
            fvc::interpolate
            (
                node.primaryWeight(),
                nei_,
                "reconstruct(weight)"
            );

        nodeNei.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                nei_,
                "reconstruct(U)"
            );
    }

    forAll(momentOrders_, mi)
    {
        momentsOwn_[mi] = dimensionedScalar("0", moments_[mi].dimensions(), 0.0);
        momentsNei_[mi] = dimensionedScalar("0", moments_[mi].dimensions(), 0.0);

        forAll(nodes_, nodei)
        {
            const surfaceVectorNode& nodeNei = nodesNei_()[nodei];
            const surfaceVectorNode& nodeOwn = nodesOwn_()[nodei];

            surfaceScalarField mOwn(nodeOwn.primaryWeight());
            surfaceScalarField mNei(nodeNei.primaryWeight());
            const surfaceVectorField& UOwn = nodeOwn.primaryAbscissa();
            const surfaceVectorField& UNei = nodeNei.primaryAbscissa();

            forAll(momentOrders_[mi], cmpti)
            {
                label mCmpt = momentOrders_[mi][cmpti];
                surfaceScalarField abscissaCmptOwn(UOwn.component(cmpti));
                surfaceScalarField abscissaCmptNei(UNei.component(cmpti));

                {
                    tmp<surfaceScalarField> mPow
                    (
                        mOwn*pow(abscissaCmptOwn, mCmpt)
                    );
                    mOwn.dimensions().reset(mPow().dimensions());
                    mOwn == mPow;
                }
                {
                    tmp<surfaceScalarField> mPow
                    (
                        mNei*pow(abscissaCmptNei, mCmpt)
                    );
                    mNei.dimensions().reset(mPow().dimensions());
                    mNei == mPow;
                }
            }

            momentsOwn_[mi] += mOwn;
            momentsNei_[mi] += mNei;
        }
    }
}

void Foam::velocityAdvection::AUSMPlusKinetic::updateWallCollisions()
{
    const fvMesh& mesh = own_.mesh();

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& currPatch = mesh.boundary()[patchi];
        if (isA<wallFvPatch>(currPatch))
        {
            const vectorField& bfSf(mesh.Sf().boundaryField()[patchi]);
            vectorField bfNorm(bfSf/mag(bfSf));

            forAll(nodes_, nodei)
            {
                const volVectorNode& node = nodes_[nodei];
                surfaceVectorNode& nodeNei(nodesNei_()[nodei]);
                surfaceVectorNode& nodeOwn(nodesOwn_()[nodei]);

                const volScalarField& weight = node.primaryWeight();
                surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
                surfaceScalarField& weightNei = nodeNei.primaryWeight();
                const volVectorField& U = node.primaryAbscissa();
                surfaceVectorField& UOwn = nodeOwn.primaryAbscissa();
                surfaceVectorField& UNei = nodeNei.primaryAbscissa();

                scalarField& bfwOwn = weightOwn.boundaryFieldRef()[patchi];
                scalarField& bfwNei = weightNei.boundaryFieldRef()[patchi];
                vectorField& bfUOwn = UOwn.boundaryFieldRef()[patchi];
                vectorField& bfUNei = UNei.boundaryFieldRef()[patchi];

                forAll(currPatch, facei)
                {
                    label faceCelli = currPatch.faceCells()[facei];

                    bfwOwn[facei] = weight[faceCelli];
                    bfUOwn[facei] = U[faceCelli];

                    bfwNei[facei] = bfwOwn[facei];
                    bfUNei[facei] =
                        bfUOwn[facei]
                      - (1.0 + this->ew_)*(bfUOwn[facei] & bfNorm[facei])
                       *bfNorm[facei];
                }
            }
        }
    }
}


Foam::scalar
Foam::velocityAdvection::AUSMPlusKinetic::realizableCo() const
{
    return 1.0;
}

Foam::scalar Foam::velocityAdvection::AUSMPlusKinetic::CoNum() const
{
    scalar CoNum = 0.0;
    const fvMesh& mesh = own_.mesh();
    forAll(nodes_, nodei)
    {
        CoNum =
            max
            (
                CoNum,
                0.5*gMax
                (
                    fvc::surfaceSum
                    (
                        mag(fvc::flux(nodes_[nodei].primaryAbscissa()))
                    )().primitiveField()/mesh.V().field()
                )*mesh.time().deltaTValue()
            );
    }
    return CoNum;
}

void Foam::velocityAdvection::AUSMPlusKinetic::update()
{
    const fvMesh& mesh = own_.mesh();
    dimensionedScalar zeroPhi("zero", dimVolume/dimTime, 0.0);
    surfaceVectorField Sf(mesh.Sf());
    surfaceScalarField magSf(mesh.magSf());
    surfaceVectorField normal(Sf/magSf);

    // Interpolate weights and abscissae
    interpolateNodes();

    // Set velocities at boundaries for rebounding
    updateWallCollisions();

    // Zero moment flux
    forAll(momentFluxes_, mi)
    {
        momentFluxes_[mi] =
            dimensionedScalar
            (
                "0",
                momentFluxes_[mi].dimensions(),
                0.0
            );
    }

    //- Volume fraction of own and nei faces
    const surfaceScalarField& m0Own = momentsOwn_[0];
    const surfaceScalarField& m0Nei = momentsNei_[0];

    //- Interpolated densities
    surfaceScalarField rhoOwn
    (
        fvc::interpolate(rho_, own_, interpScheme(rho_.name()))
    );
    surfaceScalarField rhoNei
    (
        fvc::interpolate(rho_, nei_, interpScheme(rho_.name()))
    );

    //- Interpolated frictional pressure
    surfaceScalarField pFricOwn
    (
        fvc::interpolate(pFric_, own_, interpScheme(pFric_.name()))
    );
    surfaceScalarField pFricNei
    (
        fvc::interpolate(pFric_, nei_, interpScheme(pFric_.name()))
    );

    //- Interpolated speed of sounds
    const volScalarField& c = moments_[0].mesh().lookupObject<volScalarField>
    (
        IOobject::groupName("c", moments_[0].group())
    );
    surfaceScalarField cOwn
    (
        fvc::interpolate(c, own_, interpScheme(c.name()))
    );
    surfaceScalarField cNei
    (
        fvc::interpolate(c, nei_, interpScheme(c.name()))
    );
    surfaceScalarField c12
    (
        sqrt
        (
            (m0Own*sqr(cOwn) + m0Nei*sqr(cNei))/max(m0Own + m0Nei, 1e-6)
        ) + residualU_
    );

    surfaceScalarField zeta
    (
        max
        (
            (max(m0Own, m0Nei) - alphaMinFriction_)
           /(alphaMax_ - alphaMinFriction_),
            0.0
        )
    );
    surfaceScalarField G(max(2.0*(1.0 - D_*sqr(zeta)), 0.0));

    surfaceScalarField Kp(0.25 + 0.75*(1.0 - G/2.0));
    surfaceScalarField Ku(0.75 + 0.25*(1.0 - G/2.0));
    surfaceScalarField sigma(0.75*G/2.0);

    forAll(nodes_, nodei)
    {
        const surfaceVectorNode& nodeNei(nodesNei_()[nodei]);
        const surfaceVectorNode& nodeOwn(nodesOwn_()[nodei]);

        const surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        const surfaceScalarField& weightNei = nodeNei.primaryWeight();

        const surfaceVectorField& UOwn = nodeOwn.primaryAbscissa();
        const surfaceVectorField& UNei = nodeNei.primaryAbscissa();
        surfaceScalarField UvOwn(UOwn & normal);
        surfaceScalarField UvNei(UNei & normal);

        surfaceScalarField MaOwn(UvOwn/c12);
        surfaceScalarField magMaOwn(mag(MaOwn));
        surfaceScalarField MaNei(UvNei/c12);
        surfaceScalarField magMaNei(mag(MaNei));

        surfaceScalarField MaBarSqr
        (
            (sqr(UvOwn) + sqr(UvNei))/(2.0*sqr(c12))
        );

        surfaceScalarField Ma4Own
        (
            pos0(magMaOwn - 1)*M1(MaOwn, 1)
          + neg(magMaOwn - 1)*M2(MaOwn, 1)*(1.0 - 16.0*beta_*M2(MaOwn, -1))
        );
        surfaceScalarField Ma4Nei
        (
            pos0(magMaNei - 1)*M1(MaNei, -1)
          + neg(magMaNei - 1)*M2(MaNei, -1)*(1.0 + 16.0*beta_*M2(MaNei, 1))
        );

        surfaceScalarField Ma12
        (
            Ma4Own + Ma4Nei
          - 2.0*Kp/fa_*max(1.0 - sigma*MaBarSqr, 0.0)*(pFricNei - pFricOwn)
           /((weightOwn*rhoOwn + weightNei*rhoOwn + residualRho_)*sqr(c12))
        );

        surfaceScalarField p5Own
        (
            pos0(magMaOwn - 1)*pos(MaOwn)
          + neg(magMaOwn - 1)
           *(M2(MaOwn, 1)*(2.0 - MaOwn) - 16.0*xi_*MaOwn*M2(MaOwn, -1))
        );

        surfaceScalarField p5Nei
        (
            pos0(magMaNei - 1)*neg(MaNei)
          + neg(magMaNei - 1)
           *(M2(MaNei, -1)*(-2.0 - MaNei) + 16.0*xi_*MaNei*M2(MaNei, 1))
        );

        surfaceScalarField F
        (
            (c12 - residualU_)
           *(
                1.0
              + mag(Ma12)*(1.0 - G/2.0)
            )
           *max(m0Own, m0Nei)/alphaMax_
           *(weightOwn - weightNei)/2.0
        );

        surfaceScalarField mDot
        (
            F
          + c12*Ma12
           *(
                pos(Ma12)*weightOwn
              + neg0(Ma12)*weightNei
            )
        );

        surfaceVectorField gradp
        (
            (
               -Ku*fa_*(c12 - residualU_)*p5Own*p5Nei
               *(weightOwn + weightNei)*(UvNei - UvOwn)
              + p5Own*pFricOwn/rhoOwn  + p5Nei*pFricNei/rhoNei
            )*Sf
        );

        surfaceScalarField m0Flux(magSf*mDot);
        forAll(momentFluxes_, mi)
        {
            const labelList& momentOrder = momentOrders_[mi];
            if (max(momentOrder) == 0)
            {
                momentFluxes_[mi] += m0Flux;
                continue;
            }
            surfaceScalarField momentOwnFlux(m0Flux);
            surfaceScalarField momentNeiFlux(m0Flux);
            surfaceVectorField pCmptOwn(weightOwn*gradp);
            surfaceVectorField pCmptNei(weightNei*gradp);

            forAll(momentOrder, cmpti)
            {
                const label cmptMomentOrder = momentOrder[cmpti];
                surfaceScalarField absCmptOwn
                (
                    pow(UOwn.component(cmpti), cmptMomentOrder)
                );
                surfaceScalarField absCmptNei
                (
                    pow(UNei.component(cmpti), cmptMomentOrder)
                );
                tmp<surfaceScalarField> mOwnPow = momentOwnFlux*absCmptOwn;
                momentOwnFlux.dimensions().reset(mOwnPow().dimensions());
                momentOwnFlux == mOwnPow();

                tmp<surfaceScalarField> mNeiPow = momentNeiFlux*absCmptNei;
                momentNeiFlux.dimensions().reset(mNeiPow().dimensions());
                momentNeiFlux == mNeiPow;

//                 if (cmptMomentOrder > 0)
//                 {
//                     forAll(momentOrder, cmptj)
//                     {
//                         tmp<surfaceScalarField> pOwnPow;
//                         tmp<surfaceScalarField> pNeiPow;
//
//                         if (cmpti == cmptj)
//                         {
//                             pOwnPow =
//                                 pCmptOwn.component(cmpti)
//                                *pow
//                                 (
//                                     UOwn.component(cmpti),
//                                     max(cmptMomentOrder - 1, 0)
//                                 );
//
//                             pNeiPow =
//                                 pCmptNei.component(cmpti)
//                                *pow
//                                 (
//                                     UNei.component(cmpti),
//                                     max(cmptMomentOrder - 1, 0)
//                                 );
//
//                         }
//                         else
//                         {
//                             pOwnPow = pCmptOwn.component(cmpti)*absCmptOwn;
//                             pNeiPow = pCmptNei.component(cmpti)*absCmptNei;
//                         }
//                         pCmptOwn.replace(cmptj, pOwnPow);
//                         pCmptNei.replace(cmptj, pNeiPow);
//                     }
//                 }
//                 pCmptOwn.dimensions().reset(momentOwnFlux.dimensions());
//                 pCmptNei.dimensions().reset(momentNeiFlux.dimensions());
            }

            momentFluxes_[mi] +=
                momentOwnFlux*pos0(mDot) + momentNeiFlux*neg0(mDot);
//             divMoments_[divi] +=
//                 fvc::surfaceIntegrate
//                 (
//                     (momentOwnFlux + (pCmptOwn & vector::one))*pos0(mDot)
//                   + (momentNeiFlux + (pCmptNei & vector::one))*neg0(mDot)
//                 );
        }
    }
}

void Foam::velocityAdvection::AUSMPlusKinetic::update
(
    const surfaceScalarField& phi,
    const bool wallCollisions
)
{
    NotImplemented;
}

void Foam::velocityAdvection::AUSMPlusKinetic::update
(
    const mappedPtrList<volVectorField>& Us,
    const bool wallCollisions
)
{
    NotImplemented;
}

// ************************************************************************* //
