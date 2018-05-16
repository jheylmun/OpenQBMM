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

#include "velocityPDFTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::velocityPDFTransportModel::velocityPDFTransportModel
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh,
    const word& support
)
:
    PDFTransportModel(name, dict, mesh, support),
    momentAdvection_
    (
        velocityMomentAdvection::New
        (
            quadrature_.subDict("momentAdvection"),
            quadrature_,
            support
        )
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::velocityPDFTransportModel::~velocityPDFTransportModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::velocityPDFTransportModel::solve()
{
    momentAdvection_().update();

    // List of moment transport equations
    PtrList<fvScalarMatrix> momentEqns(this->quadrature_.nMoments());

    // Solve moment transport equations
    forAll(this->quadrature_.moments(), momenti)
    {
        volMoment& m = this->quadrature_.moments()[momenti];
        momentEqns.set
        (
            momenti,
            new fvScalarMatrix
            (
                fvm::ddt(m)
              + momentAdvection_().divMoments()[momenti]
            )
        );
    }

    if (collision())
    {
        if (this->solveODESource_)
        {
            explicitMomentSource();
        }
        else
        {
            updateImplicitCollisionSource();
        }

        forAll(this->quadrature_.moments(), mEqni)
        {
            volMoment& m = this->quadrature_.moments()[mEqni];

            if (this->solveODESource_)
            {
                momentEqns[mEqni] -= fvc::ddt(m);
            }
            else
            {
                // Solve moment transport excluding collisions
                momentEqns[mEqni].relax();
                momentEqns[mEqni].solve();

                //  Set moments.oldTime to moments transport is not neglected due to
                //  large collision source terms
                m.oldTime() = m;

                // Solve collisions
                momentEqns.set
                (
                    mEqni,
                    new fvScalarMatrix
                    (
                        fvm::ddt(m)
                     ==
                        implicitCollisionSource(m)
                    )
                );
            }
        }
    }

    forAll(this->quadrature_.moments(), mEqni)
    {
        momentEqns[mEqni].relax();
        momentEqns[mEqni].solve();
    }

    this->quadrature_.updateQuadrature();
}

void Foam::PDFTransportModels::velocityPDFTransportModel::meanTransport
(
    const surfaceScalarField& phi,
    const bool wallCollisions
)
{
    Info<< "Solving mean transport" << endl;

    momentAdvection_().update(phi, wallCollisions);

    // Solve moment transport equations
    forAll(this->quadrature_.moments(), momenti)
    {
        volScalarField& m = this->quadrature_.moments()[momenti];
        fvScalarMatrix mEqn
        (
            fvm::ddt(m)
          - fvc::ddt(m)
          + momentAdvection_().divMoments()[momenti]
        );

        mEqn.relax();
        mEqn.solve();
    }
}

void Foam::PDFTransportModels::velocityPDFTransportModel::relativeTransport
(
    const mappedPtrList<volVectorField>& Vs,
    const bool wallCollisions
)
{
    Info<< "Solving relative transport" << endl;

    momentAdvection_().update(Vs, wallCollisions);

    // Solve moment transport equations
    forAll(this->quadrature_.moments(), momenti)
    {
        volScalarField& m = this->quadrature_.moments()[momenti];
        fvScalarMatrix mEqn
        (
            fvm::ddt(m)
          + momentAdvection_().divMoments()[momenti]
        );

        mEqn.relax();
        mEqn.solve();
    }
}

// ************************************************************************* //
