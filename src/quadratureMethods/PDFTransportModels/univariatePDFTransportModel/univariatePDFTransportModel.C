/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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

#include "univariatePDFTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariatePDFTransportModel
::univariatePDFTransportModel
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh,
    const surfaceScalarField& phi,
    const word& support
)
:
    PDFTransportModel(name, dict, mesh, support),
    momentAdvection_
    (
        univariateMomentAdvection::New
        (
            this->quadrature_.subDict("momentAdvection"),
            this->quadrature_,
            phi,
            support
        )
    )
{
    if (quadrature_.momentOrders()[0].size() != 1)
    {
        FatalErrorInFunction
            << "Only one dimensional distributions can be used with" << nl
            << "    univariatePDFTransportModel, but "
            << this->quadrature_.momentOrders()[0].size()
            << " dimensions have" << nl
            << "    been specified." << nl
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariatePDFTransportModel
::~univariatePDFTransportModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::univariatePDFTransportModel::solve()
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
              - momentDiffusion(m)
              ==
                implicitMomentSource(m)
            )
        );
    }

    if (this->solveODESource_)
    {
        explicitMomentSource();
    }

    forAll (momentEqns, mEqni)
    {
        volMoment& m = this->quadrature_.moments()[mEqni];

        if (this->solveODESource_)
        {
            momentEqns[mEqni] -= fvc::ddt(m);
        }

        momentEqns[mEqni].relax();
        momentEqns[mEqni].solve();
    }

    this->quadrature_.updateQuadrature();
}


// ************************************************************************* //
