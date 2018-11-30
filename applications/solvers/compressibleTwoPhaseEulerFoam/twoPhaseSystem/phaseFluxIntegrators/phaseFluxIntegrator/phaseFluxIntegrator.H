/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Jeff Heylmun
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

Class
    Foam::phaseFluxIntegrator

Description
    Class used to handle encoding/decoding of primative and conservative
    variables used in compressible flows. Also includes funcitons to calculate
    the conservative fluxes

SourceFiles
    phaseFluxIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef phaseFluxIntegrator_H
#define phaseFluxIntegrator_H

#include "IOdictionary.H"
#include "phaseModel.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fluxIntegrator Declaration
\*---------------------------------------------------------------------------*/

class phaseFluxIntegrator
{
protected:

    // Protected data

        //- phase 1
        phaseModel& phase1_;

        //- phase 2
        phaseModel& phase2_;

        //- Reference to the dominate volume fraction gradient
        const volVectorField& gradAlpha_;

        //- Number of steps
        virtual label nSteps() const = 0;

        //- Conserved variable coeffients
        virtual List<scalarList> coeffs() const = 0;

        //- Flux coefficients
        virtual List<scalarList> Fcoeffs() const = 0;

        //- Return fields that need to be set
        void setCoeffs(boolList& setFields, boolList& setDeltas);

public:

    //- Runtime type information
    TypeName("phaseFluxIntegrator");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseFluxIntegrator,
            dictionary,
            (
                phaseModel& phase1,
                phaseModel& phase2
            ),
            (phase1, phase2)
        );


    // Constructor
    phaseFluxIntegrator(phaseModel& phase1, phaseModel& phase2);

    // Selector
    static autoPtr<phaseFluxIntegrator> New
    (
        phaseModel& phase1,
        phaseModel& phase2
    );


    //- Destructor
    virtual ~phaseFluxIntegrator();


    // Member Functions

        //- Transport conservative variables
        virtual void integrateFluxes
        (
            const dimensionedVector& g,
            volVectorField& Ui,
            volScalarField& pi
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
