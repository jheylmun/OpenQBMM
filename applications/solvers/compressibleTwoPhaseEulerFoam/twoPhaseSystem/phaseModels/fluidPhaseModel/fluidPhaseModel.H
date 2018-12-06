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

Description
    Phase model for a slave pressure, based on a 6 equation model

Class
    Foam::liquidPhaseModel

SourceFiles
    liquidPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef fluidPhaseModel_H
#define fluidPhaseModel_H

#include "phaseModel.H"
#include "fluidFluxFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class fluidPhaseModel
:
    public phaseModel
{
    // Private member data

        //- Phase specific thermodynamic pressure
        volScalarField p_;

        // Conserved variables

            //- Phasic density
            volScalarField alphaRho_;

            //- Phasic momentum
            volVectorField alphaRhoU_;

            //- Total phasic energy
            volScalarField alphaRhoE_;


        //- Conservative fluxes

            //- Mass flux
            surfaceScalarField massFlux_;

            //- Momentum flux
            surfaceVectorField momentumFlux_;

            //- Total energy flux
            surfaceScalarField energyFlux_;


        //- Fields stored for partial ode steps

            //- Changes in alpha
            PtrList<volScalarField> alphas_;
            PtrList<volScalarField> deltaAlphas_;

            //- Changes in mass
            PtrList<volScalarField> alphaRhos_;
            PtrList<volScalarField> deltaAlphaRhos_;

            //- Changes in momentum
            PtrList<volVectorField> alphaRhoUs_;
            PtrList<volVectorField> deltaAlphaRhoUs_;

            //- Changes in energy
            PtrList<volScalarField> alphaRhoEs_;
            PtrList<volScalarField> deltaAlphaRhoEs_;

        //- Model to calculate hyperbolic phase flux
        autoPtr<fluidFluxFunction> fluxFunction_;

        //- update temperture
        void correctP();

public:

    // Constructors
        fluidPhaseModel
        (
            const twoPhaseSystem& fluid,
            const dictionary& dict,
            const word& phaseName
        );

    //- Destructor
    virtual ~fluidPhaseModel();


    // Member Functions

        //- Return that the phase uses a slave pressure
        virtual bool slavePressure() const
        {
            return false;
        }

        //- Is the phase granular
        virtual bool granular() const
        {
            return false;
        }

        //- Const access to the phase pseudo thermal energy
        virtual tmp<volScalarField> Theta() const
        {
            NotImplemented;
            return *this;
        }

        //- Non-const access to the phase pseudo thermal energy
        virtual volScalarField& ThetaRef()
        {
            NotImplemented;
            return *this;
        }

        virtual const volScalarField& p() const
        {
            return p_;
        }

        //- Non-const access to the phase pressure
        virtual volScalarField& pRef()
        {
            return p_;
        }

        //- Speed of sound of the phase
        virtual tmp<volScalarField> c() const;

        //- Return the phase volumentric flux
        virtual tmp<surfaceScalarField> phi(const label i = -1) const
        {
            return fluxFunction_->phi();
        }

        //- Return the phase volumentric flux
        virtual tmp<surfaceScalarField> alphaPhi(const label i = -1) const
        {
            return fluxFunction_->alphaPhi();
        }

        //- Const access to phasic density
        virtual tmp<volScalarField> alphaRho(const label i = -1) const
        {
            return alphaRho_;
        }

        //- Non-const access to phasic density
        virtual volScalarField& alphaRhoRef(const label i = -1)
        {
            return alphaRho_;
        }

        //- Const access to phasic momentum
        virtual tmp<volVectorField> alphaRhoU(const label i = -1) const
        {
            return alphaRhoU_;
        }

        //- Non-const access to phasic momentum
        virtual volVectorField& alphaRhoURef(const label i = -1)
        {
            return alphaRhoU_;
        }

        //- Const access to total phasic energy
        virtual tmp<volScalarField> alphaRhoE(const label i = -1) const
        {
            return alphaRhoE_;
        }

        //- Non-const access to total phasic energy
        virtual volScalarField& alphaRhoERef(const label i = -1)
        {
            return alphaRhoE_;
        }

        //- Const access to phase internal energy
        virtual tmp<volScalarField> alphaRhoPTE(const label i = -1) const
        {
            NotImplemented;
            return alphaRhoE_;
        }

        //- Non-const access to phase internal energy
        virtual volScalarField& alphaRhoPTERef(const label i = -1)
        {
            NotImplemented;
            return alphaRhoE_;
        }

        //- Const access to phasic mass flux
        virtual tmp< surfaceScalarField> massFlux() const
        {
            return massFlux_;
        }

        //- Const access to phasic momentum flux
        virtual tmp<surfaceVectorField> momentumFlux() const
        {
            return momentumFlux_;
        }

        //- Const access to total phasic energy flux
        virtual tmp<surfaceScalarField> energyFlux() const
        {
            return energyFlux_;
        }

        //- Const access to phase internal energy flux
        virtual tmp<surfaceScalarField> PTEFlux() const
        {
            NotImplemented;
            return energyFlux_;
        }

        //- Return the gradient of pressure
        virtual tmp<surfaceScalarField> alphaf(const label i = -1) const
        {
            return fluxFunction_->alphaf();
        }

        //- Return the gradient of velocity
        virtual tmp<volTensorField> gradU(const label i = -1) const
        {
            return fluxFunction_->gradU();
        }

        //- Return the gradient of pressure
        virtual tmp<volVectorField> gradp() const
        {
            return fluxFunction_->gradp();
        }

        //- Advect with conservative fluxes
        virtual void advect
        (
            const label stepi,
            const scalarList& coeffs,
            const scalarList& Fcoeffs,
            const dimensionedScalar& deltaT,
            const dimensionedVector& g,
            const volVectorField& F,
            const volVectorField& Ui,
            const volScalarField& pi
        );

        //- Update hyperbolic fluxes
        virtual void updateFluxes();

        //- Update hyperbolic fluxes using a set interplated volume fraction
        virtual void updateFluxes(const surfaceScalarField& alphaf);

        //- Correct primatives based on internal energy/enthalpy
        virtual void correctThermo();

        //- Encode primative variables
        virtual void encode();

        //- Decode conservative variables
        virtual void decode();

        //- Set old time values to current
        virtual void store();

        //- Set the number of steps used in ode solution to hyperbolic terms
        virtual void setNSteps
        (
            const boolList& storeCoeffs,
            const boolList& storeFcoeffs
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
