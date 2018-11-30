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

Class
    Foam::phaseModel

SourceFiles
    phaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "rhoThermo.H"
#include "fvMatrix.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class twoPhaseSystem;
class diameterModel;

template<class Phase>
class PhaseCompressibleTurbulenceModel;


/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField
{
protected:
    // Protected data

        //- Reference to the twoPhaseSystem to which this phase belongs
        const twoPhaseSystem& fluid_;

        //- Name of phase
        word name_;

        // Phase dictionary
        dictionary phaseDict_;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        dimensionedScalar residualAlpha_;

        //- Return the residual phasic density for given phase
        dimensionedScalar residualRho_;

        //- Rho thermo class
        autoPtr<rhoThermo> thermoPtr_;


        // Primative variables

            //- Phase density
            volScalarField& rho_;

            //- Phase velocity
            volVectorField U_;

            //- Phase internal energy
            volScalarField& he_;

            //- Phase total energy
            volScalarField E_;


        //- Gradient of alpha
        volVectorField gradAlpha_;

        //- Diameter model
        autoPtr<diameterModel> dPtr_;

        //- Turbulence model
        autoPtr<PhaseCompressibleTurbulenceModel<phaseModel>> turbulence_;

        //- fields to store for ode time stepping
        labelList storedFieldIndexes_;
        labelList storedDeltaIndexes_;

protected:

    void setTurbulenceModel();

public:

    //- Runtime type information
    ClassName("phaseModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseModel,
            dictionary,
            (
                const twoPhaseSystem& fluid,
                const dictionary& dict,
                const word& phaseName
            ),
            (fluid, dict, phaseName)
        );

    // Constructors
        phaseModel
        (
            const twoPhaseSystem& fluid,
            const dictionary& dict,
            const word& phaseName
        );


    // Selectors

        static autoPtr<phaseModel> New
        (
            const twoPhaseSystem& fluid,
            const dictionary& dict,
            const word& phaseName
        );


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        //- Return the name of this phase
        inline const word& name() const;

        //- Return the twoPhaseSystem to which this phase belongs
        inline const twoPhaseSystem& fluid() const;

        //- Return the other phase in this two-phase system
        const phaseModel& otherPhase() const;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        inline const dimensionedScalar& residualAlpha() const;

        //- Return that the phase uses a slave pressure
        virtual bool slavePressure() const = 0;

        //- Is the phase granular
        virtual bool granular() const = 0;

        //- Return the residual phase-density for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        inline const dimensionedScalar& residualRho() const;

        //- Return the residual alphaRho
        inline dimensionedScalar residualAlphaRho() const;

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        virtual scalar alphaMax() const
        {
            return 1.0;
        }

        //- Return the turbulence model
        const PhaseCompressibleTurbulenceModel<phaseModel>&
            turbulence() const;

        //- Return non-const access to the turbulence model
        PhaseCompressibleTurbulenceModel<phaseModel>&
            turbulence();

        //- Return the thermodynamic model
        inline const rhoThermo& thermo() const;

        //- Return non-const access to the thermodynamic model
        inline rhoThermo& thermo();

        //- Const access to the phase density
        inline const volScalarField& rho() const;

        //- Non-const access to the phase density
        inline volScalarField& rho();

        //- Const access to the phase velocity
        inline const volVectorField& U() const;

        //- Non-const access to the phase velocity
        inline volVectorField& U();

        //- Return the phase volumentric flux
        virtual const surfaceScalarField& phi() const = 0;

        //- Return the phase volumentric flux
        virtual tmp<surfaceScalarField> alphaPhi() const = 0;

        //- Const access to the phase total energy
        inline const volScalarField& E() const;

        //- Non-const access to the phase total energy
        inline volScalarField& E();

        //- Const access to the phase pseudo thermal energy
        virtual const volScalarField& Theta() const = 0;

        //- Non-const access to the phase pseudo thermal energy
        virtual volScalarField& Theta() = 0;

        //- Const access to the phase pressure
        virtual const volScalarField& p() const = 0;

        //- Non-const access to the phase pressure
        virtual volScalarField& p() = 0;

        //- Return the Sauter-mean diameter
        tmp<volScalarField> d() const;

        //- Speed of sound of the phase
        virtual tmp<volScalarField> c() const = 0;

        //- Return the laminar viscosity
        inline tmp<volScalarField> nu() const;

        //- Return the laminar viscosity for patch
        inline tmp<scalarField> nu(const label patchi) const;

        //- Return the laminar dynamic viscosity
        inline tmp<volScalarField> mu() const;

        //- Return the laminar dynamic viscosity for patch
        inline tmp<scalarField> mu(const label patchi) const;

        //- Return the thermal conductivity on a patch
        inline tmp<scalarField> kappa(const label patchi) const;

        //- Return the thermal conductivity
        inline tmp<volScalarField> kappa() const;

        //- Thermal diffusivity for energy of mixture [kg/m/s]
        inline tmp<volScalarField> alphahe() const;

        //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        inline tmp<scalarField> alphahe(const label patchi) const;

        //- Return the laminar thermal conductivity
        inline tmp<volScalarField> kappaEff(const volScalarField& alphat) const;

        //- Return the laminar thermal conductivity on a patch
        inline tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;

        //- Return the laminar thermal diffusivity for enthalpy
        inline tmp<volScalarField> alpha() const;

        //- Return the laminar thermal diffusivity for enthalpy on a patch
        inline tmp<scalarField> alpha(const label patchi) const;

        //- Return the effective thermal diffusivity for enthalpy
        inline tmp<volScalarField> alphaEff
        (
            const volScalarField& alphat
        ) const;

        //- Return the effective thermal diffusivity for enthalpy on a patch
        inline tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;

        //- Return the specific heat capacity
        inline tmp<volScalarField> Cp() const;

        //- Return the specific heat capacity
        inline tmp<volScalarField> Cv() const;

        //- Return the specific heat ratio
        inline tmp<volScalarField> gamma() const;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const;

        //- Return the effective viscosity for patch
        virtual tmp<scalarField> nuEff(const label patchi) const;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<volScalarField> pPrime() const;

        //- Return the face-phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<surfaceScalarField> pPrimef() const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Const access to phasic density
        virtual const volScalarField& alphaRho() const = 0;

        //- Non-const access to phasic density
        virtual volScalarField& alphaRho() = 0;

        //- Const access to phasic momentum
        virtual const volVectorField& alphaRhoU() const = 0;

        //- Non-const access to phasic momentum
        virtual volVectorField& alphaRhoU() = 0;

        //- Const access to total phasic energy
        virtual const volScalarField& alphaRhoE() const = 0;

        //- Non-const access to total phasic energy
        virtual volScalarField& alphaRhoE() = 0;

        //- Const access to phase internal energy
        virtual const volScalarField& alphaRhoPTE() const = 0;

        //- Non-const access to phase internal energy
        virtual volScalarField& alphaRhoPTE() = 0;

        //- Const access to phasic mass flux
        virtual const surfaceScalarField& massFlux() const = 0;

        //- Non-const access to phasic mass flux
        virtual surfaceScalarField& massFlux() = 0;

        //- Const access to phasic momentum flux
        virtual const surfaceVectorField& momentumFlux() const = 0;

        //- Non-const access to phasic momentum flux
        virtual surfaceVectorField& momentumFlux() = 0;

        //- Const access to total phasic energy flux
        virtual const surfaceScalarField& energyFlux() const = 0;

        //- Non-const access to total phasic energy flux
        virtual surfaceScalarField& energyFlux() = 0;

        //- Const access to phase internal energy flux
        virtual const surfaceScalarField& PTEFlux() const = 0;

        //- Non-const access to phase internal energy flux
        virtual surfaceScalarField& PTEFlux() = 0;

        //- Return the gradient of pressure
        virtual const surfaceScalarField& alphaf() const = 0;

        //- Return the gradient of volume fraction
        virtual const volVectorField& gradAlpha() const
        {
            return gradAlpha_;
        }

        //- Return non-const access to the gradient of volume fraction
        virtual volVectorField& gradAlpha()
        {
            return gradAlpha_;
        }

        //- Return the gradient of velocity
        virtual tmp<volTensorField> gradU() const = 0;

        //- Return the gradient of pressure
        virtual tmp<volVectorField> gradp() const = 0;

        //- Ensure that the flux at inflow/outflow BCs is preserved
        void correctInflowOutflow(surfaceScalarField& alphaPhi) const;

    // Edit

        //- Correct the phase properties
        //  other than the thermodynamics and turbulence
        //  which have special treatment
        void correct();

        //- Set the number of steps used in ode solution to hyperbolic terms
        virtual void setNSteps
        (
            const boolList& storeCoeffs,
            const boolList& storeFcoeffs
        ) = 0;

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
        ) = 0;

        //- Return granular dissipation coefficient
        virtual tmp<volScalarField> dissipationCoeff() const;

        //- Return granular production coefficient
        virtual tmp<volScalarField> productionCoeff() const;

        //- Decode conservative variables
        virtual void updateFluxes() = 0;

        //- Decode conservative variables
        virtual void updateFluxes(const surfaceScalarField& alphaf) = 0;

        //- Correct primatives based on internal energy/enthalpy
        virtual void correctThermo() = 0;

        //- Encode primative variables
        virtual void encode() = 0;

        //- Decode conservative variables
        virtual void decode() = 0;

        //- Set old time values to current
        virtual void store() = 0;


    // Read
        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);

        //- Dummy Read for transportModel
        virtual bool read()
        {
            return true;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
