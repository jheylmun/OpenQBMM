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
    Foam::vdfPhaseModel

SourceFiles
    vdfPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef vdfPhaseModel_H
#define vdfPhaseModel_H

#include "phaseModel.H"
#include "dragModel.H"
#include "granularViscosityModel.H"
#include "radialModel.H"
#include "granularPressureModel.H"
#include "frictionalStressModel.H"
#include "populationBalanceModel.H"
#include "quadratureApproximations.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class vdfPhaseModel Declaration
\*---------------------------------------------------------------------------*/

class vdfPhaseModel
:
    public phaseModel
{
protected:
    // Protected data

    IOdictionary populationBalanceProperties_;

        // Sub-models

            //- Run-time selected viscosity model
            autoPtr<kineticTheoryModels::viscosityModel>
                viscosityModel_;

                //- Run-time selected radial distribution model
            autoPtr<kineticTheoryModels::radialModel> radialModel_;

            //- Run-time selected granular pressure model
            autoPtr<kineticTheoryModels::granularPressureModel>
                granularPressureModel_;

            //- Run-time selected frictional stress model
            autoPtr<kineticTheoryModels::frictionalStressModel>
                frictionalStressModel_;


        // Kinetic Theory Model coefficients

            //- Coefficient of restitution
            dimensionedScalar e_;

            //- Maximum packing phase-fraction
            dimensionedScalar alphaMax_;

            //- Min value for which the frictional stresses are zero
            dimensionedScalar alphaMinFriction_;

            //- Maximum turbulent viscosity
            dimensionedScalar maxNut_;


        // Kinetic Theory Model Fields

            //- The granular energy/temperature
            volScalarField Theta_;

            //- The granular pressure (collisional and frictional)
            volScalarField Ps_;

            //- Frictional pressure
            volScalarField Pfric_;

            //- The granular bulk viscosity
            volScalarField lambda_;

            //- The granular radial distribution
            volScalarField gs0_;

            //- The effective viscosity
            volScalarField nut_;

            //- The frictional viscosity
            volScalarField nuFric_;

            //- The speed of sound
            volScalarField c_;

            //- Total momentum
            volVectorField alphaRhoU_;

            //- List of node momentums
            PtrList<volVectorField> alphaRhoUs_;

            //- Total phasic energy
            volScalarField alphaRhoE_;


        //- Conservative fluxes

            //- volumetric flux
            surfaceScalarField phi_;

            //- Mass flux
            surfaceScalarField massFlux_;

            //- Total energy flux
            surfaceScalarField energyFlux_;


    autoPtr<populationBalanceModel> populationBalance_;

    velocityQuadratureApproximation& quadrature_;

        //- Fields stored for partial ode steps

            //- Changes in mass
            PtrList<PtrList<volScalarField>> moments_;
            PtrList<PtrList<volScalarField>> deltaMoments_;

            PtrList<volScalarField> alphaRhoEs_;
            PtrList<volScalarField> deltaAlphaRhoEs_;


public:

    // Constructors
        vdfPhaseModel
        (
            const twoPhaseSystem& fluid,
            const dictionary& dict,
            const word& phaseName
        );


    //- Destructor
    virtual ~vdfPhaseModel();


    // Access

        virtual label nNodes() const
        {
            return quadrature_.nodes().size();
        }

        //- Return that the phase uses a slave pressure
        virtual bool slavePressure() const
        {
            return false;
        }

        //- Is the phase granular
        virtual bool granular() const
        {
            return true;
        }

        //- Const access to the phase pseudo thermal energy
        virtual tmp<volScalarField> Theta() const
        {
            return Theta_;
        }

        //- Non-const access to the phase pseudo thermal energy
        virtual volScalarField& ThetaRef()
        {
            return Theta_;
        }

        //- Return the particle pressure
        virtual const volScalarField& p() const
        {
            return Ps_;
        }

        //- Non-const access to the phase pressure
        virtual volScalarField& pRef()
        {
            return Ps_;
        }

         //- Const access to the volume fraction
        virtual const volScalarField& alphas(const label i = -1) const
        {
            return
            (
                i == -1
              ? refCast<const volScalarField>(*this)
              : quadrature_.nodes()[i].primaryWeight()
            );
        }

        //- Non-const access to the volume fraction
        virtual volScalarField& alphasRef(const label i = -1)
        {
            return
            (
                i == -1
              ? (*this)
              : quadrature_.nodes()[i].primaryWeight()
            );
        }

        //- Const access to the phase velocity
        virtual const volVectorField& U(const label i = -1) const
        {
            return (i == -1 ? U_ : quadrature_.nodes()[i].primaryAbscissa());
        }

        //- Non-const access to the phase velocity
        virtual volVectorField& URef(const label i = -1)
        {
            return (i == -1 ? U_ : quadrature_.nodes()[i].primaryAbscissa());
        }

        //- Speed of sound in the phase
        virtual tmp<volScalarField> c() const;


        //- Return the phase volumentric flux
        virtual tmp<surfaceScalarField> phi(const label i = -1) const
        {
            return
            (
                i == -1
              ? tmp<surfaceScalarField>(phi_)
              : fvc::flux(quadrature_.nodes()[i].primaryAbscissa())
            );
        }

        //- Return the phase volumentric flux
        virtual tmp<surfaceScalarField> alphaPhi(const label i = -1) const
        {
            return
            (
                i == -1
              ? phi_*fvc::interpolate(*this)
              : fvc::flux
                (
                    quadrature_.nodes()[i].primaryWeight()
                   *quadrature_.nodes()[i].primaryAbscissa()
                )
            );
        }

        //- Const access to phasic density
        virtual tmp<volScalarField> alphaRho(const label i = -1) const
        {
            return
            (
                i == -1
              ? Foam::max((*this), residualAlpha_)*rho_
              : quadrature_.nodes()[i].primaryWeight()*rho_
            );
        }

        //- Non-const access to phasic density
        virtual volScalarField& alphaRhoRef(const label i = -1)
        {
            NotImplemented;
            return *this;
        }

        //- Const access to phasic momentum
        virtual tmp<volVectorField> alphaRhoU(const label i = -1) const
        {
            return
            (
                i == -1
              ? (*this)*rho_*U_
              : quadrature_.nodes()[i].primaryWeight()*rho_
               *quadrature_.nodes()[i].primaryAbscissa()
            );
        }

        //- Non-const access to phasic momentum
        virtual volVectorField& alphaRhoURef(const label i = -1)
        {
            return
            (
                i == -1
              ? alphaRhoU_
              : alphaRhoUs_[i]
            );
        }

        //- Const access to total phasic energy
        virtual tmp<volScalarField> alphaRhoE(const label i = -1) const
        {
            return
            (
                i == -1
              ? (*this)*rho_*E_
              : quadrature_.nodes()[i].primaryWeight()*rho_*E_
            );
        }

        //- Non-const access to total phasic energy
        virtual volScalarField& alphaRhoERef(const label i = -1)
        {
            if (i != -1)
            {
                NotImplemented;
            }
            return alphaRhoE_;
        }

        //- Const access to phase internal energy
        virtual tmp<volScalarField> alphaRhoPTE(const label i = -1) const
        {
            return (*this)*rho_*Theta_*3.0/2.0;
        }

        //- Non-const access to phase internal energy
        virtual volScalarField& alphaRhoPTERef(const label i = -1)
        {
            NotImplemented;
            return alphaRhoE_;
        }

        //- Const access to phasic mass flux
        virtual tmp<surfaceScalarField> massFlux() const
        {
            return massFlux_;
        }

        //- Const access to phasic momentum flux
        virtual tmp<surfaceVectorField> momentumFlux() const
        {
            return phi_*fvc::interpolate(U_*(*this)*rho_);
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
            return fvc::interpolate(*this);
        }

        //- Return the gradient of velocity
        virtual tmp<volTensorField> gradU(const label i = -1) const
        {
            return fvc::grad(U_);
        }

        //- Return the gradient of pressure
        virtual tmp<volVectorField> gradp() const
        {
            return fvc::grad(Ps_ + Pfric_);
        }

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


    // Member Functions

        //- Return granular dissipation coefficient
        virtual tmp<volScalarField> dissipationCoeff() const;

        //- Return granular production coefficient
        virtual tmp<volScalarField> productionCoeff() const;

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

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);

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
