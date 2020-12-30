#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"

Info<< "Reading thermophysical properties\n" << endl;

compressibleSystem fluid(mesh);

fluidThermo& thermo = fluid.thermo();

const volScalarField& rho = fluid.rho();

const volVectorField& U = fluid.U();
const surfaceScalarField& phi = fluid.phi();

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        fluid.massFlux(),
        fluid.thermo()
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermophysicalTransportModel::New(turbulence(), thermo)
);

const volScalarField& muEff(fluid.thermo().mu());
bool inviscid = (max(muEff).value() < small);

errorEstimator error(mesh);

