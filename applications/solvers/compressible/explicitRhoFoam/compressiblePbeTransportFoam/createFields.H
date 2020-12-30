#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"

Info<< "Reading thermophysical properties\n" << endl;

compressibleSystem fluid(mesh);

fluidThermo& thermo = fluid.thermo();
thermo.validate(args.executable(), "e");

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

Info<< "Reading populationBalanceProperties\n" << endl;

IOdictionary populationBalanceProperties
(
    IOobject
    (
        "populationBalanceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

autoPtr<populationBalanceModel> populationBalance
(
    populationBalanceModel::New
    (
        "populationBalance", populationBalanceProperties, phi
    )
);
