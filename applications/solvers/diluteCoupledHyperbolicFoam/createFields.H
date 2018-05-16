#include "readGravitationalAcceleration.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
const volScalarField& nuf = laminarTransport.nu();

// Particle fields

Info<< "Reading phaseProperties\n" << endl;
IOdictionary phaseProperties
(
    IOobject
    (
        "phaseProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
word particleName(phaseProperties.lookup("dilutePhase"));
scalar residualAlpha(phaseProperties.lookupType<scalar>("residualAlpha"));
dimensionedScalar rhof(phaseProperties.lookup("rhof"));
dimensionedScalar rhop(phaseProperties.lookup("rhop"));


// Field is created so it can be looked up in collision kernels
volScalarField dp
(
    IOobject
    (
        IOobject::groupName("d", particleName),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phaseProperties.lookup("d"))
);

Info<< "Reading populationBalanceProperties\n" << endl;
IOdictionary populationBalanceProperties
(
    IOobject
    (
        "populationBalanceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

autoPtr<populationBalanceModel> populationBalance
(
    populationBalanceModel::New
    (
        particleName, populationBalanceProperties, phi
    )
);
velocityQuadratureApproximation& quadrature =
    mesh.lookupObjectRef<velocityQuadratureApproximation>
    (
        IOobject::groupName
        (
            "quadratureProperties",
            particleName
        )
    );

const labelListList& momentOrders = quadrature.momentOrders();
const labelListList& nodeIndexes = quadrature.nodeIndexes();

volMomentFieldSet& moments = quadrature.moments();
mappedPtrList<volNode>& nodes = quadrature.nodes();
mappedPtrList<volVectorField> Ups(nodes.size(), nodeIndexes);
const label nDims = momentOrders[0].size();

const volScalarField& alpha = moments(labelList(nDims, 0));

volVectorField Up
(
    IOobject
    (
        IOobject::groupName("U", particleName),
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Thetap
(
    IOobject
    (
        IOobject::groupName("Theta", particleName),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", sqr(dimVelocity), 0.0)
);

forAll(Ups, nodei)
{
    const labelList& nodeIndex = nodeIndexes[nodei];
    Ups.set
    (
        nodeIndex,
        new volVectorField
        (
            IOobject
            (
                IOobject::groupName
                (
                    "U" + mappedPtrList<scalar>::listToWord(nodeIndex),
                    particleName
                ),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("Up", dimVelocity, Zero),
            Up.boundaryField().types()
        )
    );
}

#include "computeParticleFields.H"

#include "createMRF.H"
#include "createFvOptions.H"
