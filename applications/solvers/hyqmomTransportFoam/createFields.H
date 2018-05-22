//- Empty phi field to pass into poplation balance model
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("zero", dimVolume/dimTime, 0.0)
);

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
word particleName(phaseProperties.lookup("particleName"));
scalar residualAlpha(phaseProperties.lookupType<scalar>("residualAlpha"));

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
    dimensionedScalar("d", dimLength, phaseProperties)
);

volScalarField alphap
(
    IOobject
    (
        IOobject::groupName("alpha", particleName),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Zero", dimless, 0.0)
);

volVectorField Up
(
    IOobject
    (
        IOobject::groupName("U", particleName),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Zero", dimVelocity, Zero)
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
        "particles", populationBalanceProperties, phi
    )
);

const quadratureApproximation& quadrature =
    mesh.lookupObject<quadratureApproximation>
    (
        IOobject::groupName
        (
            "quadratureProperties",
            particleName
        )
    );
const labelListList& momentOrders = quadrature.momentOrders();



