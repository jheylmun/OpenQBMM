Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("0", dimVelocity, Zero)
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimVolume/dimTime, 0.0)
);

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
dimensionedScalar St
(
    "St",
    dimTime,
    phaseProperties
);

const vectorField& points = mesh.C();
const scalar twoPi = 2.0*Foam::constant::mathematical::pi;
forAll (U, celli)
{
    scalar x = points[celli].x();
    scalar y = points[celli].y();
    U[celli] =
        vector
        (
            Foam::sin(twoPi*x)*Foam::cos(twoPi*y),
            -Foam::sin(twoPi*y)*Foam::cos(twoPi*x),
            0.0
        );
}

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
volVectorMomentFieldSet& moments = quadrature.moments();
mappedPtrList<volVectorNode>& nodes = quadrature.nodes();
const label nDims = momentOrders[0].size();

const volScalarField& alpha = moments(labelList(nDims, 0));
Info<<labelList(nDims, 0)<<endl;
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
    dimensionedVector("zero", dimVelocity, Zero)
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

#include "computeParticleFields.H"
