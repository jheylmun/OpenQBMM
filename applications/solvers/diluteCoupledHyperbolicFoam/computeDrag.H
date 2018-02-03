PtrList<volScalarField> Kds(nodes.size());
volScalarField Kd
(
    IOobject
    (
        "Kd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", inv(dimTime), 0.0)
);

{
    volScalarField alphaf
    (
        max(scalar(1) - alpha, residualAlpha)
    );

    forAll(nodes, nodei)
    {
        const volVectorNode& node = nodes[nodei];
        const volScalarField& alphapi = node.primaryWeight();
        const volVectorField& Upi = node.primaryAbscissa();

        // WenYu drag
        volScalarField alphafi
        (
            max(scalar(1) - alphapi, residualAlpha)
        );
        volScalarField Re(mag(Upi - U)*dp/nuf);
        volScalarField Res(alphafi*Re);
        volScalarField CdsRes
        (
            neg(Res - 1000)*24*(1.0 + 0.15*pow(Res, 0.687))
          + pos0(Res - 1000)*0.44*max(Res, 1e-4)
        );

        volScalarField CdRe
        (
            CdsRes
           *pow(alphafi, -3.65)
           *max(1.0 - alpha, residualAlpha)
        );

        volScalarField Kdi
        (
            0.75
           *max(alphapi, residualAlpha)
           *rhof/rhop
           *CdRe
           *nuf
           /sqr(dp)
        );
        Kds.set(nodei, new volScalarField(Kdi));
        Kd += Kdi;
    }
}