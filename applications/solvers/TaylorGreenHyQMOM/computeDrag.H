PtrList<volScalarField> Kds(nodes.size());

forAll(nodes, nodei)
{
    const volVectorNode& node = nodes[nodei];
        const volScalarField& alphapi = node.primaryWeight();
    const volVectorField& Upi = node.primaryAbscissa();

    volScalarField Re(mag(Upi - U)*dp/nu);
    volScalarField CdRe
    (
        24.0*(1.0 + 0.15*pow(Re, 0.687))
    );

    Kds.set
    (
        nodei,
        new volScalarField
        (
            max(alphapi, residualAlpha)*rhof/rhop*0.75*CdRe*nu/sqr(dp)
        )
    );
}