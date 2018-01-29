forAll(nodes, nodei)
{
    const volVectorNode& node = nodes[nodei];
    const volScalarField& alphapi = node.primaryWeight();
    volVectorField Upi
    (
        IOobject
        (
            "U.particles",
            runTime.timeName(),
            mesh
        ),
        node.primaryAbscissa(),
        U.boundaryField().types()
    );

    fvVectorMatrix UpiEqn
    (
        fvm::ddt(alphapi, Upi)
      - alphapi*fvc::ddt(Upi)
     ==
        g*alphapi*rhop/rhof

      + Kd*U
      - fvm::Sp(Kd, Upi)
    );

    UpiEqn.relax();
    UpiEqn.solve();

    nodes[nodei].primaryAbscissa() == Upi;
}
forAll(moments[0], celli)
{
    quadrature.updateLocalMoments(celli);
}