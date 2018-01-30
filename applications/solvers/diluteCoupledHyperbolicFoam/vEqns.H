forAll(nodes, nodei)
{
    const volVectorNode& node = nodes[nodei];
    const volScalarField& alphapi = node.primaryWeight();
    volVectorField& Upi = nodes[nodei].primaryAbscissa();


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
}
forAll(moments[0], celli)
{
    quadrature.updateLocalMoments(celli);
}