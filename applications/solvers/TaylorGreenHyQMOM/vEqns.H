forAll(nodes, nodei)
{
    const volVectorNode& node = nodes[nodei];
    const volScalarField& alphapi = node.primaryWeight();
    volVectorField& Upi = nodes[nodei].primaryAbscissa();


    fvVectorMatrix UpiEqn
    (
        alphapi*fvm::ddt(Upi)
      - alphapi*fvc::ddt(Upi)
     ==
        Kds[nodei]*U
      - fvm::Sp(Kds[nodei], Upi)
    );

    UpiEqn.relax();
    UpiEqn.solve();
}
forAll(moments[0], celli)
{
    quadrature.updateLocalMoments(celli);
}