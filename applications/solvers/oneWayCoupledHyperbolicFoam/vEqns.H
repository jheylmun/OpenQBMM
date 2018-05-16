
forAll(nodes, nodei)
{
    Ups[nodei] = nodes[nodei].velocityAbscissae();

    fvVectorMatrix UpEqn
    (
        fvm::ddt(Ups[nodei])
      - fvc::ddt(Ups[nodei])
     ==
        1.0/St*U
      - fvm::Sp(1.0/St, Ups[nodei])
    );

    UpEqn.relax();
    UpEqn.solve();

    nodes[nodei].setVelocityAbscissae(Ups[nodei]);
}

//- Update moments to include drag source
quadrature.updateMoments();