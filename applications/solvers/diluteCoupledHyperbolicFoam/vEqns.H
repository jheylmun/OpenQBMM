forAll(nodes, nodei)
{
    const volNode& node = nodes[nodei];
    Ups[nodei] = node.velocityAbscissae();
    volScalarField alphaRho(node.primaryWeight()*rhop);


    fvVectorMatrix UpEqn
    (
        alphaRho*fvm::ddt(Ups[nodei])
      - alphaRho*fvc::ddt(Ups[nodei])
     ==
        g*alphaRho

      + Kds[nodei]*U
      - fvm::Sp(Kds[nodei], Ups[nodei])
    );

    UpEqn.relax();
    UpEqn.solve();

    nodes[nodei].setVelocityAbscissae(Ups[nodei]);
}

// Update moments to include drag and gravity
quadrature.updateMoments();
