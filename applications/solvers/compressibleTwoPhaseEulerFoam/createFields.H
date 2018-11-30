#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating twoPhaseSystem\n" << endl;
twoPhaseSystem fluid(mesh, g);
volVectorField& U = fluid.U();
volScalarField& p = fluid.p();

phaseModel& phase1 = fluid.phase1();
phaseModel& phase2 = fluid.phase2();

volScalarField& alpha1 = phase1;
volScalarField& alpha2 = phase2;

volVectorField& U1 = phase1.U();
const surfaceScalarField& phi1 = phase1.phi();
volVectorField& U2 = phase2.U();
const surfaceScalarField& phi2 = phase2.phi();

#include "gh.H"