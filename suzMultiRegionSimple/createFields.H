Info<< "Reading tansportProperties\n" << endl;
IOdictionary transportProperties
(
     IOobject
     (
    "transportProperties",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
     )
 );
  
Info << "Reading air density rho" << endl;
dimensionedScalar rho
(
    "rho",
    dimensionSet(1, -3, 0, 0, 0, 0, 0),
    transportProperties.lookup("rho")
);
Info << "    rho = " << rho << endl;

Info << "Reading kinematic viscosity nu" << endl;
dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);
Info << "    nu = " << nu << endl;

Info << "Reading relative permitivity, epsilon0" << endl;
dimensionedScalar epsilon0
(
    "epsilon0",
    transportProperties.lookup("epsilon0")
);
Info << "    epsilon0 = " << epsilon0 << endl;

Info << "Reading Debye length lambda" << endl;
dimensionedScalar lambda
(
    "lambda",
    dimensionSet(0, 1, 0, 0, 0, 0, 0),
    transportProperties.lookup("lambda")
);
Info << "    lambda = " << lambda << endl;

Info << "Reading Max electric potential epMax" << endl;
dimensionedScalar epMax
(
     "epMax",
     dimensionSet(0, 0, 0, 0, 0, 0, 0),
     transportProperties.lookup("epMax")
);
Info << "    epMax = " << epMax << endl;

Info << "Reading Max charge density rhoCMax" << endl;
dimensionedScalar rhoCMax
(
     "rhoCMax",
     dimensionSet(0, 0, 0, 0, 0, 0, 0),
     transportProperties.lookup("rhoCMax")
 );
Info << "    rhoCMax = " << rhoCMax << endl;

Info<< "Reading field electric potential ep0\n" << endl;
volScalarField epA0
(
    IOobject
    (
        "ep0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Calculating field electric potential ep1\n" << endl;
volScalarField epA1
(
    IOobject
    (
        "ep1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Calculating field epTot\n" << endl;
volScalarField epATot
(
    IOobject
    (
        "epTot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(1,2,-3,0,0,-1,0),scalar(0))
);

Info<< "Calculating electric field E0\n" << endl;
volVectorField EA0
(
    IOobject
    (
        "E0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(epA0)
);

Info<< "Calculating electric field Echarge\n" << endl;
volVectorField EAcharge
(
    IOobject
    (
        "Echarge",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(epA1)
);

Info<< "Calculating field electric field Etot\n" << endl;
volVectorField EATot
(
    IOobject
    (
        "ETot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    //mesh
    -fvc::grad(epATot)
);


Info<< "Reading field relative permitivity epsR\n" << endl;
volScalarField epsR
(
     IOobject
     (
        "epsR",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
     ),
     mesh
);

Info<< "Reading field charge density rhoC\n" << endl;
volScalarField rhoC
(
      IOobject
      (
        "rhoC",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
      mesh
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
