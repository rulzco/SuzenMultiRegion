Info<< "Reading field electric potential ep0\n" << endl;
volScalarField epS0
(
    IOobject
    (
        "ep0",
        runTime.timeName(),
        meshSolid,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSolid
);

Info<< "Reading field electric potential epS1\n" << endl;
volScalarField epS1
(
    IOobject
    (
        "ep1",
        runTime.timeName(),
        meshSolid,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshSolid
);

Info<< "Reading field epTot\n" << endl;
volScalarField epSTot
(
    IOobject
    (
        "epTot",
        runTime.timeName(),
        meshSolid,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshSolid,
    dimensionedScalar("temp",dimensionSet(1,2,-3,0,0,-1,0),scalar(0))
);

Info<< "Calculating electric field E0\n" << endl;
volVectorField ES0
(
    IOobject
    (
        "E0",
        runTime.timeName(),
        meshSolid,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(epS0)
);

Info<< "Calculating electric field Echarge\n" << endl;
volVectorField EScharge
(
    IOobject
    (
        "Echarge",
        runTime.timeName(),
        meshSolid,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(epS1)
);


Info<< "Calculating field electric field Etot\n" << endl;
volVectorField ESTot
(
    IOobject
    (
        "ETot",
        runTime.timeName(),
        meshSolid,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(epSTot)
);

Info<< "Reading field relative permitivity epsD\n" << endl;
volScalarField epsD
(
     IOobject
     (
        "epsD",
        runTime.timeName(),
        meshSolid,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
     ),
     meshSolid
);