Foam::fvMesh meshSolid
(
    IOobject
    (
        "dielectric",
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);
