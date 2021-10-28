/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    suzMultiSimple

Description
    Steady-state multiregion solver for the DBD Suzen model coupled with an steady-state solver for 
	incompressible, turbulent flow, using the SIMPLE algorithm.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    
     #include "setRootCaseLists.H"
     #include "createTime.H"
     #include "createMesh.H"
     #include "createDielectricMesh.H"
     #include "createControl.H"
     #include "createFields.H"
     #include "createFieldsSolid.H"
     #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    //Info<< "\nCalculating Electric Potential\n" << endl;
    //solve ( fvm::laplacian(epsR,ElPot) );
    
    scalar          epA0InitialResidual = 1.0;
    scalar          epS0InitialResidual = 1.0;
    int             Iter = 0;
    
    //scalar eqnInitialResidual = 1;
    //label  eqnNoIterations = 0;
    //lduMatrix::solverPerformance solverPerf;

    //solverPerf = Eqn.solve();
    //eqnNoIterations = solverPerf.nIterations();
    //eqnInitialResidual = solverPerf.initialResidual();

    Info<< "\nCalculating electric potential ep0\n" << endl;
    
    while ((epA0InitialResidual > 1.0e-5 || epS0InitialResidual > 1.0e-5) && (Iter < 1000))
    {
        Iter ++;
        Info<< "Iteration = " << Iter << nl << endl;

        //solve electric potential field for air-region
        Foam::solverPerformance sp_epA0 = solve( fvm::laplacian(epsR,epA0) );
        
        epA0InitialResidual = sp_epA0.initialResidual(); //retrive initial residual for epA0 
        
        //solve electric potential field for solid-region
        Foam::solverPerformance sp_epS0 = solve( fvm::laplacian(epsD,epS0) );

        epS0InitialResidual = sp_epS0.initialResidual(); //retrive initial residual for epS0 
        
    }
        
    if (Iter == 1000000)
    {
        cerr << "Electric potential solution diverged!" << endl;
        exit(EXIT_FAILURE);
    }

    Info<< "\nCalculating charge density\n" << endl;
    solve ( fvm::laplacian(epsR,rhoC) == fvm::Sp(1. / (lambda * lambda),rhoC) ); //poisson eq. for the charge density
    
    scalar          epA1InitialResidual = 1.0;
    scalar          epS1InitialResidual = 1.0;
    int             Iter1 = 0;
    
    Info<< "\nCalculating electric potential ep1\n" << endl;
    
    while ((epA1InitialResidual > 1.0e-5 || epS1InitialResidual > 1.0e-5) && (Iter1 < 1000))
    {
        Iter1 ++;
        Info<< "Iteration = " << Iter << nl << endl;

        //solve electric potential field for air-region
        Foam::solverPerformance sp_epA1 = solve( fvm::laplacian(epsR,epA1) + rhoC/epsilon0 );
                
        epA1InitialResidual = sp_epA1.initialResidual(); //retrive initial residual for epA1 

        //solve electric potential field for solid-region
        Foam::solverPerformance sp_epS1 = solve( fvm::laplacian(epsD,epS1) );

        epS1InitialResidual = sp_epS1.initialResidual(); //retrive initial residual for epS1 
    }
    
    if (Iter1 == 1000000)
    {
        cerr << "Electric potential solution diverged!" << endl;
        exit(EXIT_FAILURE);
    }

    //total electric potential
    epATot = epA0 + epA1;
    epSTot = epS0 + epS1;
	
	EA0 = - fvc::grad(epA0);
    EAcharge = - fvc::grad(epA1);
    EATot = - fvc::grad(epATot);

    ES0 = - fvc::grad(epS0);
    EScharge = - fvc::grad(epS1);
    ESTot = - fvc::grad(epSTot);
	
    
    Info<< "\nCalculating time-independent force bForce\n" << endl;
    volVectorField bForce
    (
        IOobject
        (
            "bForce",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rhoCMax * epMax * rhoC * EATot / rho
    );
      
    runTime++;
    runTime.write();
    //ElPot.write();
    //rhoC.write();
    //Efield.write();
    //bForce.write();

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
