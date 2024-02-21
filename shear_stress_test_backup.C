/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    scalarTransportFoam

Description
    Solves the steady or transient transport equation for a passive scalar.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    //#include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        //fvModels.correct();
		
      
	  // wallShearStress is tu/rho   unit is [0 2 -2 0 0 0 0]
	  
	

      volScalarField wallShearStressMag = mag(wallShearStress);

      //volVectorField u_star_direction = wallShearStress / wallShearStressMag;

      //u_star = sqr(wallShearStressMag) * u_star_direction;
	  
	  //const fvMesh& mesh = runTime.constant().mesh();
	  
	  //const label mesh_size = mesh().cells().size();

	  /*
	  //forAll(wallShearStressMag, cellI)
	  forAll(mesh().cells(),cellI)

    {
		
	Info<< "cellI = " << cellI << "  wallShearStressMag[cellI] = " << wallShearStressMag[cellI] <<  nl << endl;
	  
  
	  //int chert_pert;
	  //std::cin >> chert_pert;
	  
	  
	  
        if (wallShearStressMag[cellI] > 0.0)
			
        {
         
         //vector u_star_direction = wallShearStress[cellI] / wallShearStressMag[cellI];

         u_star[cellI] = sqr(wallShearStressMag[cellI]) * (wallShearStress[cellI] / wallShearStressMag[cellI]);
		 
		 // Print the cell index and u_star vector components
		 
		Info<< "u_star = " <<  nl << endl;

        std::cout << "Cell " << cellI << " u_star: ";
        for (int component = 0; component < u_star.size(); ++component)
        {
            std::cout << u_star[cellI][component] << " ";
        }
        std::cout << std::endl;
		 
		 //int chert_pert;
		 //std::cin >> chert_pert;
		 
        }
		
		/*else
		{
		
		u_star[cellI] = vector::zero;
		
		}*/
		
		//volVectorField u_star_1 = u_star;
		
		
/*IOobject wallShearStressMagOutput(
    "wallShearStressMag",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE  // No need to write during the simulation, as we'll write it once outside the loop
);

IOobject u_star_1_output(
    "u_star_1",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
);*

		
		
//wallShearStressMag.write();

//u_star_1.write();


    }
	
	int chert_pert;
	std::cin >> chert_pert;
	
	 
      Info<< "u_star is calculated = " << u_star() << nl << endl; */
	  
	  
	  //#include "polyMesh.H"


//const polyMesh& mesh = runTime.io().mesh();
//const polyMesh& mesh = runTime.io().mesh();
//const polyMesh& mesh = runTime.mesh();
//const polyMesh& mesh = runTime.constant().lookupObject<polyMesh>();

const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

forAll(mesh.boundary(), patchI)
{
    //const polyPatch& currentPatch = mesh.boundary()[patchI];
	
	const word& patchName = mesh.boundary()[patchI].name();
	
	
	
	if (patchName == "ground")
	{
     const faceList& faces = mesh.faces();
	 

            forAll(faces, faceI)
            {
              const label& patch = boundaryMesh.whichPatch(faceI);       // Boundary patch index
              const label& facei = boundaryMesh[patch].whichFace(faceI); // Local boundary face index
           
		   
		   scalar wallShearStress_Mag = mag(wallShearStress[patch][facei]);

			u_star[patch][facei] = sqr(wallShearStress_Mag) * (wallShearStress[patch][facei]/wallShearStress_Mag);
            }
			
     }
	
	
		
	//const labelList& faces = currentPatch.faceLabels();
    //const label& bCell = boundaryMesh[patchI].faceCells()[facei];    

    
        //label faceLabel = faces[faceI];
        // Access and process boundary face with label faceLabel
		
		
		//u_star[bCellI] = sqr(wallShearStressMag[bCellI]) * (wallShearStress[bCellI] / wallShearStressMag[bCellI]);
    
	}

     
	 /*fvVectorField u_star_2 
	 (
     fvc::u_star
	 );*/

        /*while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + fvm::div(phi, T)
              - fvm::laplacian(DT, T)
             ==
                fvModels.source(T)
            );

            TEqn.relax();
            fvConstraints.constrain(TEqn);
            TEqn.solve();
            fvConstraints.constrain(T);
        }*/
		
		
		

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
