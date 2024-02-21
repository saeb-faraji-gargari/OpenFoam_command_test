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

#include "wallShearStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields_play_test.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    //#include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        //fvModels.correct();
		
      
	  // wallShearStress is tu/rho   unit is [0 2 -2 0 0 0 0]

      volScalarField wallShearStressMag = mag(wallShearStress);

      //const volVectorField& u_star_1 = mesh().lookupObject<volVectorField>("u_star");

      // Cast u_star to a non-const reference
      //volVectorField& nonConstU_star = const_cast<volVectorField&>(u_star_1);

      //volVectorField u_star_new = std::move(u_star);


  
	  //int chert_pert;
	  //std::cin >> chert_pert;

    // Create an instance of the wallShearStress class

    /*Foam::dictionary dict_shearStress
{
    {
        "functionObjectLibs", Foam::List<Foam::word>("libfieldFunctionObjects.so"),
        "type", "wallShearStress",
        "writeInterval", 2,
        "patches", Foam::List<Foam::word>("upperWall")
    }
};*/


    //Foam::functionObjects::wallShearStress wallFunction("wallFunction", runTime.userTimeName(), U);   -------------Error


    /*tmp<volVectorField> twallShearStress              -------------Error
    (
        volVectorField::New
        (
            //type(),
            //mesh_,
            mesh,
            dimensionedVector(dimensionSet(0,1,-1,0,0,0,0), Zero)
        )
    );
    */

   //volVectorField::Boundary& wallShearStressBf = wallShearStress.boundaryField();            -------------Error
   //volVectorField::Boundary& wallShearStressBf = wallShearStress.ref().boundaryFieldRef();   -------------Error
   //volVectorField::Boundary& wallShearStressBf = wallShearStress.ref();                      -------------Error


tmp<volVectorField> tmp_volVectorField_1
(new volVectorField
(
IOobject
(
"tmp_volVectorField_1",
mesh.time().timeName(),
mesh,
IOobject::NO_READ,
IOobject::NO_WRITE
),
mesh,
dimensionedVector
(
"zero",
dimensionSet(0,1,-1,0,0,0,0),
//vector::zero
vector(10,20,50)
)
)
);

volVectorField& tmp_volVectorField_1_ref = tmp_volVectorField_1.ref();

u_star = tmp_volVectorField_1.ref();

Info << "tmp_volVectorField_1_ref : " << tmp_volVectorField_1_ref << endl;


volVectorField::Boundary& tmp_volVectorField_1_Bf =
        tmp_volVectorField_1.ref().boundaryFieldRef();



/*patchSet_ =
        mesh.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );
        */

forAll(tmp_volVectorField_1_Bf, patchi)
    {
        //label patchi = iter.key();

        const vectorField& Sfp = mesh.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];

        tmp_volVectorField_1_Bf[patchi] = (-Sfp/magSfp);
        //tmp_volVectorField_1_Bf[patchi] = (-Sfp/magSfp) & tau.boundaryField()[patchi];

        Info << "patchi: " << patchi   <<   "tmp_volVectorField_1_Bf: " << tmp_volVectorField_1_Bf[patchi] << endl;

        double chert_pert_2;

        std::cin >> chert_pert_2;


    }

volVectorField& u_star_ref_ref = u_star_ref;

volVectorField::Boundary& u_star_ref_ref_Bf =
        u_star_ref_ref.boundaryFieldRef();

forAll(u_star_ref_ref_Bf, patchi)
    {
        //label patchi = iter.key();
        for (int i=1; i<10; i++)
        
        {
        //const vectorField& Sfp = mesh.Sf().boundaryField()[patchi][i];
        const Foam::Vector<double>& Sfp = mesh.Sf().boundaryField()[patchi][i];

        //const scalarField& magSfp = mesh.magSf().boundaryField()[patchi][i];

        u_star_ref_ref_Bf[patchi][i] = (-Sfp);
        //tmp_volVectorField_1_Bf[patchi] = (-Sfp/magSfp) & tau.boundaryField()[patchi];

        Info << "patchi: " << patchi  << "i :" << i <<   "u_star_ref_ref_Bf[patchi][i]: " << u_star_ref_ref_Bf[patchi][i] << endl;
        
        //u_star_ref[patchi][i] = (-Sfp);
        Info << "patchi: " << patchi  << "i :" << i <<   "u_star_ref[patchi][i]: " << u_star_ref[patchi][i] << endl;

        // Print the size of u_star_ref
        Info << "Size of u_star_ref: " << u_star_ref.size() << endl;

        // Print the dimensions of u_star_ref
        Info << "Dimensions of u_star_ref: " << u_star_ref.dimensions() << endl;

        // Print the type of u_star_ref
        Info << "Type of u_star_ref: " << typeid(u_star_ref).name() << endl;
        
        }

    }




tmp<volVectorField> tmp_volVectorField_2
(new volVectorField
(
IOobject
(
"tmp_volVectorField_2",
mesh.time().timeName(),
mesh,
IOobject::NO_READ,
IOobject::NO_WRITE
),
mesh,
dimensionedVector
(
"zero",
dimensionSet(0,1,-1,0,0,0,0),
//vector::zero
vector(10,20,50)
)
)
);

volVectorField& tmp_volVectorField_2_ref = tmp_volVectorField_2.ref();

u_star = tmp_volVectorField_2.ref();

Info << "tmp_volVectorField_2_ref : " << tmp_volVectorField_2_ref << endl;


volVectorField::Boundary& tmp_volVectorField_2_Bf =
        tmp_volVectorField_2.ref().boundaryFieldRef();


forAll(tmp_volVectorField_2_Bf, patchi)
    {
        //label patchi = iter.key();

        const vectorField& Sfp = mesh.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];

        tmp_volVectorField_2_Bf[patchi] = wallShearStress.boundaryField()[patchi];
        //tmp_volVectorField_1_Bf[patchi] = (-Sfp/magSfp) & tau.boundaryField()[patchi];

        Info << "patchi: " << patchi   <<   "tmp_volVectorField_2_Bf: " << tmp_volVectorField_2_Bf[patchi] << endl;


        double chert_pert_2;

        std::cin >> chert_pert_2;


    }


const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

label patch_index_ground;

forAll(mesh.boundary(), patchI)
{
    //const polyPatch& currentPatch = mesh.boundary()[patchI];
	
	const word& patchName = mesh.boundary()[patchI].name();

	 if (patchName == "ground")
	 {
		
		 patch_index_ground = patchI;
		
		Info<< "patchName = " << patchName << " patch_index_ground: " << patch_index_ground << nl << endl;
		

		
		// index fac avale patchI ra neshan modahad labate in index face general hast ahr patch local face index ham darad
       const label& faceStart = mesh.boundaryMesh()[patchI].start();
       
	   const label localFaceStart = mesh.boundaryMesh()[patchI].whichFace(faceStart);
	   
		
		Info<< "faceStart = " << faceStart << " localFaceStart: " << localFaceStart << nl << endl;


		int chert_pert;
	    std::cin >> chert_pert;
		

           
	 }
	 
}


// obtain patchID or index certain name of patches
label patchID_ground = mesh.boundaryMesh().findPatchID("ground");

Info << "patchID_ground: "<< patchID_ground << endl;
	 

const Foam::fvPatch& mesh_boundary_patch_ground = mesh.boundary()[patchID_ground];

Info << "mesh_boundary_patch_ground: "<< mesh_boundary_patch_ground.name() << endl;

int chert_pert1;
std::cin >> chert_pert1;

vector u_star_vector_i;

volVectorField& u_star_nonconst = const_cast<volVectorField&>(u_star);
	 
	 // Loop over all faces of boundary patch
  forAll(mesh.boundary()[patchID_ground], facei)
  {
    const label& bCell = boundaryMesh[patchID_ground].faceCells()[facei];    // Boundary cell index
    const label& face = boundaryMesh[patchID_ground].start() + facei;        // Face index
	
	scalar wallShearStress_Mag = mag(wallShearStress.boundaryField()[patchID_ground][facei]);
	
	Info << "bCell: " << bCell << " face: " << face << endl;
	
	Info << "wallShearStress_Mag: " << wallShearStress_Mag << endl;
	
	vector u_star_vector_i = sqr(wallShearStress_Mag) * (wallShearStress.boundaryField()[patchID_ground][facei]/wallShearStress_Mag);
	
	Info << "u_star_vector_i: " << u_star_vector_i << endl; 


	//u_star.boundaryField()[patchID_ground][facei] =  vector::zero;   //u_star_vector_i;
	//U.boundaryField()[patchID_ground][facei] =  vector::zero;   //u_star_vector_i;


	//int chert_pert;
	//std::cin >> chert_pert;
        
    // Do your calculations e.g.
    // U.boundaryField()[patch][facei] = vector::zero;
	
	/*scalar wallShearStress_Mag = mag(wallShearStress[patchID_ground][facei]);

	u_star[patchID_ground][facei] = sqr(wallShearStress_Mag) * (wallShearStress[patchID_ground][facei]/wallShearStress_Mag);
	
	
	Info << "patchID_ground: " << patchID_ground << " facei: " << facei << " wallShearStress_Mag : " << wallShearStress_Mag << " u_star: " << u_star << endl; 
	
	int chert_pert;
	std::cin >> chert_pert;*/
	
  }
  
  forAll(u_star.boundaryField()[patchID_ground], facei)
  {
	
	scalar wallShearStress_Mag = mag(wallShearStress.boundaryField()[patchID_ground][facei]);
	
	vector u_star_vector_i = sqr(wallShearStress_Mag) * (wallShearStress.boundaryField()[patchID_ground][facei]/wallShearStress_Mag);
    

    std::cout << "typeid(u_star.boundaryField()[patchID_ground][facei]).name(): " << typeid(u_star.boundaryField()[patchID_ground][facei]).name() << std::endl;
    std::cout << "typeid(u_star_vector_i).name(): " << typeid(u_star_vector_i).name() << std::endl;
    
    std::cout << "vector(10,20,50): " << typeid(vector(10,20,50)).name() << std::endl;
    
    std::cout << "variable u_star.boundaryField()[patchID_ground][facei]: " << "\n";

    //std::cout << "u_star.boundaryField()[patchID_ground][facei]: " << u_star.boundaryField()[patchID_ground][facei];

    //std::cout << "u_star.boundaryField()[patchID_ground][facei].ref(): " << u_star.boundaryField()[patchID_ground][facei].ref();
            
            if (std::is_reference<decltype(u_star.boundaryField()[patchID_ground][facei])>::value)
            {
                std::cout << "Variable is a reference" << "\n";
            }
            else if (std::is_pointer<decltype(u_star.boundaryField()[patchID_ground][facei])>::value)
            {
                std::cout << "Variable is a pointer" << "\n";
            }
            else
            { 
                std::cout << "Variable is neither a reference nor a pointer" << "\n";
            }

    std::cout << "variable u_star_vector_i: " << "\n";
            
            if (std::is_reference<decltype(u_star_vector_i)>::value)
            {
                std::cout << "Variable is a reference" << "\n";
            }
            else if (std::is_pointer<decltype(u_star_vector_i)>::value)
            {
                std::cout << "Variable is a pointer" << "\n";
            }
            else
            { 
                std::cout << "Variable is neither a reference nor a pointer" << "\n";
            }

	//u_star.boundaryField()[patchID_ground][facei] = vector(10,20,50);   //u_star_vector_i;    //  vector::zero;   //vector(10,20,50);   //u_star_vector_i;
	
    //Foam::Vector<double>& vector_ref = u_star.boundaryField()[patchID_ground][facei];
    //vector_ref = Foam::Vector<double>(10, 20, 50);  // This modifies the vector object directly

    //u_star.boundaryField()[patchID_ground][facei] = wallShearStress.boundaryField()[patchID_ground][facei];

    //const vector u_star_value(1, 2, 3);
    //u_star.boundaryField()[patchID_ground][facei] = u_star_value;


    //u_star.boundaryField()[patchID_ground][facei] = Foam::Vector<double>(10, 20, 50);

    /*u_star.boundaryField()[patchID_ground][facei].x() = 10;
    u_star.boundaryField()[patchID_ground][facei].y() = 20;
    u_star.boundaryField()[patchID_ground][facei].z() = 50;
    */

    //u_star_nonconst.boundaryField()[patchID_ground][facei] = vector(10,20,50);

    //tmp_volVectorField_1_ref.boundaryField()[patchID_ground][facei] = vector(3,3,3);      --------------Error

    // tmp_volVectorField_1.boundaryField()[patchID_ground][facei] = vector(3,3,3); 

	u_star.correctBoundaryConditions();
	
	Info << "u_star.boundaryField()[patchID_ground][facei]: " << u_star.boundaryField()[patchID_ground][facei] << endl;

	Info << "u_star_nonconst.boundaryField()[patchID_ground][facei]: " << u_star_nonconst.boundaryField()[patchID_ground][facei] << endl;
    
    Info << "mp_volVectorField_1.ref().boundaryField()[patchID_ground][facei]: " << tmp_volVectorField_1.ref().boundaryField()[patchID_ground][facei] << endl;


    //Info << "tmp_volVectorField : " << tmp_volVectorField << endl; 
    //Info << "tmp_volScalarField : " << tmp_volScalarField << endl; 


	
	//vector myVector = vector(0, 0, 0); // Create a vector and initialize it with (0, 0, 0)
	//vector& myVectorRef = myVector.ref(); // Create a reference to myVector

	//u_star.boundaryField()[patchID_ground][facei].write(vector::zero);   // myVector;   //vector::zero;   //vector(10,20,50);   //u_star_vector_i;
	//u_star.boundaryField()[patchID_ground][facei] (vector::zero);
	
  
  }
  
  
  //u_star.write();
  
    //Foam::fvPatch& patch_u_star = u_star.boundaryField()[patchID_ground];

    // This now loops over the faces of the corresponding boundary
    /*forAll (u_star.boundaryField()[patchID_ground], faceI)
    {
        u_star.boundaryField()[patchID_ground][faceI].x() = 0;
		u_star.boundaryField()[patchID_ground][faceI].y() = 10;
		u_star.boundaryField()[patchID_ground][faceI].z() = 800;

    }
	*/

	 
	
		
		
          /*const faceList& faces = mesh.faces();
	 

            forAll(faces, faceI)
            {
              const label& patch = boundaryMesh.whichPatch(faceI);       // Boundary patch index
              const label& facei = boundaryMesh[patch].whichFace(faceI); // Local boundary face index
           
		   
		   scalar wallShearStress_Mag = mag(wallShearStress[patch][facei]);

			u_star[patch][facei] = sqr(wallShearStress_Mag) * (wallShearStress[patch][facei]/wallShearStress_Mag);
			
			int chert_pert;
	        std::cin >> chert_pert;
			
	  
            }
			*/
			
     //}
	

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
		

        //u_star_ref = tmp_volVectorField_1.ref();

		
		// Write u_star to file
        u_star.write();

        tmp_volVectorField_1.ref().write();

        u_star_ref.write();

        tmp_volVectorField_2.ref().write();

        //u_star_ref_ref_Bf.write();

        //tmp_volVectorField_1.write();

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
