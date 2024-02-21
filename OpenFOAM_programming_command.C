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

    #include "createFields_OpenFOAM_programming_command.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    //#include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

      

      volScalarField wallShearStressMag = mag(wallShearStress);


// 1)  define temporary variable --- tarfi variable temporary 

tmp<volVectorField> u_star
(new volVectorField
(
IOobject
(
"u_star_result",
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
vector::zero
)
)
);

// volVectorField& u_star_ref = u_star.ref();

// 2) Create reference to boundary of u_star  

volVectorField::Boundary& u_star_Bf =
        u_star.ref().boundaryFieldRef();


// /*  BLOCK A

// compute and update u_star_Bf which automatically update u_star too since u_star_Bf is referred to u_star

forAll(u_star_Bf, patchi)    // loop allover patches mesh ID 
    {
        //label patchi = iter.key();

        // 3) mesh face normal vector on a boundary field

        const vectorField& Sfp = mesh.Sf().boundaryField()[patchi];

        // 4) magnitude of mesh face normal vector on a boundary field
        const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];

        Info << "Sfp: " << Sfp << endl;
        Info << "magSfp: " << magSfp << endl;

        //const vectorField& Sfp = mesh.Sf().boundaryField()[patchi];
        //const scalarField& sqr_mag_wallShearStress = sqr(mag(wallShearStress.boundaryField()[patchi])); 
        //Info << "sqr_mag_wallShearStress: " << sqr_mag_wallShearStress << endl;

        //auto wallShearStress_Mag = mag(wallShearStress.boundaryField()[patchi]);
        //const scalarField& wallShearStress_Mag = mag(wallShearStress.boundaryField()[patchi]);  Error in run segmentationfault since if it is defined as reference 
        //it will be deallocated after first time calling noting that mag(wallShearStress.boundaryField()[patchi]) return a temporary object
        

        // 5) define magnitude of vector field
        const scalarField wallShearStress_Mag = mag(wallShearStress.boundaryField()[patchi]);
	    Info << " patchi: " << patchi << endl;
	    Info << "wallShearStress_Mag: " << mag(wallShearStress.boundaryField()[patchi]) << endl;
	
	   //u_star_Bf[patchi] = (-Sfp/magSfp);
    
        u_star_Bf[patchi] = wallShearStress.boundaryField()[patchi];

        // 5) find the id number of a patch with its name
        label patchID_ground = mesh.boundaryMesh().findPatchID("ground");

        Info << "mesh.boundaryMesh().whichPatchID(\"ground\"): " << patchID_ground << endl;


        forAll(u_star_Bf[patchi], facei)
        {
            // 9) find a variable vector 
            Info << "patchi, componentIndex: "
                 << " (" << patchi << "," << facei << ") "
                 << "wallShearStress.boundaryField()[patchi][facei]: "
                 << wallShearStress.boundaryField()[patchi][facei] << endl;
            
            // 10) find a variable x component 
            Info << "patchi, componentIndex: "
                 << " (" << patchi << "," << facei << ") "
                 << "wallShearStress.boundaryField()[patchi][facei].x(): "
                 << wallShearStress.boundaryField()[patchi][facei].x() << endl;
            
            // 11) find a mesh boundary name, size of boundar patches and start number of boundary patches and end  number of boundary patches
            Info << "mesh.boundary()[patchi].name(): " << mesh.boundary()[patchi].name() 
                 << " mesh.boundary()[patchi].Cf().size() :"  << mesh.boundary()[patchi].Cf().size()
                 << " mesh.boundary()[patchi].start(): " << mesh.boundary()[patchi].start() 
                 << " mesh.boundary()[patchi].start() + mesh.boundary()[patchi].size(): "<< mesh.boundary()[patchi].start() + mesh.boundary()[patchi].size() << endl;
            
            
            //Info << " wallShearStress.boundaryField()[patchi]: "<< wallShearStress.boundaryField()[patchi] << endl;
            // in dorost nist : Info << " wallShearStress.boundary()[patchi].start() + wallShearStress.boundary()[patchi].size(): "<< wallShearStress.boundaryField()[patchi].start() + wallShearStress.boundaryField()[patchi].size() << endl;

            // 6) find the id number of a patch that includes the input face id 
            //label patch_id= mesh.boundaryMesh().whichPatch(faceid);  return -1 if faceid is not involved in boundary patches faces
            label patch_name_ground_face_164654 = mesh.boundaryMesh().whichPatch(164654);
            label patch_name_ground = mesh.boundaryMesh().whichPatch(facei);

            Info << "mesh.boundaryMesh().whichPatch(facei): " << patch_name_ground << endl;
            Info << "mesh.boundaryMesh().whichPatch(164654): " << patch_name_ground_face_164654 << endl;


            //Info << "wallShearStress.boundaryField().name(): " << wallShearStress.boundaryField().name() <<
            //   " with "
            //     << wallShearStress.boundaryField().size() << " faces. Starts at total face "
            //     << wallShearStress.boundaryField().start() << endl;    //  -- Error

            double chert_pert_2;

            //std::cin >> chert_pert_2;
        }

        Info << "patchi: " << patchi   <<   "u_star_Bf: " << u_star_Bf[patchi] << endl;

        double chert_pert_2;

        //std::cin >> chert_pert_2;


    }

//End BLOCK A */


/* BLOCK B

const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

label patch_index_ground;

forAll(mesh.boundary(), patchI)
{
    //const polyPatch& currentPatch = mesh.boundary()[patchI];
	
	const word& patchName = mesh.boundary()[patchI].name();

        if (patchName == "ground")
        {

            patch_index_ground = patchI;

            Info << "patchName = " << patchName << " patch_index_ground: " << patch_index_ground << nl
                 << endl;

            // index face avale patchI ra neshan midahad albate in index face general hast ahr patch local
            // face index ham darad
            // 7) start id global for faces patchI
            const label& faceStart = mesh.boundaryMesh()[patchI].start();
            Info << "mesh.boundaryMesh()[patchI].start(): " << faceStart <<  nl << endl;


            // 8) start id local (inside mesh) for faces patchI
            const label localFaceStart = mesh.boundaryMesh()[patchI].whichFace(faceStart);
            Info << "mesh.boundaryMesh()[patchI].whichFace(faceStart): " << localFaceStart << nl << endl;
            
            label size_face_in_patch = mesh.boundary()[patchI].size();
            Info << "mesh.boundary()[patchI].size(): " << size_face_in_patch << nl << endl;

            double chert_pert;
            std::cin >> chert_pert;
        }
}

end BLOCK B */

/*  Start : -------- block test part   method 2 --- it is not perfect but can be used in future code 
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

    //u_star.boundaryField()[patchID_ground][facei].x() = 10;
    //u_star.boundaryField()[patchID_ground][facei].y() = 20;
    //u_star.boundaryField()[patchID_ground][facei].z() = 50;
    

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
    End : -------- block test part   method 2 --- it is not perfect but can be used in future code */
 
 
        //u_star.write();   // Error 

        u_star.ref().write();

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
