/*IOobject wallShearStressMag(
    "wallShearStressMag",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE  // No need to write during the simulation, as we'll write it once outside the loop
);

IOobject u_star_1(
    "u_star_1",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
);*/




Info<< "Creating field u_star\n" << endl;

volVectorField u_star
(
    IOobject
    (
        "u_star",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Creating field u_star\n" << endl;

volVectorField u_star_ref
(
    IOobject
    (
        "u_star_ref",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//------------- this part is not needed ----------------------



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

//---------------------------------------------------------------*/

/*------- reading wallShearStress ---------------------------*/

Info<< "Reading field wallShearStress\n" << endl;

volVectorField wallShearStress
(
    IOobject
    (
        "wallShearStress",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE    //AUTO_WRITE
    ),
    mesh
);

tmp<volScalarField> tmp_volScalarField
 (
 new volScalarField
 (
 IOobject
 (
 "tmp_volScalarField",
 runTime.timeName(),
 mesh,
 IOobject::NO_READ,
 IOobject::NO_WRITE
 ),
 mesh,
 dimensionedScalar("0", dimMass/dimLength/dimTime, 0.0)
 )
 );


tmp<volVectorField> tmp_volVectorField
(
new volVectorField
(
IOobject
(
"tmp_volVectorField",
mesh.time().timeName(),
mesh,
IOobject::NO_READ,
IOobject::NO_WRITE
),
mesh,
dimensionedVector
(
"zero",
dimensionSet(1, -2, -2, 0, 0),
vector::zero
)
)
);

/*------- End reading wallShearStress ---------------------------*/


Info<< "Reading physicalProperties\n" << endl;

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading diffusivity rho\n" << endl;

dimensionedScalar rho
(
    physicalProperties.lookup("rho")
);

Info<< "Reading diffusivity gravity_mag\n" << endl;

dimensionedScalar gravity_mag
(
    physicalProperties.lookup("gravity_mag")
);


// #include "createPhi.H"

// #include "createFvModels.H"
// #include "createFvConstraints.H"
