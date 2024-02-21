/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 76 "/home/saeb/OpenFoam_run/20240214_wall_ustar_shear_stress_sahand/20240215_run/0/U/boundaryField/inlet/#codeStream"
#include "fvCFD.H"
			#include <math.h>   ///////////
			#include <cmath>    //////////
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    void codeStream_f381c06b0140f78f3b4f53873d3c9294fb13bab7
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 92 "/home/saeb/OpenFoam_run/20240214_wall_ustar_shear_stress_sahand/20240215_run/0/U/boundaryField/inlet/#codeStream"
const IOdictionary& d = static_cast<const IOdictionary&>
			(
				dict.parent().parent()
			);
			const fvMesh& mesh = refCast<const fvMesh>(d.db());
			const label id = mesh.boundary().findPatchID("inlet");
			const fvPatch& patch = mesh.boundary()[id];
			vectorField U(patch.size(), vector(0, 0, 0));
			//const scalar  pi = constant::mathematical::pi;
			//const scalar  U_0 = 2.; //maximum velocity
			//const scalar  p_ctr = 8.; //patch center     
			//const scalar  p_r = 8.; //patch radius
			
			
			const scalar  PI = constant::mathematical::pi;

			
			
			forAll(U, i)	//equivalent to for (int i=0; patch.size()<i; i++)
			{
				//const scalar y = patch.Cf()[i][1];
				
				//double aa1 = ((y+(5.000000e-04)) / (5.000000e-04)) +10.;
				
				const double y = patch.Cf()[i][1];
				
				//double y_01 = -10;
				//double y_elevation_from0 = y + (-y_01);
				double y_elevation_from0 = y + (-1*-20);

								
				//double aa1 = 12.+ y;   //( (5.000000e-04+ y) / (5.000000e-04))+10.;    ;


				//double u_x =  log(aa1));
				
				//double aa =log(10.1);
				
				//double aa2 =log(aa1);
				
				// U[i] = vector(((8.500000e-01) / (4.200000e-01))* (log((y+(5.000000e-04)) / (5.000000e-04))), 0., 0.);   // vector(($u_star/4.200000e-01)*log((y+5.000000e-04)/5.000000e-04), 0., 0.);  // vector (20, 0 ,0);   

                // U[i] = vector(((8.500000e-01) / (4.200000e-01) )*($a*fabs((pow(y,3))))+($b*(pow(y,2))), 0., 0.);

				// U[i] = vector(2, 0., 0.);
				
				std::ofstream write_y;
	            write_y.open("write_y_inlet_result_elevation_from0", std::ofstream::out | std::ofstream::app);  
	            write_y.precision(20);
	            write_y << y_elevation_from0  << "\n";
	            write_y.close();
				
								
				/*if (y < 10) {


                double taylor_series_x =(y/z_01);
				
				//scalar u_x = ((u_star1) / (k_van_karman1))* ((taylor_series_x)-(pow(taylor_series_x,2.)/2.)+(pow(taylor_series_x,3)/3.));
                
				double u_x=0;

			    }
			    else {
				
				//scalar u_x = ((u_star1) / (k_van_karman1))* (log((y+(z_01)) / (z_01)));
                
				double u_x=10;

			    }*/
				
				double teta_wind_velocity_radian = (0*PI)/180;
				
				 scalar u_x = ((8.500000e-01) / (4.200000e-01))* (log((y_elevation_from0+(5.000000e-04)) / (5.000000e-04))) * cos(teta_wind_velocity_radian);
				 scalar u_z = ((8.500000e-01) / (4.200000e-01))* (log((y_elevation_from0+(5.000000e-04)) / (5.000000e-04))) * sin(teta_wind_velocity_radian);
				 
				 U[i] = vector(u_x, 0., u_z);
			}
			
			//U.writeEntry("", os);  // Error
			
			writeEntry(os, "", U);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

