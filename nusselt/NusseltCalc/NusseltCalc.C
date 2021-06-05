/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    wallHeatFlux

Description
    Calculates and writes the heat flux for all patches as the boundary field
    of a volScalarField and also prints the integrated flux for all wall
    patches.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "wallFvPatch.H"
#include "basicThermo.H"
//#include "readRefValues.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        #include "createFields.H"

        surfaceScalarField heatFlux
        (
            fvc::interpolate
            (
                (
                    turbulence.valid()
                  ? turbulence->alphaEff()()
                  : thermo->alpha()
                )
            )*fvc::snGrad(h)
        );

        const surfaceScalarField::Boundary& patchHeatFlux =
            heatFlux.boundaryField();

        const volScalarField::Boundary& patchRadHeatFlux =
            Qr.boundaryField();

        const surfaceScalarField::Boundary& magSf =
            mesh.magSf().boundaryField();

        Info<< "\nWall heat fluxes [W]" << endl;
        forAll(patchHeatFlux, patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                scalar convFlux = gSum(magSf[patchi]*patchHeatFlux[patchi]);
                scalar radFlux = -gSum(magSf[patchi]*patchRadHeatFlux[patchi]);

                Info<< mesh.boundary()[patchi].name() << endl
                    << "    convective: " << convFlux << endl
                    << "    radiative:  " << radFlux << endl
                    << "    total:      " << convFlux + radFlux << endl;
            }
        }
        Info<< endl;

        volScalarField wallHeatFlux
        (
            IOobject
            (
                "wallHeatFlux",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("wallHeatFlux", heatFlux.dimensions(), 0.0)
        );

        volScalarField::Boundary& wallHeatFluxBf =
            wallHeatFlux.boundaryFieldRef();

        forAll(wallHeatFluxBf, patchi)
        {
            wallHeatFluxBf[patchi] = patchHeatFlux[patchi];
        }

        wallHeatFlux.write();

//////////
    Info << "Reading refValues\n" << endl;

    IOdictionary refValues
    (
        IOobject
        (
            "refValues",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar k (readScalar(refValues.lookup("k")));
    Info << "Conductivity is:"<< k << endl;
    scalar T_initial(readScalar(refValues.lookup("T_initial")));
    Info << "Initial temperature is:"<< T_initial << endl;
    scalar T_hot(readScalar(refValues.lookup("T_hot")));
    Info << "Hot wall temperature:"<< T_hot << endl;
    scalar length(readScalar(refValues.lookup("length")));
    Info << "Length scale is set to:"<< length << endl;    

//////////
/////////////////////////////
volScalarField NusseltNumber
        (
                IOobject
                (
                        "NusseltNumber",
                        runTime.timeName(),
                        mesh
                ),
                mesh,
                dimensionedScalar("NusseltNumber", heatFlux.dimensions(), 0.0)
        );

        forAll(NusseltNumber.boundaryFieldRef(), patchi)
        {
                NusseltNumber.boundaryFieldRef()[patchi] = length*
                patchHeatFlux[patchi]/((T_hot-T_initial)*k);
        }
        NusseltNumber.write();

/////////////////////////////
        // Write the total heat-flux including the radiative contribution
        // if available
        if (Qr.headerOk())
        {
            volScalarField totalWallHeatFlux
            (
                IOobject
                (
                    "totalWallHeatFlux",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar
                (
                    "totalWallHeatFlux",
                    heatFlux.dimensions(),
                    0.0
                )
            );

            volScalarField::Boundary& totalWallHeatFluxBf =
                totalWallHeatFlux.boundaryFieldRef();

            forAll(totalWallHeatFluxBf, patchi)
            {
                totalWallHeatFluxBf[patchi] =
                    patchHeatFlux[patchi] - patchRadHeatFlux[patchi];
            }

            totalWallHeatFlux.write();
        }
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
