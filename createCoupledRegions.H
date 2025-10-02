/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
		Copyright (C) 2023 Jakub Cranmer <jakub@posteo.net>
-------------------------------------------------------------------------------
License
    This file is part of interMultiKFoam.

    interMultiKFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    interMultiKFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with interMultiKFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
fvSolution solutionDict(runTime);

bool coupled(solutionDict.getOrDefault("coupledEnergyField", false));

autoPtr<fvMatrix<scalar>> fvMatrixAssemblyPtr;

forAll(fluidRegions, i)
{
    const rhoThermo& thermo = refCast<const rhoThermo>(thermoFluid[i]);
    const auto& bpsi = thermo.T().boundaryField();

    forAll(bpsi, patchI)
    {
        if (bpsi[patchI].useImplicit())
        {
            coupled = true;
        }
    }
}

forAll(solidRegions, i)
{
    solidThermo& thermo = thermos[i];
    const auto& bpsi = thermo.T().boundaryField();

    forAll(bpsi, patchI)
    {
        if (bpsi[patchI].useImplicit())
        {
            coupled = true;
        }
    }
}

forAll(fluidRegions, i)
{
    const rhoThermo& thermo = refCast<const rhoThermo>(thermoFluid[i]);
    if (coupled)
    {
        Info<< "Create fvMatrixAssembly." << endl;
        fvMatrixAssemblyPtr.reset
        (
            new fvMatrix<scalar>
            (
                thermo.he(),
                dimEnergy/dimTime
            )
        );
        break;
    }
}

if (coupled && !fvMatrixAssemblyPtr)
{
    forAll(solidRegions, i)
    {
        solidThermo& thermo = thermos[i];
        Info<< "Create fvMatrixAssembly." << endl;
        fvMatrixAssemblyPtr.reset
        (
            new fvMatrix<scalar>
            (
                thermo.he(),
                dimEnergy/dimTime
            )
        );
        break;
    }
}
