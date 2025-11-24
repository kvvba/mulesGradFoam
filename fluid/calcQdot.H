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
    This file is part of mulesGradFoam.

    mulesGradFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mulesGradFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mulesGradFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
boolList interfaceCells(alpha1.mesh().nCells(), false);
labelList VapourCellsLabels;
labelList LiquidCellsLabels;
boolList vapourCells(alpha1.mesh().nCells(), false);
boolList liquidCells(alpha1.mesh().nCells(), false);

scalarField ap
(
 volPointInterpolation::New(mesh).interpolate(mixture.alpha1())
);

cutCellIso cutCell(mesh,ap);

forAll (alpha1, celli)
{
	label status = cutCell.calcSubCell(celli,0.5);
	if (status == 0) // cell is cut
		{
			// interfaceLabels.append(celli);
			interfaceCells[celli] = true;
			normal[celli] = cutCell.faceArea();
			centre[celli] = cutCell.faceCentre();
			SArea[celli] =
				mag(cutCell.faceArea())/mesh.V()[celli];
			kappa[celli] = k1.value() * alpha1[celli] + k2.value() * alpha2[celli];
		}
	else
		{
			normal[celli] = Zero;
			centre[celli] = Zero;
			interfaceCells[celli] = false;
		}
	
	if (status == -1)							// liquid
		{
			liquidCells[celli] = true;
			LiquidCellsLabels.append(celli);
			kappa[celli] = k1.value();
		}
	else if (status == 1)					// vapour
		{
			vapourCells[celli] = true;
			VapourCellsLabels.append(celli);
			kappa[celli] = k2.value();
		}
}

implicitInterfaceDiffFlux impDiffFlux(mesh);

tmp<volScalarField> energyFlux1
(
 impDiffFlux.diffusiveFlux
 (
	interfaceCells,
	normal,
	centre,
	TL,
	kappaL,
	TsatField,
	false
 )
);

tmp<volScalarField> energyFlux2
(
 impDiffFlux.diffusiveFlux
 (
	interfaceCells,
	normal,
	centre,
	TV,
	kappaV,
	TsatField,
	true
 )
);

tmp<fvScalarMatrix> TsourceL
(
 impDiffFlux.diffusiveFlux
 (
	interfaceCells,
	normal,
	centre,
	TL,
	kappaL,
	TsatField,
	false,
	dimPower
 )
);

tmp<fvScalarMatrix> TsourceV
(
 impDiffFlux.diffusiveFlux
 (
	interfaceCells,
	normal,
	centre,
	TV,
	kappaV,
	TsatField,
	true,
	dimPower
 )
);

qL = energyFlux1 * SArea;
qV = energyFlux2 * SArea;

qDot = qL + qV;

Info << "max qDot = " << max(qDot).value() << nl;
