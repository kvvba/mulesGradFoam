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

surfaceScalarField kappaLf(fvc::interpolate(kappaL));
surfaceScalarField kappaVf(fvc::interpolate(kappaV));

const labelList& owner = mesh.owner();
const labelList& neighbour = mesh.neighbour();

forAll(kappaLf, facei)
{
	const label& own = owner[facei];
	const label& nei = neighbour[facei];
	const scalar vapourOwn = vapourCells[own];
	const scalar vapourNei = vapourCells[nei];
	const scalar liquidOwn = liquidCells[own];
	const scalar liquidNei = liquidCells[nei];

	if (vapourOwn || vapourNei)
		// cell or neighbour cell is on the vapour side
		{
			kappaLf[facei] = 0.0;
		}
	if (liquidOwn || liquidNei)
		// cell or neighbour cell is on the liquid side
		{
			kappaVf[facei] = 0.0;
		}
}

volScalarField rhoCpL
(
 IOobject
 (
	"rhoCpL",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
 ),
 mesh,
 rho1 * cp1
);

volScalarField rhoCpV
(
 IOobject
 (
	"rhoCpV",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
 ),
 mesh,
 rho2 * cp2
);

forAll(rhoCpL, celli)
{
	if (interfaceCells[celli])
		{
			rhoCpL[celli] = alpha1[celli] * rho1.value() * cp1.value();
			rhoCpV[celli] = alpha2[celli] * rho2.value() * cp2.value();
		}
	else if (vapourCells[celli])
		{
			rhoCpL[celli] = rho1.value() * cp1.value() * 1.0e-7;
		}
	else if (liquidCells[celli])
		{
			rhoCpV[celli] = rho2.value() * cp2.value() * 1.0e-7;
		}
}

surfaceScalarField rhoPhiCpL
("rhoPhiCpL",
 fvc::interpolate(rhoCpL)*phi);

surfaceScalarField rhoPhiCpV
("rhoPhiCpV",
 fvc::interpolate(rhoCpV)*phi);

fvScalarMatrix TLEqn
(
 fvm::ddt(rhoCpL,TL)
 + fvm::div(rhoPhiCpL,TL)
 - fvm::Sp(fvc::ddt(rhoCpL) + fvc::div(rhoPhiCpL), TL)
 - fvm::laplacian(kappaLf,TL)
 + TsourceL
);

fvScalarMatrix TVEqn
(
 fvm::ddt(rhoCpV,TV)
 + fvm::div(rhoPhiCpV,TV)
 - fvm::Sp(fvc::ddt(rhoCpV) + fvc::div(rhoPhiCpV), TV)
 - fvm::laplacian(kappaVf,TV)
 + TsourceV
);

TLEqn.setValues(VapourCellsLabels, TsatField);
TVEqn.setValues(LiquidCellsLabels, TsatField);

TL.correctBoundaryConditions();
TV.correctBoundaryConditions();

TLEqn.relax();
TVEqn.relax();

TLEqn.solve();
TVEqn.solve();

T = TL*alpha1 + TV*alpha2;

T.correctBoundaryConditions();

mixture.correct();

Info << "Min(TL) = " << min(TL).value() << " Max(TL) = " << max(TL).value() << endl;
Info << "Min(TV) = " << min(TV).value() << " Max(TV) = " << max(TV).value() << endl;

