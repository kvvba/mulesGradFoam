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

mDot = qDot / hlv;

fvScalarMatrix mDotEqn
(
 fvm::laplacian(diff, mDot)
 - fvm::Sp(scalar(1), mDot)
 ==
 - qDot / hlv
 );

mDotEqn.solve();

volScalarField epsilonV = 0.001 - alpha1;
volScalarField vapourCutoff = pos(epsilonV);

volScalarField epsilonL = alpha1 - 0.999;
volScalarField liquidCutoff = pos(epsilonL);
	
dimensionedScalar NNum = fvc::domainIntegrate(mDot);
	
volScalarField vapourVolGen = mDot * epsilonV * vapourCutoff;
dimensionedScalar NVDen = fvc::domainIntegrate(vapourVolGen);

volScalarField liquidVolGen = mDot * epsilonL * liquidCutoff;
dimensionedScalar NLDen = fvc::domainIntegrate(liquidVolGen);

dimensionedScalar NV("NV", dimless, 2.0);
dimensionedScalar NL("NL", dimless, 2.0);

if(NVDen.value() != 0){NV = NNum/NVDen;}
if(NLDen.value() != 0){NL = NNum/NLDen;}

// dimensionedScalar NV = NNum / NVDen;
// dimensionedScalar NL = NNum / NLDen;

Info << "NV = " << NV.value() << ", NL = " << NL.value() << endl;

vapourVolGen = NV * vapourCutoff * mDot * epsilonV;
liquidVolGen = - NL * liquidCutoff * mDot * epsilonL;

mDot = vapourVolGen + liquidVolGen;

const scalar mEvap = gSum( qDot/hlv * mesh.V());
  
Info<< "Timeev= " << runTime.timeName() <<" " <<
			"mEvap= " << mEvap  <<" " <<
			nl << endl;
