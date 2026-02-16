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

// Compute vapour/liquid volume-weighted redistribution in-place
// to avoid creating 4+ full-mesh temporary fields
{
	dimensionedScalar NNum = fvc::domainIntegrate(mDot);

	// Accumulate weighted integrals without separate cutoff fields
	scalar NVDenVal = 0.0;
	scalar NLDenVal = 0.0;

	const scalarField& alpha1In = alpha1.primitiveField();
	const scalarField& mDotIn = mDot.primitiveField();
	const scalarField& Vol = mesh.V();

	forAll(alpha1In, celli)
	{
		const scalar epsilonV = 0.001 - alpha1In[celli];
		if (epsilonV > 0)
		{
			NVDenVal += mDotIn[celli] * epsilonV * Vol[celli];
		}
		const scalar epsilonL = alpha1In[celli] - 0.999;
		if (epsilonL > 0)
		{
			NLDenVal += mDotIn[celli] * epsilonL * Vol[celli];
		}
	}

	reduce(NVDenVal, sumOp<scalar>());
	reduce(NLDenVal, sumOp<scalar>());

	const scalar NV = (NVDenVal != 0) ? NNum.value() / NVDenVal : 2.0;
	const scalar NL = (NLDenVal != 0) ? NNum.value() / NLDenVal : 2.0;

	Info << "NV = " << NV << ", NL = " << NL << endl;

	// Rewrite mDot in a single pass
	scalarField& mDotIn_ref = mDot.primitiveFieldRef();
	forAll(alpha1In, celli)
	{
		scalar val = 0.0;
		const scalar epsilonV = 0.001 - alpha1In[celli];
		if (epsilonV > 0)
		{
			val += NV * mDotIn[celli] * epsilonV;
		}
		const scalar epsilonL = alpha1In[celli] - 0.999;
		if (epsilonL > 0)
		{
			val -= NL * mDotIn[celli] * epsilonL;
		}
		mDotIn_ref[celli] = val;
	}
}

const scalar mEvap = gSum( qDot/hlv * mesh.V());
  
Info<< "Timeev= " << runTime.timeName() <<" " <<
			"mEvap= " << mEvap  <<" " <<
			nl << endl;
