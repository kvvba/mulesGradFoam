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
// Initialise fluid field pointer lists
PtrList<immiscibleIncompressibleTwoPhaseMixture> mixtureFluid
(
 fluidRegions.size()
);

PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<surfaceScalarField> alphaPhi10Fluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<incompressible::turbulenceModel> turbulenceFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volScalarField> TFluid(fluidRegions.size());
PtrList<volScalarField> TLFluid(fluidRegions.size());
PtrList<volScalarField> TVFluid(fluidRegions.size());
PtrList<volVectorField> qwFluid(fluidRegions.size());

PtrList<volScalarField> kappaFluid(fluidRegions.size());
PtrList<volScalarField> mDotFluid(fluidRegions.size());
PtrList<volScalarField> qDotFluid(fluidRegions.size());
PtrList<volScalarField> qVFluid(fluidRegions.size());
PtrList<volScalarField> qLFluid(fluidRegions.size());
PtrList<volVectorField> centreFluid(fluidRegions.size());
PtrList<volVectorField> normalFluid(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());
List<scalar> cumulativeContErrFluid(fluidRegions.size());
PtrList<dimensionedScalar> ghRefFluid(fluidRegions.size());
// List<bool> frozenFlowFluid(fluidRegions.size(),false);
List<bool> alphaRestartFluid(fluidRegions.size(),false);

PtrList<pisoControl> pisoFluid(fluidRegions.size());

// PtrList <surfaceScalarField> brackFSFluid(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

scalarField maxTFluid(fluidRegions.size());

// Persistent material constants per fluid region (read once at startup)
PtrList<dimensionedScalar> k1Fluid(fluidRegions.size());
PtrList<dimensionedScalar> k2Fluid(fluidRegions.size());
PtrList<dimensionedScalar> cp1Fluid(fluidRegions.size());
PtrList<dimensionedScalar> cp2Fluid(fluidRegions.size());
PtrList<dimensionedScalar> TsatFluid(fluidRegions.size());
PtrList<dimensionedScalar> hlvFluid(fluidRegions.size());
PtrList<dimensionedScalar> diffFluid(fluidRegions.size());

// Persistent fields (avoid per-timestep re-allocation)
PtrList<volScalarField> TsatFieldFluid(fluidRegions.size());
PtrList<volScalarField> kappaLFluid(fluidRegions.size());
PtrList<volScalarField> kappaVFluid(fluidRegions.size());
PtrList<volScalarField> SAreaFluid(fluidRegions.size());
PtrList<volScalarField> rhoCpLFluid(fluidRegions.size());
PtrList<volScalarField> rhoCpVFluid(fluidRegions.size());

// Persistent calcQdot lists (avoid per-timestep allocation)
List<boolList> interfaceCellsFluid(fluidRegions.size());
List<boolList> vapourCellsFluid(fluidRegions.size());
List<boolList> liquidCellsFluid(fluidRegions.size());
List<DynamicList<label>> VapourCellsLabelsFluid(fluidRegions.size());
List<DynamicList<label>> LiquidCellsLabelsFluid(fluidRegions.size());
PtrList<implicitInterfaceDiffFlux> impDiffFluxFluid(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
	Info<< "*** Reading fluid mesh thermophysical properties for region "
			<< fluidRegions[i].name() << nl << endl;

	Info<< "    Adding to p_rghFluid\n" << endl;
	p_rghFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "p_rgh",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		 )
    );


	Info<< "    Adding to UFluid\n" << endl;
	UFluid.set
    (
		 i,
		 new volVectorField
		 (
			IOobject
			(
			 "U",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		 )
    );

	Info<< "    Adding to phiFluid\n" << endl;
	phiFluid.set
    (
		 i,
		 new surfaceScalarField
		 (
			IOobject
			(
			 "phi",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::READ_IF_PRESENT,
			 IOobject::AUTO_WRITE
			),
			linearInterpolate(UFluid[i])
			& fluidRegions[i].Sf()
		 )
    );

	Info<< "    Adding to hRefFluid\n" << endl;
	hRefFluid.set
    (
		 i,
		 new uniformDimensionedScalarField
		 (
			IOobject
			(
			 "hRef",
			 runTime.constant(),
			 fluidRegions[i],
			 IOobject::READ_IF_PRESENT,
			 IOobject::NO_WRITE
			),
			dimensionedScalar("hRef", dimLength, Zero) // uses name
		 )
    );

	Info << "    Adding to mDot\n" << endl;
	mDotFluid.set
		(
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "mDot",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		 )
		);

	Info << "    Adding to qDot\n" << endl;
	qDotFluid.set
		(i,
		 new volScalarField
		 (
			IOobject
			(
			 "qDot",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i],
			dimensionedScalar
			(
			 "qDot",
			 dimensionSet(1,-1,-3,0,0,0,0),
			 0.0
			)
		 )
		);

		qVFluid.set
		(i,
		 new volScalarField
		 (
			IOobject
			(
			 "qV",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i],
			dimensionedScalar
			(
			 "qV",
			 dimensionSet(1,-1,-3,0,0,0,0),
			 0.0
			)
		 )
		);
	qLFluid.set
		(i,
		 new volScalarField
		 (
			IOobject
			(
			 "qL",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i],
			dimensionedScalar
			(
			 "qL",
			 dimensionSet(1,-1,-3,0,0,0,0),
			 0.0
			)
		 )
		);

	centreFluid.set
		(
		 i,
		 new volVectorField
		 (
			IOobject
			(
			 "centre",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i],
			dimensionedVector(dimLength, Zero)
		 )
		);

	normalFluid.set
		(
		 i,
		 new volVectorField
		 (
			IOobject
			(
			 "normal",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i],
			dimensionedVector(dimArea, Zero)
		 )
		);

	dimensionedScalar ghRef
		(
		 mag(g.value()) > SMALL
		 ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
		 : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
		);

	ghRefFluid.set(i, new dimensionedScalar(ghRef));

	Info<< "    Adding to ghFluid\n" << endl;
	ghFluid.set
    (
		 i,
		 new volScalarField
		 (
			"gh",
			(g & fluidRegions[i].C()) - ghRef
		 )
    );

	Info<< "    Adding to ghfFluid\n" << endl;
	ghfFluid.set
    (
		 i,
		 new surfaceScalarField
		 (
			"ghf",
			(g & fluidRegions[i].Cf()) - ghRef
		 )
    );

	mixtureFluid.set
    (
		 i,
		 new immiscibleIncompressibleTwoPhaseMixture
		 (
			UFluid[i],
			phiFluid[i]
		 )
    );


	Info<< "    Adding to turbulenceFluid\n" << endl;
	turbulenceFluid.set
    (
		 i,
		 incompressible::turbulenceModel::New
		 (
			UFluid[i],
			phiFluid[i],
			mixtureFluid[i]
		 ).ptr()
    );


	immiscibleIncompressibleTwoPhaseMixture& mixture(mixtureFluid[i]);
	volScalarField& alpha1(mixture.alpha1());
	volScalarField& alpha2(mixture.alpha2());


	//- This also builds the surface reconstruction scheme

	const dimensionedScalar& rho1 = mixture.rho1();
	const dimensionedScalar& rho2 = mixture.rho2();

	dictionary phase1 = mixture.subDict(mixture.phase1Name());
	dictionary phase2 = mixture.subDict(mixture.phase2Name());

	dimensionedScalar k1
    (
		 "kappa",
		 dimensionSet(1,1,-3,-1,0,0,0),
		 readScalar(phase1.lookup("kappa"))
    );

	dimensionedScalar k2
    (
		 "kappa",
		 dimensionSet(1,1,-3,-1,0,0,0),
		 readScalar(phase2.lookup("kappa"))
    );

	Info<< "    Adding to kappaFluid\n" << endl;
	kappaFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "kappa",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			 ),
			k1*alpha1 + k2*alpha2
			)
		 );

	// Store material constants persistently (read once at startup)
	{
		dimensionedScalar cp1_init
		(
			"Cp",
			dimensionSet(0,2,-2,-1,0,0,0),
			readScalar(phase1.lookup("Cp"))
		);
		dimensionedScalar cp2_init
		(
			"Cp",
			dimensionSet(0,2,-2,-1,0,0,0),
			readScalar(phase2.lookup("Cp"))
		);
		dictionary pcProperties = mixture.subDict("PhaseChangeProperties");

		k1Fluid.set(i, new dimensionedScalar(k1));
		k2Fluid.set(i, new dimensionedScalar(k2));
		cp1Fluid.set(i, new dimensionedScalar(cp1_init));
		cp2Fluid.set(i, new dimensionedScalar(cp2_init));
		TsatFluid.set
		(
			i,
			new dimensionedScalar
			(
				"Tsat",
				dimTemperature,
				readScalar(pcProperties.lookup("Tsat"))
			)
		);
		hlvFluid.set
		(
			i,
			new dimensionedScalar
			(
				"hlv",
				dimensionSet(0,2,-2,0,0,0,0),
				readScalar(pcProperties.lookup("hlv"))
			)
		);
		diffFluid.set
		(
			i,
			new dimensionedScalar
			(
				"diff",
				dimensionSet(0,2,0,0,0,0,0),
				readScalar(pcProperties.lookup("diff"))
			)
		);

		// Persistent constant fields
		TsatFieldFluid.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"TsatField",
					runTime.timeName(),
					fluidRegions[i],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				fluidRegions[i],
				TsatFluid[i]
			)
		);
		kappaLFluid.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"kappa1Field",
					runTime.timeName(),
					fluidRegions[i],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				fluidRegions[i],
				k1
			)
		);
		kappaVFluid.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"kappa2Field",
					runTime.timeName(),
					fluidRegions[i],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				fluidRegions[i],
				k2
			)
		);

		// Persistent reusable fields
		SAreaFluid.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"SArea",
					runTime.timeName(),
					fluidRegions[i],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				fluidRegions[i],
				dimensionedScalar(dimless/dimLength, Zero)
			)
		);
		rhoCpLFluid.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"rhoCpL",
					runTime.timeName(),
					fluidRegions[i],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				fluidRegions[i],
				rho1 * cp1_init
			)
		);
		rhoCpVFluid.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"rhoCpV",
					runTime.timeName(),
					fluidRegions[i],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				fluidRegions[i],
				rho2 * cp2_init
			)
		);

		// Persistent calcQdot lists
		interfaceCellsFluid[i].setSize(fluidRegions[i].nCells(), false);
		vapourCellsFluid[i].setSize(fluidRegions[i].nCells(), false);
		liquidCellsFluid[i].setSize(fluidRegions[i].nCells(), false);
		VapourCellsLabelsFluid[i].reserve(fluidRegions[i].nCells()/2);
		LiquidCellsLabelsFluid[i].reserve(fluidRegions[i].nCells()/2);
		impDiffFluxFluid.set(i, new implicitInterfaceDiffFlux(fluidRegions[i]));
	}

	Info<< "    Adding to TFluid\n" << endl;
	TFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "T",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		 )
    );

	TLFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "TL",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		 )
    );
	TVFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "TV",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		 )
    );

		
	maxTFluid[i] = gMax(TFluid[i]);

	rhoFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "rho",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			alpha1*rho1 + alpha2*rho2
		 )
    );

	pisoFluid.set
    (
		 i,
		 new pisoControl(fluidRegions[i])
    );

	// pimpleFluid[i].dict().readIfPresent("frozenFlow", frozenFlowFluid[i]);

	pFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "p",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			p_rghFluid[i] + rhoFluid[i]*ghFluid[i]
		 )
    );

	pRefCellFluid[i] = -1;
	pRefValueFluid[i] = 0.0;
	setRefCell
    (
		 pFluid[i],
		 p_rghFluid[i],
		 pisoFluid[i].dict(),
		 pRefCellFluid[i],
		 pRefValueFluid[i]
    );



	if (p_rghFluid[i].needReference())
    {
			pFluid[i] += dimensionedScalar
        (
				 "p",
				 pFluid[i].dimensions(),
				 pRefValueFluid[i] - getRefCellValue(pFluid[i], pRefCellFluid[i])
        );
			p_rghFluid[i] = pFluid[i] - rhoFluid[i]*ghFluid[i];
    }

	fluidRegions[i].setFluxRequired(p_rghFluid[i].name());
	fluidRegions[i].setFluxRequired(alpha1.name());


	Info<< "    Adding fvOptions\n" << endl;
	fluidFvOptions.set
    (
		 i,
		 new fv::options(fluidRegions[i])
    );

	turbulenceFluid[i].validate();

	Info<< "    Adding MRF\n" << endl;
	MRFfluid.set
    (
		 i,
		 new IOMRFZoneList(fluidRegions[i])
    );

	IOobject alphaPhi10Header
    (
		 IOobject::groupName("alphaPhi0", alpha1.group()),
		 runTime.timeName(),
		 fluidRegions[i],
		 IOobject::READ_IF_PRESENT,
		 IOobject::AUTO_WRITE
    );

	alphaRestartFluid[i] =
		alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true);

	if (alphaRestartFluid[i])
    {
			Info << "Restarting alpha" << endl;
    }

	alphaPhi10Fluid.set
    (
		 i,
		 new surfaceScalarField
		 (
			alphaPhi10Header,
			phiFluid[i]*fvc::interpolate(alpha1)
		 )
    );

	// brackFSFluid.set
	//   (
	//    i,
	//    new surfaceScalarField
	//    (
	//     fvc::interpolate(mixture.sigmaK())*fvc::snGrad(alpha1)
	//     )
	//    );

	uniformDimensionedScalarField cumulativeContErrIO
    (
		 IOobject
		 (
			"cumulativeContErr",
			runTime.timeName(),
			"uniform",
			fluidRegions[i],
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		 ),
		 dimensionedScalar(dimless, Zero)
    );


	cumulativeContErrFluid[i] = cumulativeContErrIO.value();


	if(mixture.found("writeQwall"))
    {
			qwFluid.set
        (
				 i,
				 new volVectorField
				 (
					IOobject
					(
					 "qw",
					 runTime.timeName(),
					 fluidRegions[i],
					 IOobject::NO_READ,
					 IOobject::AUTO_WRITE
					),
					-kappaFluid[i]*fvc::grad(TFluid[i])
				 )
        );
    }

}
