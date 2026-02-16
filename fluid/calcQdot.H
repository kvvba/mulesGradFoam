boolList interfaceCells(alpha1.mesh().nCells(), false);
DynamicList<label> VapourCellsLabels(alpha1.mesh().nCells()/2);
DynamicList<label> LiquidCellsLabels(alpha1.mesh().nCells()/2);
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
			rhoCp[celli] = rho1.value() * cp1.value() * alpha1[celli] + rho2.value() * cp2.value() * alpha2[celli];
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
			rhoCp[celli] = rho1.value() * cp1.value();
		}
	else if (status == 1)					// vapour
		{
			vapourCells[celli] = true;
			VapourCellsLabels.append(celli);
			kappa[celli] = k2.value();
			rhoCp[celli] = rho2.value() * cp2.value();
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
