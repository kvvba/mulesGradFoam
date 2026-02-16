// Persistent lists (allocated once in createFluidFields.H, cleared each timestep)
boolList& interfaceCells = interfaceCellsFluid[i];
interfaceCells.setSize(mesh.nCells());
interfaceCells = false;
DynamicList<label>& VapourCellsLabels = VapourCellsLabelsFluid[i];
VapourCellsLabels.clear();
DynamicList<label>& LiquidCellsLabels = LiquidCellsLabelsFluid[i];
LiquidCellsLabels.clear();
boolList& vapourCells = vapourCellsFluid[i];
vapourCells.setSize(mesh.nCells());
vapourCells = false;
boolList& liquidCells = liquidCellsFluid[i];
liquidCells.setSize(mesh.nCells());
liquidCells = false;

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

implicitInterfaceDiffFlux& impDiffFlux = impDiffFluxFluid[i];

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
