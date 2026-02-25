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
dynamicFvMesh& mesh = fluidRegions[i];

immiscibleIncompressibleTwoPhaseMixture& mixture(mixtureFluid[i]);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

volScalarField& rho = rhoFluid[i];
volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];
surfaceScalarField& alphaPhi10 = alphaPhi10Fluid[i];

incompressible::turbulenceModel& turbulence = turbulenceFluid[i];
volScalarField& p_rgh = p_rghFluid[i];
volScalarField& p = pFluid[i];

volScalarField& mDot = mDotFluid[i];
volScalarField& qDot = qDotFluid[i];
volScalarField& qL = qLFluid[i];
volScalarField& qV = qVFluid[i];
volVectorField& centre = centreFluid[i];
volVectorField& normal = normalFluid[i];

volScalarField& T = TFluid[i];
volScalarField& TL = TLFluid[i];
volScalarField& TV = TVFluid[i];
volScalarField& kappa = kappaFluid[i];

const bool& alphaRestart = alphaRestartFluid[i];

scalar& cumulativeContErr = cumulativeContErrFluid[i];

volScalarField& gh = ghFluid[i];
surfaceScalarField& ghf = ghfFluid[i];
// const dimensionedScalar& ghRef = ghRefFluid[i];

// surfaceScalarField& brackFS = brackFSFluid[i];

fv::options& fvOptions = fluidFvOptions[i];

// bool frozenFlow = frozenFlowFluid[i];

const label pRefCell = pRefCellFluid[i];
const scalar pRefValue = pRefValueFluid[i];

// Use cached phase properties instead of repeated dictionary lookups
const dimensionedScalar& k1 = k1Fluid[i];
const dimensionedScalar& k2 = k2Fluid[i];
const dimensionedScalar& cp1 = cp1Fluid[i];
const dimensionedScalar& cp2 = cp2Fluid[i];
const dimensionedScalar& Tsat = TsatFluid[i];
const dimensionedScalar& hlv = hlvFluid[i];
const dimensionedScalar& diff = diffFluid[i];

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1*cp1 + alpha2*rho2*cp2,
    alpha1.boundaryField().types()
);
rhoCp.oldTime();


//Energy flux

// surfaceScalarField rhoPhiCp
// (
//     IOobject
//     (
//         "rhoPhiCp",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     fvc::interpolate(rhoCp)*phi
// );

surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

pisoControl& piso(pisoFluid[i]);

bool correctPhi
(
    piso.dict().getOrDefault("correctPhi", mesh.dynamic())
);

// bool checkMeshCourantNo
// (
//     piso.dict().getOrDefault("checkMeshCourantNo", false)
// );

// bool moveMeshOuterCorrectors
// (
//     piso.dict().getOrDefault("moveMeshOuterCorrectors", false)
// );

IOMRFZoneList& MRF = MRFfluid[i];

// surfaceScalarField& alphaPhi0(alphaPhi10Fluid[i]);

// MULES Correction
tmp<surfaceScalarField> talphaPhi1Corr0;

autoPtr<surfaceVectorField> Uf;

if (fluidRegions[i].dynamic())
{
    Info<< "Constructing face velocity Uf\n" << endl;

    Uf.reset
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(UFluid[i])
        )
    );
}

volScalarField rAU
(
    IOobject
    (
        "rAU",
        runTime.timeName(),
        fluidRegions[i],
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidRegions[i],
    dimensionedScalar("rAU", dimTime/dimDensity, 1)
);

const volScalarField TsatField
(
 IOobject
 (
	"TsatField",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
 ),
 mesh,
 Tsat
);

const volScalarField kappaL
(
 IOobject
 (
	"kappa1Field",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
 ),
 mesh,
 k1
);

const volScalarField kappaV
(
 IOobject
 (
	"kappa2Field",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
 ),
 mesh,
 k2
);

// scalar isoFaceTol = 1.0e-8;
volScalarField SArea
(
 IOobject
 (
	"SArea",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
 mesh,
 dimensionedScalar(dimless/dimLength, Zero)
 );
// --- Bubble seeding: read optional `constant/bubbleSeeding` dictionary once
//     and seed a spherical vapour region when the local temperature reaches
//     `activationTemperature`.
//
// Dictionary (optional, place in `constant/bubbleSeeding`):
// {
//     activationTemperature 373.15;    // scalar
//     coordinate (0.01 0.02 0.0);      // vector (contact point on surface)
//     surfaceNormal (0 0 1);          // vector (points into fluid)
//     contactAngle 60;                // degrees
//     initialRadius 1e-3;             // meters
//     seedOnce true;                  // optional (default true)
// }
static bool _bubbleSeedingInit = false;
static bool _bubbleSeeded = false;
static label _bubbleSeedCell = -1;
static scalar _bubbleActivationT = GREAT; // very large default = disabled
static vector _bubbleCoordinate = vector::zero;
static vector _bubbleNormal = vector::zero;
static scalar _bubbleContactAngleDeg = 90.0;
static scalar _bubbleRadius = 0.0;
static bool _bubbleSeedOnce = true;
static autoPtr<meshSearch> _bubbleMeshSearch; // Cached spatial search structure

if (!_bubbleSeedingInit)
{
    _bubbleSeedingInit = true;

    IOobject bsHeader
    (
        "bubbleSeeding",
        runTime.constant(),
        runTime,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    // Use IOobject header check for presence, then construct dictionary from file
    if (bsHeader.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary bsDict(bsHeader);

        _bubbleActivationT = bsDict.getOrDefault<scalar>("activationTemperature", GREAT);
        _bubbleCoordinate = bsDict.getOrDefault<vector>("coordinate", vector::zero);
        _bubbleNormal = bsDict.getOrDefault<vector>("surfaceNormal", vector(0,0,1));
        _bubbleContactAngleDeg = bsDict.getOrDefault<scalar>("contactAngle", 90.0);
        _bubbleRadius = bsDict.getOrDefault<scalar>("initialRadius", 0.0);
        _bubbleSeedOnce = bsDict.getOrDefault<bool>("seedOnce", true);

        if (_bubbleRadius > 0 && _bubbleActivationT < GREAT)
        {
            // Use meshSearch for efficient spatial lookup of nearest cell.
            // In parallel, only the owning processor returns a valid index.
            if (!_bubbleMeshSearch.valid())
            {
                _bubbleMeshSearch.reset(new meshSearch(mesh));
            }
            
            _bubbleSeedCell = _bubbleMeshSearch->findCell(_bubbleCoordinate);

            // Check across all ranks whether ANY processor found the cell
            label globalFound = returnReduce
            (
                (_bubbleSeedCell >= 0) ? label(1) : label(0),
                maxOp<label>()
            );

            if (globalFound == 0)
            {
                Info<< "bubbleSeeding: could not find a mesh cell for coordinate "
                    << _bubbleCoordinate << " on any processor, disabling seeding." << endl;
                _bubbleActivationT = GREAT;
            }
            else
            {
                Info<< "bubbleSeeding: configured seed"
                    << " (coord " << _bubbleCoordinate << ") with R=" << _bubbleRadius
                    << " contactAngle=" << _bubbleContactAngleDeg << " deg"
                    << " (cell local to rank "
                    << (_bubbleSeedCell >= 0 ? Pstream::myProcNo() : -1)
                    << ")" << endl;
            }
        }
        else
        {
            Info<< "bubbleSeeding: dictionary present but incomplete (must provide activationTemperature and initialRadius). Ignoring." << endl;
        }
    }
    else
    {
        // no bubbleSeeding dictionary — that's fine
			Info << "No bubbleSeeding dictionary found in " << runTime.constant() << nl;
    }
}

// Check trigger each timestep (only if configured and not yet seeded or seedOnce==false)
if (_bubbleActivationT < GREAT && (!_bubbleSeeded || !_bubbleSeedOnce))
{
    // Sample temperature at the cached/nearest cell (only on the owning
    // processor).  Then use returnReduce so that ALL processors agree on
    // whether to seed — correctBoundaryConditions() is collective MPI.
    scalar Tcell = -GREAT;  // default: not reached
    if (_bubbleSeedCell >= 0)
    {
        Tcell = T.internalField()[_bubbleSeedCell];
    }

    // Collective decision: take the maximum across all ranks
    scalar TcellGlobal = returnReduce(Tcell, maxOp<scalar>());

    if (TcellGlobal >= _bubbleActivationT)
    {
        Info << "bubbleSeeding: activation temperature reached (T=" << TcellGlobal
             << ") — inserting bubble" << endl;

        // vapour is always represented by alpha1 == 0 (user-specified)
        // compute spherical centre so the interface meets the wall at the
        // supplied contact angle: centre = contactPoint + n * (R * sin(theta))
        scalar theta = _bubbleContactAngleDeg * Foam::acos(scalar(-1.0)) / scalar(180.0);
        vector nHat = (_bubbleNormal == vector::zero) ? vector(0,0,1) : _bubbleNormal / mag(_bubbleNormal);
        vector sphereCentre = _bubbleCoordinate + nHat * (_bubbleRadius * Foam::sin(theta));

        // set alpha fields and vapour temperature inside sphere
        // — runs on ALL processors so each modifies its own local cells
        const pointField& Cc = mesh.C();
        forAll(Cc, ci)
        {
            if (mag(Cc[ci] - sphereCentre) <= _bubbleRadius)
            {
                // alpha1 == 0 => vapour, alpha2 == 1 => liquid complement
                alpha1[ci] = 0.0;
                alpha2[ci] = 1.0;

                // set vapour temperature to saturation for seeded vapour cells
                TV[ci] = TsatField[ci];

                // keep mixture/diagnostic temperature consistent
                T[ci] = TV[ci];
            }
        }

        // mark seeded to avoid repeated insertion (unless seedOnce==false)
        _bubbleSeeded = true;

        // update boundary fields — collective MPI, now called by ALL ranks
        alpha1.correctBoundaryConditions();
        alpha2.correctBoundaryConditions();
        TV.correctBoundaryConditions();
        T.correctBoundaryConditions();
        Info << "bubbleSeeding: bubble inserted (centre=" << sphereCentre
             << ", R=" << _bubbleRadius << ")" << endl;
    }
}
