fvSolution solutionDict(runTime);

const dictionary& presimControl = solutionDict.subDict("presimControl");

const Switch computeFluxes =
				presimControl.getOrDefault<Switch>("computeFluxes", false);

bool converged(false);

const Switch stopWhenConverged =
				presimControl.getOrDefault<Switch>("stopWhenConverged", false);

const scalar fluxTol(presimControl.getOrDefault<double>("fluxTolerance", 1e-8));
const scalar TTol(presimControl.getOrDefault<double>("TTolerance", 1e-8));
