fvSolution solutionDict(runTime);

const dictionary* presimControlPtr = solutionDict.findDict("presimControl");

const Switch computeFluxes =
    presimControlPtr ? presimControlPtr->getOrDefault<Switch>("computeFluxes", false) : Switch(false);

bool converged(false);

const Switch stopWhenConverged =
    presimControlPtr ? presimControlPtr->getOrDefault<Switch>("stopWhenConverged", false) : Switch(false);

const scalar fluxTol(presimControlPtr ? presimControlPtr->getOrDefault<double>("fluxTolerance", 1e-8) : 1e-8);
const scalar TTol(presimControlPtr ? presimControlPtr->getOrDefault<double>("TTolerance", 1e-8) : 1e-8);
