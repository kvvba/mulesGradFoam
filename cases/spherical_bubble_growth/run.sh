#!/bin/bash
rm bubble.csv
decomposePar -allRegions
mpirun -n 12 mulesGradFoam -parallel > log &
