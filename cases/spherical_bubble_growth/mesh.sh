#!/bin/bash
rm spherical.msh
gmsh spherical.geo -3
rm -r constant/polyMesh/ constant/*/polyMesh
gmshToFoam spherical.msh 
changeDictionary
rm -r constant/polyMesh/*Zones
rm -r constant/polyMesh/sets
mv constant/polyMesh constant/fluid/
