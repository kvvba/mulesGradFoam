length = 0.5e-3;
extend = length + 0.025e-3;
// length = 1.2e-3;
// extend = length + 0.1e-3;
// length = 2e-3;
// extend = length + 0.1e-3;
r = (extend^2 * 2)^0.5;
diff = r - length;
resolution = 1.0e-6;
n = length/resolution;
// resolution = length/n;
nout = diff / (resolution);

Point(0) = {0, 0, 0};
Point(1) = {length, 0, 0};
Point(2) = {length, 0, length};
Point(3) = {0, 0, length};

Point(10) = {r, 0, 0};
Point(11) = {0, 0, r};
Point(12) = {extend, 0, extend};

Rotate { {0, 0, 1}, {0, 0, 0}, -(2.5/360)*2*Pi}{
	Point{1, 2, 3, 10, 11, 12};
}

Line(1) = {0, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 0};


Circle(10) = {10, 0, 12};
Circle(11) = {12, 0, 11};

Line(12) = {1, 10};
Line(13) = {11, 3};
Line(14) = {2, 12};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {12, 10, -14, -2};
Curve Loop(3) = {-3, 14, 11, 13};

Surface(1) = {1};
Surface(2) = {2};
Surface(3) = {3};

Transfinite Line {1, -3} = n Using Progression 1;
Transfinite Line {4, -2} = n Using Progression 1;
Transfinite Surface {1};

Transfinite Line {2, 10} = n Using Progression 1;
Transfinite Line {12, 14} = nout Using Progression 1;
Transfinite Surface {2};

Transfinite Line {-13, 14} = nout Using Progression 1;
Transfinite Line {-3, 11} = n Using Progression 1;
Transfinite Surface {3};
Recombine Surface "*";

Extrude { {0, 0, 1}, {0, 0, 0}, (5/360)*2*Pi} {
	Surface{1, 2, 3};
	Layers {{1}, {1}};
	Recombine;
}

Mesh 3;

Physical Volume(75) = {3, 2, 1};
Physical Surface("front", 71) = {3, 1, 2};
Physical Surface("back", 72) = {53, 31, 70};
Physical Surface("edge", 73) = {44, 68};
Physical Surface("bottom", 74) = {22, 40};
Physical Curve("axis", 75) = {13, 4};

Save "spherical.msh";
