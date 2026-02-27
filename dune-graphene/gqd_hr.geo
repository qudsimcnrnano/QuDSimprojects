// Gmsh project created on Tue Sep  2 16:21:24 2025
SetFactory("OpenCASCADE");
//+
Circle(1) = {1.9, -14, 75.68, 75.68, 0, 2*Pi};
//+
Circle(2) = {1.9, -14, 75.68, 75.68, 0, 2*Pi};
//+
Circle(3) = {-1.9, -14, 0, 151.36, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1};
//+
Surface(2) = {2};
//+
Curve Loop(4) = {3};
//+
Curve Loop(5) = {1};
//+
Plane Surface(3) = {4, 5};
//+
Curve Loop(6) = {3};
//+
Curve Loop(7) = {1};
//+
Surface(4) = {6, 7};
