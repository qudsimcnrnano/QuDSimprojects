// Gmsh project created on tue sept 2 12:02:03 2025
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 75.68, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Surface(1) = {1};
//+
Circle(2) = {0, 0, 0, 151.36, 0, 2*Pi};
//+
Curve Loop(3) = {2};
//+
Curve Loop(4) = {1};
//+
Surface(2) = {3, 4};
//+
Curve Loop(5) = {2};
//+
Surface(2) = {5};
//+
Curve Loop(7) = {1};
//+
Curve Loop(8) = {2};
//+
Surface(3) = {7, 8};

//+
Curve Loop(9) = {1};
//+
Curve Loop(10) = {2};
//+
Curve Loop(11) = {1};
//+
Plane Surface(3) = {11};
//+
Curve Loop(12) = {1};
//+
Surface(4) = {12};
//+
Curve Loop(14) = {1};
//+
Curve Loop(15) = {2};
//+
Surface(5) = {14, 15};
//+
Curve Loop(16) = {2};
//+
Curve Loop(17) = {2};
//+
Surface(5) = {17};
//+
Curve Loop(19) = {2};
//+
Curve Loop(20) = {1};
//+
Plane Surface(6) = {19, 20};
//+
Curve Loop(22) = {1};
//+
Curve Loop(23) = {1};
//+
Curve Loop(24) = {1};
