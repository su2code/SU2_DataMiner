//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.0008, 0, 0, 1.0};
//+
Point(3) = {0.0008, 0.002, 0, 1.0};
//+
Point(4) = {0.0004, 0.002, 0, 1.0};
//+
Point(5) = {0.0004, 0.003, 0, 1.0};
//+
Point(6) = {0.0008, 0.003, 0, 1.0};
//+
Point(7) = {0.0008, 0.008, 0, 1.0};
//+
Point(8) = {0, 0.008, 0, 1.0};
//+
Point(9) = {0.0004, 0, 0, 1.0};
//+
Point(10) = {0, 0.002, 0, 1.0};
//+
Point(11) = {0, 0.003, 0, 1.0};
//+
Point(12) = {0.0004, 0.008, 0, 1.0};
//+
Line(1) = {1, 9};
//+
Line(2) = {9, 2};
//+
Line(3) = {10, 4};
//+
Line(4) = {4, 3};
//+
Line(5) = {11, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {8, 12};
//+
Line(8) = {12, 7};
//+
Line(9) = {1, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 8};
//+
Line(12) = {9, 4};
//+
Line(13) = {4, 5};
//+
Line(14) = {5, 12};
//+
Line(15) = {2, 3};
//+
Line(16) = {6, 7};
//+
Curve Loop(1) = {15, -4, -12, 2};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 3, -12, -1};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, 5, -13, -3};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {11, 7, -14, -5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {16, -8, -14, 6};
//+
Plane Surface(5) = {5};
//+
Physical Curve("inlet", 17) = {1, 2};
//+
Physical Curve("outlet", 18) = {8, 7};
//+
Physical Curve("wall", 19) = {4, 13, 6};
//+
Physical Curve("sides", 20) = {15, 9, 10, 11, 16};
//+
Physical Surface("fluid", 21) = {1, 2, 3, 4, 5};
//+
Transfinite Curve {13} = 100 Using Progression 1;
//+
Transfinite Curve {13, 10} = 100 Using Progression 1;
//+
Transfinite Curve {1, 2, 3, 4, 5, 6, 7, 8} = 40 Using Progression 1;
//+
Transfinite Curve {9, 12, 15} = 200 Using Progression 1;
//+
Transfinite Curve {11, 14, 16} = 500 Using Progression 1;
//+
Transfinite Surface {1} = {9, 2, 3, 4};
//+
Transfinite Surface {2} = {1, 9, 4, 10};
//+
Transfinite Surface {3} = {10, 4, 5, 11};
//+
Transfinite Surface {4} = {11, 5, 12, 8};
//+
Transfinite Surface {5} = {5, 6, 7, 12};
//+
Recombine Surface {2, 1, 3, 4, 5};
