SetFactory("OpenCASCADE");
mesh_s1 = 0.75;
mesh_s2 = 0.25;

L = 1.;
H = 0.25;

Point(1) = {0, 0, 0, mesh_s1};
Point(2) = {L, 0, 0, mesh_s1};
Point(3) = {L, H, 0, mesh_s1};
Point(4) = {0, H, 0, mesh_s1};

radius = 0.05;

x1 = 0.25;
x2 = 0.5;
x3 = 0.75;

Point(5)  = {x1 - radius, 0, 0, mesh_s2};
Point(6)  = {x1 + radius, 0, 0, mesh_s2};
Point(7)  = {x2 - radius, 0, 0, mesh_s2};
Point(8)  = {x2 + radius, 0, 0, mesh_s2};
Point(9)  = {x3 - radius, 0, 0, mesh_s2};
Point(10) = {x3 + radius, 0, 0, mesh_s2};

Point(11) = {x1, 0, 0, mesh_s1};
Point(12) = {x2, 0, 0, mesh_s1};
Point(13) = {x3, 0, 0, mesh_s1};

Point(14) = {x1 - radius, H, 0, mesh_s2};
Point(15) = {x1 + radius, H, 0, mesh_s2};
Point(16) = {x2 - radius, H, 0, mesh_s2};
Point(17) = {x2 + radius, H, 0, mesh_s2};
Point(18) = {x3 - radius, H, 0, mesh_s2};
Point(19) = {x3 + radius, H, 0, mesh_s2};

Point(20) = {x1, H, 0, mesh_s1};
Point(21) = {x2, H, 0, mesh_s1};
Point(22) = {x3, H, 0, mesh_s1};

Point(23) = {(x1+x2)/2, H/2, 0, mesh_s1};
Point(24) = {(x2+x3)/2, H/2, 0, mesh_s1};

Line(1) = {1, 5};
Line(2) = {6, 7};
Line(3) = {8, 9};
Line(4) = {10, 2};
Line(5) = {2, 3};
Line(6) = {3, 19};
Line(7) = {18, 17};
Line(8) = {16, 15};
Line(9) = {14, 4};
Line(10) = {4, 1};

Circle(11) = {5, 11, 6};
Circle(12) = {7, 12, 8};
Circle(13) = {9, 13, 10};
Circle(14) = {15, 20, 14};
Circle(15) = {17, 21, 16};
Circle(16) = {19, 22, 18};

Circle(17) = {(x1+x2)/2, H/2, 0, radius, 0, 2*Pi};
Circle(18) = {(x3+x2)/2, H/2, 0, radius, 0, 2*Pi};

Curve Loop(1) = {9, 10, 1, 11, 2, 12, 3, 13, 4, 5, 6, 16, 7, 15, 8, 14};
Curve Loop(2) = {17};
Curve Loop(3) = {18};

Plane Surface(1) = {1, 2, 3};

Physical Curve("inlet", 10) = {10};
Physical Curve("adiab_walls", 20) = {1, 2, 3, 4, 6, 7, 8, 9};
Physical Curve("coldest_walls", 30) = {14, 12, 16};
Physical Curve("hottest_walls", 40) = {11, 15, 13};
Physical Curve("hot_centre", 50) = {17};
Physical Curve("cold_centre", 60) = {18};
Physical Curve("outlet", 70) = {5};

Physical Surface("domain", 100) = {1};

// Mesh.MeshSizeFactor = 0.015;
