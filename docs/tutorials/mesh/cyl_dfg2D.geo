//Coarse
mesh_s1 = 0.02;
mesh_s2 = 0.035;
mesh_s3 = 0.005;

h1 = 0.2;
h2 = 0.21;
L  = 1.;

D  = 0.2;
r  = 0.05;


Point(1) = {0, 0,    0, mesh_s1};
Point(2) = {L, 0,    0, mesh_s1};


Point(3) = {L, h1+h2,    0, mesh_s1};
Point(4) = {0, h1+h2,    0, mesh_s1};



Point(5) = {D, D, 0, 1.0};
Point(6) = {D+r, D, 0, mesh_s3};
Point(7) = {D, D+r, 0, mesh_s3};
Point(8) = {D-r, D, 0, mesh_s3};
Point(9) = {D, D-r, 0, mesh_s3};




Circle(1) = {6, 5, 7};
Circle(2) = {7, 5, 8};
Circle(3) = {8, 5, 9};
Circle(4) = {9, 5, 6};


Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Line Loop(9) = {7, 8, 5, 6};
Line Loop(10) = {1, 2, 3, 4};

Plane Surface(11) = {9, 10};

Physical Line("inlet",  10)  = {8}; // right
Physical Line("outlet", 30)  = {6}; // left
Physical Line("walls",  20)  = {7, 5}; // top and bottom
Physical Line("walls",  20) += {1, 2, 3, 4}; // obstacle

Physical Surface("domain", 100) = {11};

// Mesh.MeshSizeFactor = 0.5;
