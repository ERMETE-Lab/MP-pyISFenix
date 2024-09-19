SetFactory("OpenCASCADE");
// Mesh.MeshSizeFactor = 0.1;

Lc1 = 1;
Lc2 = 0.4;

radius = 0.25;
side_length = .75;

channel_length = 5;

lateral_rod_radius = 0.15;

end_point_rod = channel_length * 0.4;
offset_x = 0.05;
offset_y = 0.03;

///////// LATERAL BOUNDARIES ///////////
Point(1) = {0, 0, 0, Lc1};
Point(2) = { side_length/2,  side_length/2, 0, Lc1};
Point(3) = { side_length/2, -side_length/2, 0, Lc1};
Point(4) = {-side_length/2, -side_length/2, 0, Lc1};
Point(5) = {-side_length/2,  side_length/2, 0, Lc1};

Point(6) = {offset_x, offset_y, channel_length, Lc1};
Point(7)  = { side_length/2,  side_length/2, channel_length, Lc1};
Point(8)  = { side_length/2, -side_length/2, channel_length, Lc1};
Point(9)  = {-side_length/2, -side_length/2, channel_length, Lc1};
Point(10) = {-side_length/2,  side_length/2, channel_length, Lc1};

///////// LATERAL NUCLEAR RODS ///////////

Point(11) = { side_length/2,  side_length/2 - lateral_rod_radius, 0, Lc1};
Point(12) = { side_length/2 - lateral_rod_radius,  side_length/2, 0, Lc1};

Point(13) = { -side_length/2,  side_length/2 - lateral_rod_radius, 0, Lc1};
Point(14) = { -side_length/2 + lateral_rod_radius,  side_length/2, 0, Lc1};

Point(15) = { -side_length/2,  - side_length/2 + lateral_rod_radius, 0, Lc1};
Point(16) = { -side_length/2 + lateral_rod_radius,  - side_length/2, 0, Lc1};

Point(17) = { side_length/2,  - side_length/2 + lateral_rod_radius, 0, Lc1};
Point(18) = { side_length/2 - lateral_rod_radius,  - side_length/2, 0, Lc1};

Point(19) = { side_length/2,  side_length/2 - lateral_rod_radius, channel_length, Lc1};
Point(20) = { side_length/2 - lateral_rod_radius,  side_length/2, channel_length, Lc1};

Point(21) = { -side_length/2,  side_length/2 - lateral_rod_radius, channel_length, Lc1};
Point(22) = { -side_length/2 + lateral_rod_radius,  side_length/2, channel_length, Lc1};

Point(23) = { -side_length/2,  - side_length/2 + lateral_rod_radius, channel_length, Lc1};
Point(24) = { -side_length/2 + lateral_rod_radius,  - side_length/2, channel_length, Lc1};

Point(25) = { side_length/2,  - side_length/2 + lateral_rod_radius, channel_length, Lc1};
Point(26) = { side_length/2 - lateral_rod_radius,  - side_length/2, channel_length, Lc1};

Circle(1) = {14, 5, 13};
Circle(2) = {15, 4, 16};
Circle(3) = {17, 3, 18};
Circle(4) = {12, 2, 11};
Circle(5) = {22, 10, 21};
Circle(6) = {24, 9, 23};
Circle(7) = {25, 8, 26};
Circle(8) = {20, 7, 19};

Line(9) = {17, 11};
Line(10) = {12, 14};
Line(11) = {13, 15};
Line(12) = {16, 18};
Line(13) = {25, 19};
Line(14) = {20, 22};
Line(15) = {21, 23};
Line(16) = {24, 26};
Line(17) = {25, 17};
Line(18) = {26, 18};
Line(19) = {24, 16};
Line(20) = {23, 15};
Line(21) = {19, 11};
Line(22) = {20, 12};
Line(23) = {22, 14};
Line(24) = {21, 13};

//////////////////////////// CIRCLES ////////////////////////////////////

Point(27) = {offset_x, offset_y, end_point_rod, Lc2};
Point(28) = {radius+offset_x, offset_y, end_point_rod, Lc2};
Point(29) = {offset_x, radius+offset_y, end_point_rod, Lc2};
Point(30) = {-radius+offset_x, offset_y, end_point_rod, Lc2};
Point(31) = {offset_x, -radius+offset_y, end_point_rod, Lc2};

Point(32) = {offset_x, offset_y, channel_length, Lc2};
Point(33) = {radius+offset_x, offset_y, channel_length, Lc2};
Point(34) = {offset_x, radius+offset_y, channel_length, Lc2};
Point(35) = {-radius+offset_x, offset_y, channel_length, Lc2};
Point(36) = {offset_x, -radius+offset_y, channel_length, Lc2};

Circle(25) = {28, 27, 29};
Circle(26) = {29, 27, 30};
Circle(27) = {30, 27, 31};
Circle(28) = {28, 27, 31};
Circle(29) = {33, 6, 34};
Circle(30) = {35, 6, 34};
Circle(31) = {36, 6, 35};
Circle(32) = {33, 6, 36};

Line(33) = {28, 33};
Line(34) = {29, 34};
Line(35) = {30, 35};
Line(36) = {31, 36};

//////////////////////////// SURFACES ////////////////////////////////////

Curve Loop(1) = {23, -10, -22, 14};
Curve Loop(2) = {13, 21, -9, -17};
Curve Loop(3) = {18, -12, -19, 16};
Curve Loop(4) = {20, -11, -24, 15};
Curve Loop(5) = {11, 2, 12, -3, 9, -4, 10, 1};
Curve Loop(6) = {14, 5, 15, -6, 16, -7, 13, -8};
Curve Loop(7) = {29, -30, -31, -32};
Curve Loop(8) = {20, 2, -19, 6};
Curve Loop(10) = {3, -18, -7, 17};
Curve Loop(12) = {4, -21, -8, 22};
Curve Loop(14) = {5, 24, -1, -23};
Curve Loop(16) = {34, -29, -33, 25};
Curve Loop(18) = {34, -30, -35, -26};
Curve Loop(20) = {27, 36, 31, -35};
Curve Loop(22) = {32, -36, -28, 33};
Curve Loop(24) = {26, 27, -28, 25};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6, 7};

Surface(7) = {8};
Surface(8) = {10};
Surface(9) = {12};
Surface(10) = {14};
Surface(11) = {16};
Surface(12) = {18};
Surface(13) = {20};
Surface(14) = {22};

Plane Surface(15) = {24};

///////////////////////////// VOLUMES /////////////////////////////////////

Surface Loop(1) = {1, 10, 5, 4, 7, 6, 3, 8, 2, 9, 14, 11, 15, 12, 13};
Volume(1) = {1};

///////////////////////////// PHYSICAL GROUPS /////////////////////////////////////

Physical Volume("domain", 1000) = {1};
Physical Surface("inlet", 100) = {5};
Physical Surface("outlet", 110) = {6};
Physical Surface("lateral_walls", 120) = {2, 1, 4, 3};
Physical Surface("control_rod", 130) = {11, 14, 13, 12, 15};

// The lateral rods are labelled from the top view (watching the exit)
Physical Surface("top_right_nuc", 140) = {9};
Physical Surface("top_left_nuc", 150) = {10};
Physical Surface("bottom_left_nuc", 160) = {7};
Physical Surface("bottom_right_nuc", 170) = {8};