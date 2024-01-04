lc = 0.01;

Point(1) = {0   ,0   ,0,lc};
Point(2) = {0   ,0.41,0,lc};
Point(3) = {2.2 ,0.41,0,lc};
Point(4) = {2.2 ,0   ,0,lc};

Point(5) = {0.15,0.2 ,0,lc};
Point(6) = {0.2 ,0.2 ,0,lc};
Point(7) = {0.25,0.2 ,0,lc};

Line(8)= {1,2};
Line(9)= {2,3};
Line(10)= {3,4};
Line(11)= {4,1};

Curve Loop(1) = {8,9,10,11};

Circle(12) = {5,6,7};
Circle(13) = {7,6,5};
Curve Loop(2) = {12,13};

Plane Surface(1) = {1,2};

Physical Line(1) = {8};
Physical Line(2) = {9};
Physical Line(3) = {10};
Physical Line(4) = {11};
Physical Curve(5) = {12};
Physical Curve(6) = {13};

Physical Surface(7)={1};

//Mesh.SaveAll=1;
Mesh 2;
Save "cilinder_2D_fine.msh";

Color Red{ Physical Curve{1}; }
Color Purple{ Physical Curve{2}; }
Color Green{ Physical Curve{3}; }
Color Purple{ Physical Curve{4}; }
Color Yellow{ Physical Curve{5}; }
Color Yellow{ Physical Curve{6}; }
