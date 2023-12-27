lc = 0.04;

Point(1) = {0   ,0   ,0,lc};
Point(2) = {0   ,0.41,0,lc};
Point(3) = {2.2 ,0.41,0,lc};
Point(4) = {2.2 ,0   ,0,lc};

Point(5) = {0.15,0.2 ,0,lc};
Point(6) = {0.2 ,0.2 ,0,lc};
Point(7) = {0.25,0.2 ,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Circle(5) = {5,6,7};
Circle(6) = {7,6,5};
Curve Loop(2) = {5,6};

Plane Surface(1) = {1,2};

Physical Curve(1) = {1};
Physical Curve(2) = {2};
Physical Curve(3) = {3};
Physical Curve(4) = {4};
Physical Curve(5) = {5};
Physical Curve(6) = {6};

//Physical Surface("Rectangle")={1};

Mesh.SaveAll=1;
Mesh 2;
Save "cilinder_2D_coarse.msh";

Color Red{ Physical Curve{1}; }
Color Purple{ Physical Curve{2}; }
Color Green{ Physical Curve{3}; }
Color Purple{ Physical Curve{4}; }
Color Yellow{ Physical Curve{5}; }
Color Yellow{ Physical Curve{6}; }
