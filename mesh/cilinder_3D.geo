m = 0.1; // mesh size
h = 0.41; // height in the z-direction

Point(1) = {0   ,0   ,0,m};   
Point(2) = {1.95,0   ,0,m};
Point(3) = {1.95,0.41,0,m}; 
Point(4) = {0   ,0.41,0,m};
Point(5) = {0.45,0.15,0,m};   
Point(6) = {0.55,0.15,0,m};
Point(7) = {0.55,0.25,0,m};   
Point(8) = {0.45,0.25,0,m};

Line(1)= {1, 2};  
Line(2)= {2, 3};  
Line(3)= {3, 4};  
Line(4)= {4, 1};  
Line(5)= {5, 6};
Line(6)= {6, 7}; 
Line(7)= {7, 8}; 
Line(8)= {8, 5};

Curve Loop(1) = {1,2,3,4};
Curve Loop(2) = {5,6,7,8};
Plane Surface(1) = {1, 2};

e() = Extrude {0, 0, h}{ Surface{1}; };

// Create physical groups, which are used to define the domain of the
// (co)homology computation and the subdomain of the relative (co)homology
// computation.

// Whole domain
Physical Volume(1) = {e(1)};


Mesh.SaveAll=1;
Mesh 3;
Save "cilinder_3D_coarse.msh";
Save "cilinder_3D.geo_unrolled";
