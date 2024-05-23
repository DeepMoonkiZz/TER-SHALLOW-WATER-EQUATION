// À adapter en fonction du raffinement de maillage souhaité
h = 0.06;

Point(1) = {0,0,0,h};
Point(2) = {0,1,0,h};
Point(3) = {1,0,0,h};
Point(4) = {1,1,0,h};
Line(1) = {2,1};
Line(2) = {1,3};
Line(3) = {3,4};
Line(4) = {4,2};
Line Loop(1) = {1:4};
Plane Surface(1) = {1};

Point(5) = {-1,0,0,h};
Point(6) = {-1,1,0,h};
Line(5) = {1,5};
Line(6) = {5,6};
Line(7) = {6,2};

Line Loop(2) = {5,6,7,1};
Plane Surface(2) = {2};


Physical Curve("entree", 10) = {3};
Physical Curve("haut", 11) = {4,7};
Physical Curve("sortie", 12) = {6};
Physical Curve("bas", 13) = {5,2};
Physical Curve("barrage", 14) = {1};

Physical Surface("gauche", 15) = {2};
Physical Surface("droite", 16) = {1};
