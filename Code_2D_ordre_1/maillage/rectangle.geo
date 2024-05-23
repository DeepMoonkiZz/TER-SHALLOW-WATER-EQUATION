// À adapter en fonction du raffinement de maillage souhaité
h = 0.06;

Point(1) = {0,0,0,h};
Point(2) = {0,1,0,h};
Point(3) = {2,0,0,h};
Point(4) = {2,1,0,h};
Line(1) = {2,1};
Line(2) = {1,3};
Line(3) = {3,4};
Line(4) = {4,2};
Line Loop(1) = {1:4};
Plane Surface(1) = {1};

Physical Curve("entree", 10) = {3};
Physical Curve("haut", 11) = {4};
Physical Curve("sortie", 12) = {1};
Physical Curve("bas", 13) = {2};

Physical Surface("droite", 16) = {1};
