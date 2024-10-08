Mesh.MshFileVersion = 2.2;
// definition du pas du maillage
h = 0.05;
// d�finition des points (en 3D, raison pour laquelle il y a un 0 en z)
Point(1) = {0, 0, 0, h};
Point(2) = {2, 0, 0, h};
Point(3) = {2, 2, 0, h};
Point(4) = {0, 2, 0, h};
// d�finition des segments qui relient les points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// d�finition des contours ferm�s
Line Loop(1) = {1,2,3,4};
// d�finition des surfaces � partir contours ferm�s
Plane Surface(1) = {1};
// d�finition des �l�ments physiques : pour ces �l�ments, nous pourrons r�cup�rer
//									   les r�f�rences 
Physical Point(1) = {1,2,3,4};
Physical Line(1) = {1,2,3,4};
Physical Surface(1) = {1};
