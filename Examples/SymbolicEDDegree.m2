restart
loadPackage "EuclideanDistanceDegree"

--In each example the variety is defined by F=0. 
--The input F is a system of polynomials that define a generically reduced variety X. 
--If R is a polynomial ring with a coefficient ring QQ or ZZ/p, then the output is the ED degree of X. 

--Example 4.1.
R = QQ[x0,x1,x2];
F = {x0^2*x2 - x1^2*(x1 + x2)};
determinantalUnitEDDegree F  --ED degree is 7

--Example 4.2.
R = QQ[jj]/ideal(jj^2+1)[x0,x1,x2];
F = {x0^2*x1 -(x1 - jj*x2)^2*x2};
--Here the output has a factor of 2 to account for the imaginary unit jj.
determinantalUnitEDDegree F  --ED degree is 7

--Example 4.3.
R = QQ[jj]/ideal(jj^2+1)[x0,x1,x2];
F = {x0^3 - (jj*x0^2 + x1^2)*x2};
--Here the output has a factor of 2 to account for the imaginary unit jj.
determinantalUnitEDDegree F  --ED degree is 6

--Example 4.4.
R = QQ[x0,x1,x2,x3];
F = {x0^2*x1 - x2*x3^2};
determinantalUnitEDDegree F  --ED degree is 10