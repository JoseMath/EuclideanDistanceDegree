restart
loadPackage "EuclideanDistanceDegree"

--INPUTS: 
--F 
----A list of polynomials that generate the ideal of X (assume X is generically reduced)

R = QQ[I][x0,x1,x2]; --I plays the role of the imaginary unit in Bertini. 
F = {x0^2*x2 - x1^2*(x1 + x2)};  --ED degree is 7
UED = leftKernelUnitEDDegree F

F = {x0^2*x1 -(x1 - I*x2)^2*x2};  --ED degree is 7
UED = leftKernelUnitEDDegree F

F = {x0^3 - (I*x0^2 + x1^2)*x2};  --ED degree is 6
UED = leftKernelUnitEDDegree F