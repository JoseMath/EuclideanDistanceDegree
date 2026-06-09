restart
loadPackage("EuclideanDistanceDegree",Reload=>true)
help EuclideanDistanceDegree

examples EuclideanDistanceDegree

R = QQ[x,y];
F = {x^2+y^2-1};
2 == determinantalUnitEDDegree F
2 == leftKernelUnitEDDegree F
4 == determinantalGenericEDDegree F
4 == leftKernelGenericEDDegree F

genericWeightVector = {2,3};
unitWeightVector = {1,1};
dataVector = {5,7};
4 == symbolicWeightEDDegree(F, dataVector, genericWeightVector)
2 == symbolicWeightEDDegree(F, dataVector, unitWeightVector)

R = QQ[x1,x2,x3,x4];
F = {det genericMatrix(R,2,2)};
P = (F,F);
6 == numericGenericEDDegree P
2 == numericUnitEDDegree P

R = QQ[x1,x2,x3,x4,x5,x6]
F = (minors(2, genericMatrix(R,3,2)))_*;
G = drop(F,-1);
P = (F,G);
#G == codim ideal F
10 == numericGenericEDDegree P
2 == numericUnitEDDegree P

dir = temporaryFileName();
if not fileExists dir then mkdir dir;

NCO = newNumericalComputationOptions(P, TempDirectory => dir);
NCO#"TargetWeight" = apply(#gens R, i->1);
homotopyEDDegree(NCO,"Weight",true,true)  -- recovers the unit ED degree

NCO#"TargetWeight" = (apply(-1+#gens R, i->1)) | {0};
homotopyEDDegree(NCO,"Weight",false,true)

NCO = newNumericalComputationOptions(P, TempDirectory => dir);
NCO#"TargetData" = apply(#gens R, i->1);
homotopyEDDegree(NCO,"Data",true,true)

NCO#"TargetWeight"=(apply(-1+#gens R, i->1)) | {0};
homotopyEDDegree(NCO,"Data",false,true)