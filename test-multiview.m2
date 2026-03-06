
print " This file illustrates the formulas in https://arxiv.org/abs/2512.18521"
generateData = M -> for i in support M list random(1,100)
generateObjective = (data,I)-> (
    n:=#data;
    assert(#support I)==n;
    for i to n list random(1,100))

quotientDiff = (x,f)->((diff(x,numerator f)*(denominator f) - (numerator f)*(diff(x,denominator f)))/(denominator f)^2)


singleViewVarietiesOfCubicCurves = method()
singleViewVarietiesOfCubicCurves(ZZ) := (dimPPImage)->(
    nn:=1;
    dimPPWorld := 3;
    EE :=3;
    fst:=(s,t)->for i to dimPPWorld list s^i*t^(EE-i);
    R:=QQ[s,t];
    fP1 := transpose matrix {fst(R_0,R_1)};
    camConf := for i to nn-1 list generateCameraMatrix(dimPPImage,dimPPWorld);
    multprojectiveMultiview:=for camera in camConf list camera*fP1;
    affineMultiview:= matrix for i from 1 to dimPPImage list for view in multprojectiveMultiview list  view_(i,0)/view_(0,0);
    imageData :=  matrix for i from 0 to dimPPImage-1 list for i to nn-1 list random(1,100);
    --print affineMultiview;
    --print imageData;
    objectiveFunction := transpose matrix for viewIndex to nn-1 list for entryIndex to dimPPImage-1 list ((affineMultiview)_(entryIndex,viewIndex) - (imageData)_(entryIndex,viewIndex))^2;    
    psi := sum flatten entries  objectiveFunction;
    -- QQ[t]
    -- instance(ring (t), FractionField) -- false
    -- instance(ring (1/t), FractionField)-- true
    idealCriticalPoints := if instance(ring psi, FractionField) then saturate(ideal for i in {0,1} list numerator quotientDiff(R_i, psi), denominator psi) else ideal diff(R_0,psi);
    degreeIdealCriticalPoints = degree idealCriticalPoints;
    print ("Degree of ideal of parametric critical points "|degreeIdealCriticalPoints);
    pp := symbol pp;
    S:=QQ[gens R|flatten transpose for entryIndex to dimPPImage-1 list for viewIndex to nn-1 list  pp_(entryIndex,viewIndex)];
    --print flatten transpose for entryIndex to dimPPImage-1 list for viewIndex to nn-1 list pp_(entryIndex,viewIndex);
    imageVars:=genericMatrix(S,S_2,dimPPImage,nn);
    --print imageVars;
    ratGraph:=ideal (imageVars - sub((affineMultiview),frac S));
    graph := ideal for i in ratGraph_* list numerator i;
    scan(ratGraph_*, f-> graph = saturate(graph,denominator f));
    modelImage := eliminate({S_0,S_1},saturate(graph,ideal(S_0,S_1)));
    print ("Degree of image "|toString(degree modelImage));
    modelEDDegree := determinantalUnitEDDegree flatten entries gens sub(modelImage, QQ[support modelImage]);
    print ("ED degree "|toString(modelEDDegree));    
    )
    
theoremEDRationalCurve = method()
theoremEDRationalCurve(ZZ,ZZ,FunctionClosure) := (nn,dimPPImage,fst)->(
    R:=QQ[s,t];
    dimPPWorld=#fst(s,t)-1;
    fP1 := transpose matrix {fst(R_0,R_1)};
    camConf := for i to nn-1 list generateCameraMatrix(dimPPImage,dimPPWorld);
    multprojectiveMultiview:=for camera in camConf list camera*fP1;
    affineMultiview:= matrix for i from 1 to dimPPImage list for view in multprojectiveMultiview list  view_(i,0)/view_(0,0);
    imageData :=  matrix for i from 0 to dimPPImage-1 list for i to nn-1 list random(1,100);
    print affineMultiview;
    print imageData;
    objectiveFunction := transpose matrix for viewIndex to nn-1 list for entryIndex to dimPPImage-1 list ((affineMultiview)_(entryIndex,viewIndex) - (imageData)_(entryIndex,viewIndex))^2;    
    psi := sum flatten entries  objectiveFunction;
    -- QQ[t]
    -- instance(ring (t), FractionField) -- false
    -- instance(ring (1/t), FractionField)-- true
    idealCriticalPoints := if instance(ring psi, FractionField) then saturate(ideal for i in {0,1} list numerator quotientDiff(R_i, psi), denominator psi) else ideal diff(R_0,psi);
    degreeIdealCriticalPoints = degree idealCriticalPoints;
    print degreeIdealCriticalPoints;
    pp := symbol pp;
    S:=QQ[gens R|flatten transpose for entryIndex to dimPPImage-1 list for viewIndex to nn-1 list  pp_(entryIndex,viewIndex)];
    print flatten transpose for entryIndex to dimPPImage-1 list for viewIndex to nn-1 list pp_(entryIndex,viewIndex);
    imageVars:=genericMatrix(S,S_2,dimPPImage,nn);
    print imageVars;
    ratGraph:=ideal (imageVars - sub((affineMultiview),frac S));
    graph := ideal for i in ratGraph_* list numerator i;
    scan(ratGraph_*, f-> graph = saturate(graph,denominator f));
    modelImage := eliminate({S_0,S_1},saturate(graph,ideal(S_0,S_1)));
    print ("Degree of image "|toString(degree modelImage));
    if nn==1 or nn==2 then (
	modelEDDegree := determinantalUnitEDDegree flatten entries gens sub(modelImage, QQ[support modelImage]);
	print ("ED degree "|toString(modelEDDegree));
	);
    )




end
--end
 restart
needsPackage("Bertini", Configuration=>{"BERTINIexecutable"=>"/Applications/BertiniApple_v1.7/bertini"});
path=prepend("/Users/joserodriguez/Documents/GitHub/EuclideanDistanceDegree/",path);
needsPackage("EuclideanDistanceDegree");
load"test-multiview.m2"

---Input
--generateCameraMatrix=(dimPPImage,dimPPWorld)->matrix for i to dimPPImage list for j to dimPPWorld list random(1,100)
generateCameraMatrix=(dimPPImage,dimPPWorld)->matrix for i to dimPPImage list for j to dimPPWorld list if i==j or j==dimPPWorld then random(1,100) else 0

for dimPPImage from 2 to 5 do print toString (dimPPImage=>singleViewVarietiesOfCubicCurves(dimPPImage))

EE =1
fst= (s,t)->for i to EE list s^i*t^(EE-i)

for dimPPImage from 2 to 3 do for nn from 1 to 2 do (
    print( "");
    print ("n = "|nn|"h = "|dimPPImage);
    print ("If satysfying Theorem 2.3 then ED degree is "|toString(3*EE*nn-2));
    print("Computing using probabilistic symbolic methods");
    theoremEDRationalCurve(nn,dimPPImage,fst)
    )
