
print " This file illustrates the formulas in https://arxiv.org/abs/2512.18521"
generateData = M -> for i in support M list random(1,100)
generateObjective = (data,I)-> (
    n:=#data;
    assert(#support I)==n;
    for i to n list random(1,100))

quotientDiff = (x,f)->((diff(x,numerator f)*(denominator f) - (numerator f)*(diff(x,denominator f)))/(denominator f)^2)


singleViewVarietiesOfCubicCurves = method()
singleViewVarietiesOfCubicCurves(ZZ) := (dimPPImage)->(
    generateCameraMatrix:=(dimPPImage,dimPPWorld)->matrix for i to dimPPImage list for j to dimPPWorld list if i==j or j==dimPPWorld then random(1,100) else 0;
    singleViewVarietiesOfCubicCurves(dimPPImage,generateCameraMatrix)
    )
    
singleViewVarietiesOfCubicCurves(ZZ,FunctionClosure) := (dimPPImage,generateCameraMatrix)->(
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

-*
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
*-

-- Define a generic method symbol
theoremEDRationalCurve = method()
-- Specialize the method for arguments: (nn : ZZ, dimPPImage : ZZ, fst : FunctionClosure)
theoremEDRationalCurve(ZZ,ZZ,FunctionClosure) := (nn, dimPPImage, fst) -> (
    generateCameraMatrix:=(dimPPImage,dimPPWorld)->matrix for i to dimPPImage list for j to dimPPWorld list if i==j or j==dimPPWorld then random(1,100) else 0;
    theoremEDRationalCurve(nn, dimPPImage, fst,generateCameraMatrix) 
    )

theoremEDRationalCurve(ZZ,ZZ,FunctionClosure,FunctionClosure) := (nn, dimPPImage, fst,generateCameraMatrix) -> (
    ------------------------------------------------------------------------
    -- Setup: source ring and the rational normal curve map
    -- R = QQ[s,t] is the homogeneous coordinate ring of P^1
    ------------------------------------------------------------------------
    R := QQ[s,t];

    -- The ambient projective space for the world/model is determined by fst(s,t)
    -- (fst should be a parametrization P^1 -> P^m; the image dimension is #coords - 1)
    dimPPWorld := #fst(s,t) - 1;

    -- fP1 is the parametrized column vector giving the homogeneous coordinates in P^m
    fP1 := transpose matrix { fst(R_0, R_1) };

    ------------------------------------------------------------------------
    -- Camera configuration: create nn cameras from P^m -> P^{dimPPImage}
    ------------------------------------------------------------------------
    -- generateCameraMatrix(dimPPImage, dimPPWorld) is provided by user;
    -- theoremEDRationalCurve(ZZ,ZZ,FunctionClosure) =  (nn, dimPPImage, fst)->output
    ---- automatically creates generateCameraMatrix=(dimPPImage,dimPPWorld)->matrix for i to dimPPImage list for j to dimPPWorld list if i==j or j==dimPPWorld then random(1,100) else 0

    -- camConf is a list of nn random/projective cameras
    camConf := for i to nn-1 list generateCameraMatrix(dimPPImage, dimPPWorld);

    -- Apply each camera to the parametric curve (multiview projections in homogeneous coords)
    multprojectiveMultiview := for camera in camConf list camera * fP1;

    -- Convert each projected view to affine coordinates by dividing by the first row (row 0).
    -- This yields a dimPPImage x nn matrix of rational functions in s,t.
    -- Rows 1..dimPPImage are taken as affine coordinates (row 0 is the scale).
    affineMultiview := matrix for i from 1 to dimPPImage list for view in multprojectiveMultiview list view_(i,0) / view_(0,0);

    ------------------------------------------------------------------------
    -- Synthetic image data: random integers to simulate measured points
    ------------------------------------------------------------------------
    -- imageData has the same shape as affineMultiview: dimPPImage x nn
    -- (uses integers in [1..100])
    imageData := matrix for i from 0 to dimPPImage-1 list for i to nn-1 list random(1,100);

    -- Diagnostics
    --print affineMultiview;
    --print imageData;

    ------------------------------------------------------------------------
    -- Least-squares objective: sum of squared residuals across all entries
    -- psi(s,t) = sum_{views,coords} (affineMultiview - imageData)^2
    ------------------------------------------------------------------------
    objectiveFunction := transpose matrix for viewIndex to nn-1 list for entryIndex to dimPPImage-1 list (
        (affineMultiview)_(entryIndex, viewIndex) - (imageData)_(entryIndex, viewIndex)
    )^2;

    psi := sum flatten entries objectiveFunction;

    ------------------------------------------------------------------------
    -- Critical equations: ∂psi/∂s = 0 and ∂psi/∂t = 0
    -- If psi lives in a fraction field, we use quotientDiff and clear denominators
    -- to get polynomial equations. Otherwise, we differentiate directly.
    ------------------------------------------------------------------------
    -- QQ[t]
    -- instance(ring (t), FractionField) -- false
    -- instance(ring (1/t), FractionField) -- true
    idealCriticalPoints := if instance(ring psi, FractionField)
        then saturate(
                ideal for i in {0,1} list numerator quotientDiff(R_i, psi),
                denominator psi)
        else ideal diff(R_0, psi);

    -- Degree of the (zero-dimensional) critical locus (if finite)
    degreeIdealCriticalPoints = degree idealCriticalPoints;
    print "Working parametrically we find: ";
    print ("    the ED degree is "|degreeIdealCriticalPoints);

    ------------------------------------------------------------------------
    -- Build a larger ring S = QQ[s,t, {pp_(i,j)}], where pp_(i,j) are variables
    -- representing image coordinates for entryIndex=0..dimPPImage-1 and viewIndex=0..nn-1.
    ------------------------------------------------------------------------
    print"Working implicitly we find: ";
    pp := symbol pp;

    -- Adjoin the pp_(i,j) to gens R (which are s,t) to make S
    S := QQ[ gens R | flatten transpose for entryIndex to dimPPImage-1 list for viewIndex to nn-1 list pp_(entryIndex, viewIndex) ];

    -- Show the variables {pp_(i,j)}
    --print flatten transpose for entryIndex to dimPPImage-1 list for viewIndex to nn-1 list pp_(entryIndex, viewIndex);

    -- Create a dimPPImage x nn generic matrix of pp-variables over S
    -- NOTE: The use of S_2 here indexes a generator, not a symbol; see the optional fix below.
    imageVars := genericMatrix(S, S_2, dimPPImage, nn);
    --print imageVars;

    ------------------------------------------------------------------------
    -- Graph of the rational map (s,t) ↦ affineMultiview, with coordinates pp_(i,j):
    -- ratGraph ⊂ Spec(S): equations pp_(i,j) - affineMultiview_(i,j) = 0
    -- We first set them over Frac(S), then clear denominators to obtain a polynomial ideal.
    ------------------------------------------------------------------------
    ratGraph := ideal (imageVars - sub((affineMultiview), frac S));

    -- Clear denominators: take numerators and saturate by each denominator factor
    graph := ideal for i in ratGraph_* list numerator i;
    scan(ratGraph_*, f -> graph = saturate(graph, denominator f));

    ------------------------------------------------------------------------
    -- Eliminate (s,t) to get the model of the image in the pp-variables.
    -- Here S_0,S_1 correspond to s,t because S was built as QQ[s,t,pp...].
    ------------------------------------------------------------------------
    modelImage := eliminate({S_0, S_1}, saturate(graph, ideal(S_0, S_1)));

    print ("    the degree of the affine multiview variety is " | toString(degree modelImage));
    modelEDDegree := determinantalUnitEDDegree flatten entries gens sub(modelImage, QQ[support modelImage]);
    print ("    the ED degree is " | toString(modelEDDegree));
)

end
--end
 restart
--needsPackage("Bertini", Configuration=>{"BERTINIexecutable"=>"/Applications/BertiniApple_v1.7/bertini"});
path=prepend("/Users/joserodriguez/Documents/GitHub/EuclideanDistanceDegree/",path);
needsPackage("EuclideanDistanceDegree");
load"test-multiview.m2"

---Input
--generateCameraMatrix=(dimPPImage,dimPPWorld)->matrix for i to dimPPImage list for j to dimPPWorld list random(1,100)
generateCameraMatrix=(dimPPImage,dimPPWorld)->matrix for i to dimPPImage list for j to dimPPWorld list if i==j or j==dimPPWorld then random(1,100) else 0

--for dimPPImage from 2 to 5 do print toString (dimPPImage=>singleViewVarietiesOfCubicCurves(dimPPImage))

EE =1
fst= (s,t)->for i to EE list s^i*t^(EE-i)

fst= (s,t)->for i to EE list s^i*t^(EE-i)
dimPPWorld = 3

fst = (s,t)->for i to dimPPWorld+1 list sum for i to EE list random(1,100)*s^i*t^(EE-i)
-- testR=QQ[A,B]; fst(A,B)

for dimPPImage from 2 to 3 do for nn from 1 to 2 do (
    print( "");
    print ("Assuming n = "|nn|" and  h = "|dimPPImage|".");
    print ("If satysfying Theorem 2.3 then ED degree is "|toString(3*EE*nn-2)|".");
    print("Computing using probabilistic symbolic methods...");
    theoremEDRationalCurve(nn,dimPPImage,fst)
)