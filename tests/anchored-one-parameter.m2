
print " This file illustrates the formulas in https://arxiv.org/abs/2512.18521"
quotientDiff = (x,f)->(
    --print f;
    x =sub(x,ring numerator f);
    (diff(x,numerator f)*(denominator f) - (numerator f)*(diff(x,denominator f)))/(denominator f)^2)


-- Define a generic method 
theoremEDRationalCurve = method()
theoremEDRationalCurve(List,Matrix,List) := (arrangement, fP1,stBoth) -> (
    nn := #arrangement;
    -- fP1 is the parametrized column vector giving the homogeneous coordinates in P^m
    if (numcols fP1) =!=1 then error" incorrect number of columns.";
    dimPPWorld := numrows fP1-1;
    firstC := first arrangement;
    dimPPImage := numrows firstC-1;
    camConf:=arrangement;
    ------------------------------------------------------------------------
    -- Setup: source ring and the rational normal curve map
    -- R = QQ[s,t] is the homogeneous coordinate ring of P^1
    ------------------------------------------------------------------------
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


oneParameterAnchoredMultiviewVarieties = method()
oneParameterAnchoredMultiviewVarieties(
    List,Matrix,List) := (
    arrangement,M,stBoth) -> (
    nn:=#arrangement;
    if (numcols M) >= (numrows M) then error" Take the transpose of the matrix M.";
    kGr := numcols M-1;
    dimPPWorld := numrows M-1;
    firstC := first arrangement;
    dimPPImage := numrows firstC-1;
    assert(numcols firstC==numrows M);
    wedgeCameraMatrix := (camera)->(
	rowIndex := subsets(numrows camera, kGr+1);
	colIndex := subsets(numcols camera, kGr+1);
	dMatrix := matrix for i in rowIndex list for j in colIndex list det submatrix(camera,i,j)
	);
    dArrangement := arrangement/(camera->wedgeCameraMatrix camera);
    s:=stBoth_0;
    t:=stBoth_1;    
    fP1 := matrix for i in subsets(dimPPWorld+1,kGr+1) list {det submatrix(M,i,)};
    theoremEDRationalCurve(dArrangement,fP1,stBoth); 
    )

end


restart 
needsPackage("Bertini", Configuration=>{"BERTINIexecutable"=>"/Applications/BertiniApple_v1.7/bertini"});
path=prepend("/Users/joserodriguez/Documents/GitHub/EuclideanDistanceDegree/",path);
needsPackage("EuclideanDistanceDegree");
load"anchored-one-parameter.m2"

---Input
--generateCameraMatrix=(dimPPImage,dimPPWorld)->matrix for i to dimPPImage list for j to dimPPWorld list random(1,100)
generateCameraMatrix=(dimPPImage,dimPPWorld)->matrix for i to dimPPImage list for j to dimPPWorld list if i==j or j==dimPPWorld then random(1,100) else 0


R=QQ[s,t]
EE =1
dimPPWorld = 3
stBoth= gens R
generateCameraMatrix(2,dimPPWorld)

kGr = 1
for dimPPImage from 3 to 3 do for nn from 1 to 2 do (
    arrangement := for i to nn-1 list generateCameraMatrix(dimPPImage,dimPPWorld);
    print arrangement;
    fP1 = matrix for i to dimPPWorld list for j to kGr list random({EE},R);
    print( "");
    print("f(P1)");
    print fP1;
    print ("Assuming n = "|nn|" and  h = "|dimPPImage|".");
    print ("If satysfying Theorem 2.3 then ED degree is "|toString(3*EE*nn-2)|".");
    print("Computing using probabilistic symbolic methods...");
    oneParameterAnchoredMultiviewVarieties(arrangement,fP1,stBoth)
    )


