-- Extend theorem 2.3 from https://arxiv.org/abs/2512.18521 to rational surfaces
needsPackage "EuclideanDistanceDegree"

generateCameraMatrix = method()
generateCameraMatrix(ZZ, ZZ) := (dimPPImage, dimPPWorld) -> (
    matrix for i to dimPPImage list for j to dimPPWorld list if i==j or j==dimPPWorld then random(1,100) else 0
)

theoremEDRationalCurve = method()
theoremEDRationalCurve(ZZ, ZZ, FunctionClosure) := (nn, dimPPImage, fst) -> (
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
    camConf := for i to nn-1 list generateCameraMatrix(dimPPImage, dimPPWorld);

    -- Apply each camera to the parametric curve (multiview projections in homogeneous coords)
    multprojectiveMultiview := for camera in camConf list camera * fP1;

    -- Convert each projected view to affine coordinates by dividing by the first row (row 0).
    -- This yields a dimPPImage x nn matrix of rational functions in s,t.
    -- Rows 1..dimPPImage are taken as affine coordinates (row 0 is the scale).
    affineMultiview := matrix for i from 1 to dimPPImage list for view in multprojectiveMultiview list view_(i,0) / view_(0,0);

    ------------------------------------------------------------------------
    -- Build a larger ring S = QQ[s,t, {pp_(i,j)}], where pp_(i,j) are variables
    -- representing image coordinates for entryIndex=0..dimPPImage-1 and viewIndex=0..nn-1.
    ------------------------------------------------------------------------
    pp := symbol pp;

    -- Adjoin the pp_(i,j) to gens R (which are s,t) to make S
    S := QQ[ gens R | flatten transpose for entryIndex to dimPPImage-1 list for viewIndex to nn-1 list pp_(entryIndex, viewIndex) ];

    -- Show the variables {pp_(i,j)}
    --print flatten transpose for entryIndex to dimPPImage-1 list for viewIndex to nn-1 list pp_(entryIndex, viewIndex);

    -- Create a dimPPImage x nn generic matrix of pp-variables over S
    -- NOTE: The use of S_2 here indexes a generator, not a symbol
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

    print ("    The degree of the affine multiview variety is " | toString(degree modelImage));
    modelEDDegree := determinantalUnitEDDegree flatten entries gens sub(modelImage, QQ[support modelImage]);
    print ("    The ED degree is " | toString(modelEDDegree));
)

theoremEDRationalSurface = method()
theoremEDRationalSurface(ZZ, ZZ, FunctionClosure) := (nn, dimPPImage, f) -> (
    R := QQ[s,t,u];
    dimPPWorld := #f(s,t,u) - 1;
    fP2 := transpose matrix { f(R_0, R_1, R_2) };

    camConf := for i to nn-1 list generateCameraMatrix(dimPPImage, dimPPWorld);
    multprojectiveMultiview := for camera in camConf list camera * fP2;
    affineMultiview := matrix for i from 1 to dimPPImage list for view in multprojectiveMultiview list view_(i,0) / view_(0,0);

    pp := symbol pp;
    S := QQ[ gens R | flatten transpose for entryIndex to dimPPImage-1 list for viewIndex to nn-1 list pp_(entryIndex, viewIndex) ];
    -- print flatten transpose for entryIndex to dimPPImage-1 list for viewIndex to nn-1 list pp_(entryIndex, viewIndex);

    imageVars := genericMatrix(S, S_3, dimPPImage, nn);
    -- print imageVars;

    ratGraph := ideal(imageVars - sub((affineMultiview), frac S));
    graph := ideal for i in ratGraph_* list numerator i;
    scan(ratGraph_*, f -> graph = saturate(graph, denominator f));
    modelImage := eliminate({S_0, S_1, S_2}, saturate(graph, ideal(S_0, S_1, S_2)));
    
    -*
    ratMap = flatten entries affineMultiview;
    T = QQ[ flatten transpose for entryIndex to dimPPImage-1 list for viewIndex to nn-1 list pp_(entryIndex, viewIndex) ];
    phi = map(frac R, T, ratMap);
    modelImage = kernel phi;
    *-

    print ("    The degree of the affine multiview variety is " | toString(degree modelImage));
    modelEDDegree := determinantalUnitEDDegree flatten entries gens sub(modelImage, QQ[support modelImage]);
    print ("    The ED degree is " | toString(modelEDDegree));
)

end

----- rational curve case -----
EE = 1;
dimPPWorld = 3;
curve = (s, t) -> for i to dimPPWorld+1 list sum for i to EE list random(1,100) * s^i * t^(EE-i);
for dimPPImage from 2 to 3 do for nn from 1 to 2 do (
    print("For n = "|nn|" and h = "|dimPPImage|", and f a curve");
    theoremEDRationalCurve(nn, dimPPImage, curve)
)

----- rational surface case -----
EE = 1;
dimPPWorld = 3;
surface = (s, t, u) -> for i to dimPPWorld+1 list (
    sum for j from 0 to EE list sum for k from 0 to EE-j list random(1,100) * s^j * t^k * u^(EE-j-k)
)

for dimPPImage from 2 to 4 do for nn from 2 to 2 do (
    print("For n = "|nn|" and h = "|dimPPImage|", and f a surface");
    theoremEDRationalSurface(nn, dimPPImage, surface)
)