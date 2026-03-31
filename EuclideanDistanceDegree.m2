-- TODO: if parameteryHomotopy works, check that we get the same answer as EDHomotopy, but it's probably faster

newPackage(
  "EuclideanDistanceDegree",
  Version => "1.1", 
  Date => "Feb 2026",
  Authors => {
    {Name => "Jose Israel Rodriguez",
    Email => "Jose@Math.wisc.edu",
    HomePage => "http://www.math.wisc.edu/~jose/"},
    {Name => "Will Huang",
    Email => "williamhuang5120@gmail.com",
    HomePage => ""}
  },
  Headline => "Produces equations and computes ED degrees. ",
  --DebuggingMode => false, -- turn off for release builds
  DebuggingMode => true,
  AuxiliaryFiles => false,
  PackageImports => {"SimpleDoc","Bertini","NumericalAlgebraicGeometry","Elimination"},
  PackageExports => {"Bertini","NumericalAlgebraicGeometry"},
  Configuration => {
    --"RandomCoefficients"=>CC,
    "Continuation"=>Bertini },
  CacheExampleOutput => false
)

--path=prepend("/Users/joserodriguez/Documents/GitHub/EuclideanDistanceDegree",path)
--loadPackage("EuclideanDistanceDegree",Reload=>true)
--restart
 
randomCC=()->random CC
randCC=()->random CC
randomRR=()->((-1)^(random(1,2)) *random RR)
randomZZ=()->random(1,30103)
randomValue=(kk)-> if kk===CC then randomCC() else if kk===RR then randomRR() else randomZZ() 
randomVector=method(Options=>{		})
randomVector(ZZ,Thing):= o->(n,R) ->apply(n,i->randomValue(R))--list of length n of randomValue

load"./EDD_Determinantal.m2"
load"./EDD_LeftKernel.m2"
load"./EDD_Numerical.m2"

export {
  "TempDirectory",
  "FileDirectory",
  "ReturnCriticalIdeal",
  "symbolicWeightEDDegree",
  "determinantalUnitEDDegree",
  "determinantalGenericEDDegree",
  "leftKernelWeightEDDegree",
  "leftKernelUnitEDDegree",
  "leftKernelGenericEDDegree",
  "newNumericalComputationOptions",
  "NumericalComputationOptions",
  "homotopyEDDegree",
  "numericWeightEDDegree",
  "numericGenericEDDegree",
  "numericUnitEDDegree",
  "vanishTally"
}

----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------

--##########################################################################--
-- INTERNAL METHODS
--##########################################################################--
parString = (aString) -> ("("|toString(aString)|")");
addSlash = (aString) -> (
  if #aString =!= 0 then (
    if aString_-1 === " " then error(aString | " cannot end with whitespace.");
    if aString_-1 =!= "/" then aString = aString | "/";
  );
  aString
);
makeJac = (system,unknowns) -> (
  --it is a list of lists of partial derivatives of a polynomial
  for i in system list for j in unknowns list diff(j,i)
)
checkZero = (aSol, eps) -> if aSol/abs//min < eps then false else true
sortPointFunction = (aSol) -> (if not (apply(aSol,i->{realPart i,imaginaryPart i}/abs//max)//min<1e-8) then true else false);

--##########################################################################--
-- DOCUMENTATION
--##########################################################################--

beginDocumentation()

doc /// --EuclideanDistanceDegree Package 
  Key
    EuclideanDistanceDegree 
  Headline
    a package to compute Euclidean distance degrees
  Description
    Text
      This package provides several routines for determining the (generic or unit) Euclidean distance degree of an algebraic variety.
      These routines include symbolic methods and numerical methods for determining these degrees. 
    Text
      Using symbolic computation, this code computes the (unit) ED degree of a circle. 		 
    Example
      R = QQ[x,y]
      F = {x^2 + y^2 - 1}
      2 == determinantalUnitEDDegree(F)
    Text
      Using numeric computation, this code computes the (unit) ED degree of a circle. 		 
    Example
      R = QQ[x,y];
      F = {x^2 + y^2 - 1};
      2 == leftKernelUnitEDDegree(F)
    Text
      This package also computes generic ED degrees. The generic ED degree of a variety $X$ is always greater than or equal to the unit 
      ED degree of $X$.
    Example
      R = QQ[x,y];
      F = {x^2 + y^2 - 1};
      4 == determinantalGenericEDDegree(F)
      2 == determinantalUnitEDDegree(F)
    Text
      The most general method for computing ED degrees with symbolic computation is the `symbolicWeightEDDegree` method.
    Example
      R = QQ[x, y];
      F = {x^2 + y^2 - 1};
      genericWeightVector = {2, 3};
      unitWeightVector = {1, 1};
      dataVector = {5, 7};
      4 == symbolicWeightEDDegree(F, dataVector, genericWeightVector)
      2 == symbolicWeightEDDegree(F, dataVector, unitWeightVector)		
    Text
      When the variety is an affine cone, one is able to compute ED degrees using ED degree homotopies, i.e., a structured parameter 
      homotopy. The easiest case is when the variety is a hypersurface (or more generally, a complete intersection)  
    Example
      R = QQ[x1,x2,x3,x4]
      F = G = {det genericMatrix(R,2,2)}
      6 == numericGenericEDDegree(F, G)
      2 == numericUnitEDDegree(F, G)
    Text
      When a $V(F)$ is not a complete intersection we incorporate a membership test to filter out residual critical points. Here $V(F)$
      is an irreducible component of $V(G)$ (a reducible variety) and `#G===codim ideal F`.  These methods employ an equation by equation
      method called regeneration. 
    Example
      R = QQ[x1,x2,x3,x4,x5,x6]
      F = (minors(2, genericMatrix(R,3,2)))_*
      G = drop(F, -1)
      #G == codim ideal F
      10 == numericGenericEDDegree(F, G)
      2 == numericUnitEDDegree(F, G)
    Text
      One may also determine (Unit) ED degrees using a parameter homotopy called a Weight-ED Degree Homotopy. 
    Example
      R = QQ[x1,x2,x3,x4,x5,x6];
      F = (minors(2, genericMatrix(R,3,2)))_*;
      G = drop(F, -1);
      NCO = newNumericalComputationOptions(F, G)
      NCO#"TargetWeight" = apply(#gens R, i->1)
      2 == homotopyEDDegree(NCO, "Weight", true, true)

      NCO#"TargetWeight" = apply(#gens R, i->random RR)
      10 == homotopyEDDegree(NCO, "Weight", false, true)
    Text
      One may also compute critical points for different data using a parameter homotopy called a Data-ED Degree Homotopy. 
    Example
      R = QQ[x1,x2,x3,x4,x5,x6];
      F = (minors(2, genericMatrix(R,3,2)))_*;
      G = drop(F, -1);
      NCO = newNumericalComputationOptions(F, G);
      NCO#"TargetData" = apply(#gens R, i->1)
      10 == homotopyEDDegree(NCO, "Data", true, true)
      importSolutionsFile(NCO#"Directory", NameSolutionsFile => "member_points")

      NCO#"TargetWeight" = {0} | (apply(-1 + #gens R, i->1))
      homotopyEDDegree(NCO, "Data", false, true)
      importSolutionsFile(NCO#"Directory", NameSolutionsFile => "member_points")	
///

doc /// --symbolicWeightEDDegree
  Key
    symbolicWeightEDDegree
    (symbolicWeightEDDegree, List, List, List)
    determinantalUnitEDDegree
    (determinantalUnitEDDegree, List)
    determinantalGenericEDDegree
    (determinantalGenericEDDegree, List)
    ReturnCriticalIdeal
  Headline
    compute Euclidean distance degrees of affine varieties using minors of the augmented Jacobian
  Usage
    UED = determinantalUnitEDDegree(F)
    GED = determinantalGenericEDDegree(F)
    GED = symbolicWeightEDDegree(F, U, W)
  Inputs
    F:List
      a system of polynomials (system need not be square)
    U:List
      a (generic) data vector 
    W:List
      a (generic) weight vector
  Outputs
    GED:ZZ
      a generic Euclidean distance degree 
    UED:ZZ
      a unit Euclidean distance degree 
    ICP:Ideal
      an ideal for the variety of critical points
  Description
    Text
      This method computes Euclidean distance (ED) degrees for the variety defined by the system $F$ via symbolic computation. The critical
      ideal is formed by saturating the ideal defined by $F$ with the ideals of minors of the Jacobian and Augmented Jacobian. The degree of
      the critical ideal is the ED degree of the variety.
    Text
      The default is to return a degree (GED or UED) but the option ReturnCriticalIdeal=>true will change the method to return the critical
      ideal instead. The `determinantalUnitEDDegree` method computes an ED degree using random (integer) data and unit weights, whereas
      `determinantalGenericEDDegree` will use random data and random weights.
    Example
      R = QQ[x,y]
      F = {x^2 + y^2 - 1}
      (U,W) = ({12, 23}, {15, 331})
      UED = determinantalUnitEDDegree F
      GED = determinantalGenericEDDegree F
      ICP = symbolicWeightEDDegree(F, U, W, ReturnCriticalIdeal => true)
///

doc /// --leftKernel
  Key
    leftKernelWeightEDDegree
    (leftKernelWeightEDDegree, List, List, List)
    leftKernelUnitEDDegree
    (leftKernelUnitEDDegree, List)
    leftKernelGenericEDDegree
    (leftKernelGenericEDDegree, List)
    TempDirectory
  Headline
    compute Euclidean distance degrees of affine varieties that are complete intersections using numerical computation
  Usage
    GED = leftKernelWeightEDDegree(F,U,W)
    GED = leftKernelGenericEDDegree(F)
    UED = leftKernelUnitEDDegree(F)
  Inputs
    F:List
      a system of polynomials (need not be square) defining an affine variety that is a complete intersection
    U:List
      a (generic) data vector
    W:List
      a (generic) weight vector
  Outputs
    GED:ZZ
      a generic Euclidean distance degree
    UED:ZZ
      a unit Euclidean distance degree 
  Description
    Text
      This method computes Euclidean distance (ED) degrees for the variety defined by the system $F$ numerically using Lagrange 
      multipliers. The left kernel of the augmented Jacobian matrix is used to define a critical ideal which is passed to Bertini. The 
      Bertini input files are written in a temporary directory and then Bertini is ran to compute critical points.
    Text
      By default, this method creates a temporary directory to store the input files. A specific directory can be specified as a string 
      using the `TempDirectory` option, a directory will be created if it does not exist. The `leftKernelUnitEDDegree` method computes an 
      ED degree using random (complex) data and unit weights, whereas `leftKernelGenericEDDegree` will use random data and random weights.
    Example
      R = QQ[x,y]
      F = {x^2 + y^2 - 1}
      (U,W) = ({.12, .23}, {.15, .331})
      dir = temporaryFileName()
      GED = leftKernelWeightEDDegree(F, U, W, TempDirectory => dir)
      GED = leftKernelGenericEDDegree F
      UED = leftKernelUnitEDDegree F
///

-- Unused keys, maybe for a later version?
--   `OutputType` one of "Standard" or "TestHomotopyConjectureGEDvUED"
--   `BertiniMembershipTestConfiguration`
--   `BertiniSubstitute`
--   `BertiniConstants`
--   `BertiniStartFiberSolveConfiguration`
--   `TrackSolutions`
doc ///
  Key
    NumericalComputationOptions
    newNumericalComputationOptions
    (newNumericalComputationOptions, String, List, List, List)
    (newNumericalComputationOptions, String, List, List)
    (newNumericalComputationOptions, List, List, List)
    (newNumericalComputationOptions, List, List)
  Headline
    define homotopy options and configurations
  Usage
    NCO = newNumericalComputationOptions(dir, F, G, L)
    NCO = newNumericalComputationOptions(dir, F, G)
    NCO = newNumericalComputationOptions(F, G, L)
    NCO = newNumericalComputationOptions(F, G)
  Inputs
    dir:String
      directory to write Bertini files to, will be created if it does not exist
    F:List
      a system of polynomials (need not be square) defining the variety
    G:List
      a system of polynomials (complete intersection) such that V(F) is an irreducible component of V(G).
    L:List
      a system of linear polynomials
  Outputs
    NCO:NumericalComputationOptions
      a MutableHashTable to keeps track of the options and configurations for the homotopy methods
  Description
    Text
      A `NumericalComputationOptions` object stores all the options needed to define a homotopy run with Bertini. At mimumum it requires
      the model $F$ for which the ED Degree will be computed and a witness model $G$. A submodel $L$ may be passed in, by default it will
      be the empty set. A String indicating the temporary directory from which Bertini will read/write files during the run may also be
      passed in as an optional argument, by default it will be a random temporary file name. The temporary directory will be created
      if it does not exist.
    Text
      Keys available to customize the homotopy include: 
      `Directory` to change the directory Bertini is run from, the directory need not exist.
      `StartData` and `TargetData` if executing a Data homotopy, defaults to random complex data. 
      `StartWeight` and `TargetWeight` if executing a Weight homotopy, defaults to random complex and unit weights respectively. 
      `Infinity` to add a hyperplane at infinity, one is create by default in @TO homotopyEDDegree@ if not present. 
      `PrimalCoordinates` to add additional variables to the ambient space. By default this field will already contain the variables of 
      `ring F`, thus it is safer to append extra variables rather than overwriting the less, otherwise subsequence methods may fail. 
    Text
      Keys are also available to customize the Bertini run: 
      `HomogeneousVariableGroups` and `AffineVariableGroups` to modify variable groups. By default the variables of `ring F` are treated
      as Homogeneous. 
      `FinerRestriction` a list of polynomials to filter down critical points. Critical points will only be kept if they vanish
      for every polynomial in this list.
    Example
      R = QQ[x,y]
      F = G = {x^2 + y^2 - 1}
      dir = temporaryFileName()
      NCO = newNumericalComputationOptions(dir, F, G)
///

doc /// --numeric
  Key
    homotopyEDDegree   
    (homotopyEDDegree, NumericalComputationOptions, String, Boolean, Boolean)
    numericWeightEDDegree
    (numericWeightEDDegree, List, List, List, List)
    numericUnitEDDegree
    (numericUnitEDDegree, List, List)
    numericGenericEDDegree
    (numericGenericEDDegree, List, List)
  Headline
    compute Euclidean distance degrees of projective varieties (affine cones) numerically using homotopy continuation
  Usage
    GED = numericGenericEDDegree(F,G)
    UED = numericUnitEDDegree(F,G)
    GED = homotopyEDDegree(NCO, "Weight", true, false)
  Inputs
    NCO:NumericalComputationOptions
      a MutableHashTable to keeps track of the options and configurations for the homotopy methods
    F:List
      a system of polynomials (need not be square) defining the variety
    G:List
      a system of polynomials (complete intersection) such that V(F) is an irreducible component of V(G).
    L:List
      a system of linear polynomials
    U:List
      a (generic) data vector
    W:List
      a (generic) weight vector
    ht:String
      one of "Weight", "Data", or "Submodel" which defines the type of homotopy to execute
    isStageOne:Boolean
      indicates the stage of the homotopy continuation
    isStageTwo:Boolean
      indicates the stage of the homotop continuation
  Outputs
    GED:ZZ
      a generic Euclidean distance degree
    UED:ZZ
      a unit Euclidean distance degree
  Description
    Text
      Uses a weight homotopy by default, though data and submodel homotopies are available. The Bertini input files are written in dir
      and then Bertini is ran. The `numericWeightEDDegree`, `numericUnitEDDegree`, and `numericGenericEDDegree` methods are all weight
      homotopy methods. Finer control over the homotopy is accomplished using the `homotopyEDDegree` method by modifying the
      @TO NumericalComputationOptions@ object.
    Text
      By default, this method creates a temporary directory to store the input files. A specific directory can be specified as a string 
      using the `TempDirectory` option, a directory will be created if it does not exist. The `numericUnitEDDegree` method computes an 
      ED degree using random (complex) data and unit weights, whereas `numericGenericEDDegree` will use random data and random weights.
      The `homotopyEDDegree` method uses the directory defined in the @TO NumericalComputationOptions@ object.
    Example
      R = QQ[x,y]
      F = G = {x^2+y^2-1}
      GED = numericGenericEDDegree(F, G)
      UED = numericUnitEDDegree(F, G)

      NCO = newNumericalComputationOptions(F, G)
      NCO#"TargetWeight"
      GED = homotopyEDDegree(NCO, "Weight", true, false)
      UED = homotopyEDDegree(NCO, "Weight", true, true)

      dir = temporaryFileName()
      U = {1, 2}
      W1 = {1, 1}
      WS = {.7, 1.2}
      GED = numericWeightEDDegree(F, G, U, WS, TempDirectory => dir)
      UED = numericWeightEDDegree(F, G, U, W1, TempDirectory => dir)
///

doc /// --vanishTally
  Key
    vanishTally
    (vanishTally, NumericalComputationOptions, Ideal, RR)
    (vanishTally, NumericalComputationOptions, List, RR)
    (vanishTally, NumericalComputationOptions, Ideal)
    (vanishTally, NumericalComputationOptions, List)
  Headline
    validate a Bertini homotopy continuation run
  Usage
    counts = vanishTally(NCO, Z, 1e-8)
  Inputs
    NCO:NumericalComputationOptions
      a MutableHashTable storing configuration for a Bertini run
    Z:Ideal
      test ideal
    F:List
      a list of polynomials defining a test ideal
    setTolerance:ZZ
      a numerical tolerance (default 1e-8)
  Outputs
    counts:Tally
      A tally of paths which vanish (up to setTolerance) on the test ideal
  Description
    Text
      Reads a data file from a Bertini run and tallies the points which lie on the given test ideal up to the given tolerance. Points are
      tallied according to the path number assigned by Bertini.
    Example
      R = QQ[x,y]
      F = G = {x^2 + y^2 - 1}
      NCO = newNumericalComputationOptions(F, G)
      GED = homotopyEDDegree(NCO, "Weight", true, true)
      vanishTally(NCO, F)
  Caveat
    This method assumes that a Bertini run was completed in the directory specified by NCO.
///

--##########################################################################--
-- END DOCUMENTATION
--##########################################################################--

TEST /// -- basic examples
  R = QQ[x0, x1, x2];
  F = {x0^2*x2 - x1^2*(x1 + x2)};
  assert(determinantalGenericEDDegree(F) === 7);
  assert(determinantalUnitEDDegree(F) === 7);

  R = QQ[jj]/ideal(jj^2+1)[x0,x1,x2];
  F = {x0^2*x1 -(x1 - jj*x2)^2*x2};
  -- The output has a factor of 2 to account for the imaginary unit jj.
  assert(determinantalGenericEDDegree(F) === 2*7);
  assert(determinantalUnitEDDegree(F) === 2*7);

  R = QQ[x0,x1,x2,x3];
  F = {x0^2*x1-x2*x3^2};
  assert(determinantalGenericEDDegree(F) === 10);
  assert(determinantalUnitEDDegree(F) === 10);

  R = QQ[x,y];
  F = {x^2 + y^2 - 1};
  (U, W) = ({12, 23}, {15, 331});
  assert(symbolicWeightEDDegree(F, U, W) === 4);
  assert(determinantalGenericEDDegree F === 4);
  assert(determinantalUnitEDDegree F === 2);
  
  R = QQ[x,y];
  F = G = {x^2+y^2-1};
  (U, W) = ({.12, .23}, {.15, .331})
  assert(leftKernelWeightEDDegree(F, U, W), === 4);
  assert(leftKernelGenericEDDegree F === 4);
  assert(leftKernelUnitEDDegree F === 2);

  assert(numericWeightEDDegree(F, G, U, W) === 4);
  assert(numericGenericEDDegree(F, G) === 4);
  assert(numericUnitEDDegree(F, G) === 2);
///

-- TODO: add remaining surfaces (73 total)
-- TODO: random seed issues?
-- https://homepage.univie.ac.at/herwig.hauser/bildergalerie/gallery.html (check all 3 give the same result)
-- https://www.imaginary.org/gallery/herwig-hauser-classic
TEST ///  -- Herwig Hauser's algebraic surfaces gallery
  R = QQ[x,y,z];
  surfaces = {
    x^2 + y^2*z^3 - z^4,
    x^2 + y^2*z - z^2,
    x^3*y + x*z^3 + y^3*z + z^3 + 7*z^2 + 5*z,
    x^6 + y^6 + z^6 - 1,
    3*x^2 + 3*y^2 + z^2 - 1,
    (x^2 - y^3)^2 - (z^2 - y^2)^3,
    x^2 + y^2 + z^3 - z^2,
    x^2 + y^2 + z^2 + 1000 * (x^2 + y^2) * (x^2 + z^2) * (y^2 + z^2) - 1
  };

  for surface in surfaces do (
    F = {surface};

    UED_symb = determinantalUnitEDDegree F;
    GED_symb = determinantalGenericEDDegree F;
    UED_left = leftKernelUnitEDDegree F;
    GED_left = leftKernelGenericEDDegree F;
    UED_numeric = numericUnitEDDegree(F, F);
    GED_numeric = numericGenericEDDegree(F, F);

    assert(UED_symb === UED_left);
    assert(UED_left === UED_numeric);
    assert(GED_symb === GED_left);
    assert(GED_left === UED_numeric);
  )
///

TEST ///  -- PNN function space
  d = (3,1,1);
  r = 2;

  -- Parameter space: weights
  R = QQ[W_0..W_(d_1 * d_0), V_0..V_(d_2 * d_1)];
  W = genericMatrix(R, W_0, d_1, d_0);
  V = genericMatrix(R, V_0, d_2, d_1);

  -- Function space
  T = R[x_0..x_(d_0 - 1)];
  X = transpose matrix{apply(d_0, i -> x_i)};  -- input
  Z = W * X;
  A = matrix table(d_1, 1, (i, j) -> (Z_(i,j))^r);
  Phi = V * A;  -- output function

  -- Create the ambient space (space of symmetric tensors)
  mons = flatten entries basis(r, T);  -- homogeneous monomials of degree r
  S = QQ[c_0..c_(d_2 * #mons - 1)];

  -- Get image of the parameterization map
  im = flatten apply(d_2, i -> (
    f = Phi_(i,0);
    apply(mons, m -> coefficient(m, f))
  ));

  -- Get the defining ideal
  paramMap = map(R, S, im);
  I = kernel paramMap;
  c = codim I;

  -- Prune variety down to a complete intersection
  gensI = flatten entries gens I;
  F = apply(c, i -> sum apply(#gensI, j -> (random(QQ) * gensI_j)));
  assert(leftKernelUnitEDDegree(F) === numericUnitEDDegree(F, F));
  assert(leftKernelGenericEDDegree(F) === numericGenericEDDegree(F, F));
///

TEST ///  -- self-attention network function space
  d = (1,2,1);
  l = #d - 1;
  t = 1;
  
  -- Parameter space: weights
  ARing = QQ[A_0..A_(d_0 * d_0), B_0..B_(d_1 * d_1)];  -- ring for attention matrices
  VRing = QQ[V_0..V_(d_1 * d_0), W_0..W_(d_2 * d_1)];  -- ring for value matrices
  R = ARing ** VRing;

  As = (genericMatrix(R, A_0, d_0, d_0), genericMatrix(R, B_0, d_1, d_1));
  Vs = (genericMatrix(R, V_0, d_1, d_0), genericMatrix(R, W_0, d_2, d_1));

  -- Function space
  T = R[x_0..x_(d_0 * t)];
  X = genericMatrix(T, x_0, d_0, t);  -- input
  X = Vs_0 * X * transpose X * As_0 * X;
  Phi = Vs_1 * X * transpose X * As_1 * X;  -- output

  -- Create the ambient space (space of symmetric tensors)
  mons = flatten entries basis(3^l, T);  -- homogeneous monomials of degree 3^l
  S = QQ[c_0..c_(d_2 * t * #mons - 1)];

  -- Get image of the parameterization map
  PhiFlat = flatten entries Phi;
  im = flatten apply(d_2 * t, i -> (
    f = PhiFlat_i;
    apply(mons, m -> coefficient(m, f))
  ));

  -- Get the defining ideal
  paramMap = map(R, S, im);
  I = kernel paramMap;
  c = codim I;

  gensI = flatten entries gens I;
  F = apply(c, i -> sum apply(#gensI, j -> (random(QQ) * gensI_j)));
  assert(leftKernelUnitEDDegree(F) === numericUnitEDDegree(F, F));
  assert(leftKernelGenericEDDegree(F) === numericGenericEDDegree(F, F));
///

-*
TEST ///
--load concatenate(MultiprojectiveWitnessSets#"source directory","./AEO/TST/Example1.tst.m2")
///
*-

end

restart
loadPackage "EuclideanDistanceDegree"
installPackage "EuclideanDistanceDegree"
check "EuclideanDistanceDegree"


-- Debugging surface tests
loadPackage "EuclideanDistanceDegree"
R = QQ[x,y,z];
surfaces = {
  x^2 + y^2*z^3 - z^4,
  x^2 + y^2*z - z^2,
  x^3*y + x*z^3 + y^3*z + z^3 + 7*z^2 + 5*z,
  x^6 + y^6 + z^6 - 1,
  3*x^2 + 3*y^2 + z^2 - 1,
  (x^2 - y^3)^2 - (z^2 - y^2)^3,
  x^2 + y^2 + z^3 - z^2,
  x^2 + y^2 + z^2 + 1000 * (x^2 + y^2) * (x^2 + z^2) * (y^2 + z^2) - 1  -- this one (left kernel disagrees on unit)
};

UED1 = {};
UED2 = {};
UED3 = {};

GED1 = {};
GED2 = {};
GED3 = {};

for surface in surfaces do (
  F = {surface};
  UED1 = UED1 | {determinantalUnitEDDegree F};
  UED2 = UED2 | {leftKernelUnitEDDegree F};
  UED3 = UED3 | {numericUnitEDDegree(F, F)};

  GED1 = GED1 | {determinantalGenericEDDegree F};
  GED2 = GED2 | {leftKernelGenericEDDegree F};
  GED3 = GED3 | {numericGenericEDDegree(F, F)};
)

print("Comparing unit degrees")
scan(#UED1, i -> (
  if UED1_i =!= UED2_i then print("symb-left F: " | toString surface | "; " | UED1_i | ", " | UED2_i);
  if UED2_i =!= UED3_i then print("left-num F: " | toString surface | "; " | UED2_i | ", " | UED3_i);
  if UED1_i =!= UED3_i then print("sym-num F: " | toString surface | "; " | UED1_i | ", " | UED3_i);
))

print("Comparing generic degrees")
scan(#GED1, i -> (
  if GED1_i =!= GED2_i then print("symb-left F: " | toString surface | "; " | GED1_i | ", " | GED2_i);
  if GED2_i =!= GED3_i then print("left-num F: " | toString surface | "; " | GED2_i | ", " | GED3_i);
  if GED1_i =!= GED3_i then print("sym-num F: " | toString surface | "; " | GED1_i | ", " | GED3_i);
))