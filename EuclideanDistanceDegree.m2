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

-- TODO: double check examples
doc /// --EuclideanDistanceDegree  
  Key
    EuclideanDistanceDegree 
  Headline
    a package to determine Euclidean distance degrees
  Description
    Text
      This package provides several routines for determining the (generic or unit) Euclidean distance degree of an algebraic variety.
      These routines include symbolic methods and numerical methods for determining these degrees. 
    Text
      Using symbolic computation, this code computes the (unit) ED degree of a circle. 		 
    Example
      R=QQ[x,y];
      F={x^2+y^2-1};
      2==determinantalUnitEDDegree(F)
    Text
      Using numeric computation, this code computes the (unit) ED degree of a circle. 		 
    Example
      R=QQ[x,y];
      F={x^2+y^2-1};
      2==leftKernelUnitEDDegree(F)
    Text
      This package also computes generic ED degrees. The generic ED degree of X is always greater than or equal to the unit ED degree X.
    Example
      R=QQ[x,y];
      F={x^2+y^2-1};
      genericWeightVector={2,3};
      unitWeightVector={1,1};
      dataVector={5,7};
	    4==symbolicWeightEDDegree(F,dataVector,genericWeightVector)
	    2==symbolicWeightEDDegree(F,dataVector,unitWeightVector)
    Text
      The most general method for computing ED degrees with symbolic computation is symbolicWeightEDDegree
    Example
   	  R=QQ[x,y];
      F={x^2+y^2-1};
      genericWeightVector={2,3};
      unitWeightVector={1,1};
      dataVector={5,7};
      4==symbolicWeightEDDegree(F,dataVector,genericWeightVector)
      2==symbolicWeightEDDegree(F,dataVector,unitWeightVector)		
    Text
      When the variety is an affine cone, one is able to compute ED degrees using ED degree homotopies, i.e., a structured parameter homotopy.
	  The easiest case is when the variety is a hypersurface (or more generally, a complete intersection)  
    Example
      R=QQ[x1,x2,x3,x4]
	    F={det genericMatrix(R,2,2)};
      P=(F,F)
      6==numericEDDegree(P,"Generic")
	    2==numericEDDegree(P,"Unit")
    Text
      When a V(F) is not a complete intersection we incorporate a membership test to filter out residual critical points.
      Here V(F) is an irreducible component of V(G) (a reducible variety) and #G===codim ideal F.
      These methods employ an equation by equation method called regeneration. 
    Example
      R=QQ[x1,x2,x3,x4,x5,x6]
	    F=(minors(2,genericMatrix(R,3,2)))_*;
      G=drop(F,-1);	
      P=(F,G)
      #G==codim ideal F;
      10==numericEDDegree(P,"Generic")
      2==numericEDDegree(P,"Unit")
    Text
      One may also determine (Unit) ED degrees using a parameter homotopy called a Weight-ED Degree Homotopy. 
    Example
      dir=temporaryFileName();if not fileExists dir then mkdir dir;
      R=QQ[x1,x2,x3,x4,x5,x6]
	    F=(minors(2,genericMatrix(R,3,2)))_*;
      G=drop(F,-1);	
      #G==codim ideal F;
      P=(F,G)
	    NCO=newNumericalComputationOptions(dir,P)
      NCO#"TargetWeight"=apply(#gens R,i->1)
	    2==homotopyEDDegree(NCO,"Weight",true,true)
      NCO#"TargetWeight"=(apply(#gens R,i->random RR))
	    10==homotopyEDDegree(NCO,"Weight",false,true)
    Text
      One may also compute crtiical points for different data using a parameter homotopy called a Data-ED Degree Homotopy. 
    Example
      dir=temporaryFileName();if not fileExists dir then mkdir dir;
      R=QQ[x1,x2,x3,x4,x5,x6]
	    F=(minors(2,genericMatrix(R,3,2)))_*;
      G=drop(F,-1);	
      #G==codim ideal F;
      P=(F,G)
	    NCO=newNumericalComputationOptions(dir,P)
      NCO#"TargetData"=apply(#gens R,i->1)
      10==homotopyEDDegree(NCO,"Data",true,true)
      importSolutionsFile(NCO#"Directory",NameSolutionsFile=>"member_points")	
      NCO#"TargetWeight"={0}|(apply(-1+#gens R,i->1))
      homotopyEDDegree(NCO,"Data",false,true)
      importSolutionsFile(NCO#"Directory",NameSolutionsFile=>"member_points")	
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
    GED = symbolicWeightEDDegree(F,U,W)
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

      The default is to return a degree (GED or UED) but the option ReturnCriticalIdeal=>true will change the method to return the critical
      ideal instead. The `determinantalUnitEDDegree` method computes an ED degree using random (integer) data and unit weights, whereas
      `determinantalGenericEDDegree` will use random data and random weights.
    Example
      R = QQ[x, y];     
      F = {x^2 + y^2 - 1}
      (U,W) = ({12, 23}, {15, 331})
      UED = determinantalUnitEDDegree F 
      GED = determinantalGenericEDDegree F 
      ICP = symbolicWeightEDDegree(F, U, W, ReturnCriticalIdeal => true)
  Caveat
    none
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
    compute Euclidean distance degrees of affine varieties using the left kernel of the augmented Jacobian
  Usage
    GED = leftKernelWeightEDDegree(F,U,W)
    GED = leftKernelGenericEDDegree(F)
    UED = leftKernelUnitEDDegree(F)
  Inputs
    F:List
      a system of polynomials (need not be square)
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
      This method computes Euclidean distance (ED) degrees for the variety defined by the system $F$ via Lagrange multipliers. The augmented
      matrix is constructed and its left kernel is used to construct a critical ideal which is passed to Bertini. The Bertini input files 
      are written in a temporary directory and then Bertini is ran to compute critical points.

      By default, this method creates a temporary directory to store the input files. A specific directory can be specified as a string 
      using the `TempDirectory` option, a directory will be created if it does not exist. The `leftKernelUnitEDDegree` method computes an 
      ED degree using random (complex) data and unit weights, whereas `leftKernelGenericEDDegree` will use random data and random weights.
    Example
      R = QQ[x,y];
      F = {x^2 + y^2 - 1}
      (U,W) = ({1, 2}, {.15, .331})
      dir = temporaryFileName();
      GED = leftKernelWeightEDDegree(F, U, W, FileDirectory => dir);
      GED = leftKernelGenericEDDegree F
      UED = leftKernelUnitEDDegree F
  Caveat
    none
///

-- TODO: Split off NumericalComputationOptions (?)
-*
doc ///
  Key
    NumericalComputationOptions
    newNumericalComputationOptions
    (newNumericalComputationOptions, Sequence)
    TempDirectory
  Headline
  Usage
  Inputs
  Outputs
  Description
    Text
    Example
  Caveat
    None
///
*-

-- TODO: Document the new methods numericUnitEDDegree, numericGenericEDDegree
doc /// --numeric
  Key
    homotopyEDDegree   
    (homotopyEDDegree,NumericalComputationOptions,String,Boolean,Boolean)
    numericWeightEDDegree
    (numericWeightEDDegree, Sequence, List, List)
    numericUnitEDDegree
    (numericUnitEDDegree, Sequence)
    numericGenericEDDegree
    (numericGenericEDDegree, Sequence)
  Headline
    compute Euclidean distance degrees of projective varieties (affine cones) using numerical computation
  Usage
    GED = numericEDDegree((F,G),"Generic")
    UED = numericEDDegree((F,G),"Unit")
    NCO = newNumericalComputationOptions(dir,(F,G))
    (NCO = newNumericalComputationOptions(dir,(F,G)); GED = homotopyEDDegree(NCO, "Weight", true, false))
    (NCO = newNumericalComputationOptions(dir,(F,G)); UED = homotopyEDDegree(NCO, "Weight", true, true))
  Inputs
    F:List
      polynomials 
    G:List
      polynomials (complete intersection) such that V(F) is an irreducible component of V(G).
    dir:String
      a directory 
    NCO:NumericalComputationOptions
      a MutableHashTable to keeps track of the options and configurations for the homotopy methods
  Outputs
    GED:ZZ
      a generic Euclidean distance degree 
    UED:ZZ
      a unit Euclidean distance degree 
    NCO:NumericalComputationOptions
      a MutableHashTable to keeps track of the options and configurations for the homotopy methods
    WS:List
      a (generic) start weight vector
  Description
    Text
      The Bertini input files are written in dir and then Bertini is ran.  
    Example
      R = QQ[x,y];     
      F = G = {x^2+y^2-1}
      W1 = {1,1}
      WS = {.7,1.2}
      dir=temporaryFileName(); mkdir dir;
      GED = numericEDDegree((F,G),"Generic")
      UED = numericEDDegree((F,G),"Unit")
      NCO = newNumericalComputationOptions(dir,(F,G))
      NCO#"TargetWeight"
      GED = homotopyEDDegree(NCO, "Weight", true, false)
      UED = homotopyEDDegree(NCO, "Weight", true, true)
      GED = numericWeightEDDegree((F,G),WS)
      UED = numericWeightEDDegree((F,G),W1)
  Caveat
    none
///

-- TODO: finish documenting this
doc /// --vanishTally
  Key
    vanishTally
    (vanishTally, NumericalComputationOptions, Ideal, RR)
    (vanishTally, NumericalComputationOptions, List, RR)
    (vanishTally, NumericalComputationOptions, Ideal)
    (vanishTally, NumericalComputationOptions, List)
  Headline
    compute Euclidean distance degrees using symbolic computation
  Usage
    counts = vanishTally(NCO, Z, 1e-8)
  Inputs
    NCO:NumericalComputationOptions
      ll
    Z:Ideal
      an ideal for the variety of critical points
    setTolerance:ZZ
      ll
    fegZ:List
      ll
  Outputs
    counts:Tally
      ll
  Description
    Text
      The default is to return a degree (GED or UED) but the option ReturnCriticalIdeal=>true will change the option to ICP 
    Example
      R = QQ[x,y,z]
  Caveat
    none
///

--##########################################################################--
-- END DOCUMENTATION
--##########################################################################--

TEST ///
  R = QQ[x0,x1,x2];
  F = {x0^2*x2 - x1^2*(x1 + x2)};
  assert(determinantalGenericEDDegree(F) === 7);
  assert(determinantalUnitEDDegree(F) === 7);

  R = QQ[jj]/ideal(jj^2+1)[x0,x1,x2];
  F = {x0^2*x1 -(x1 - jj*x2)^2*x2};
  -- The output has a factor of 2 to account for the imaginary unit jj.
  assert(determinantalGenericEDDegree(F) === 2*7)
  assert(determinantalUnitEDDegree(F) === 2*7);

  R = QQ[x0,x1,x2,x3];
  F = {x0^2*x1-x2*x3^2};
  assert(determinantalGenericEDDegree(F) === 10)
  assert(determinantalUnitEDDegree(F) === 10)
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
    UED_numeric = numericUnitEDDegree F;
    GED_numeric = numericGenericEDDegree F;

    assert(UED_symb === UED_left);
    assert(UED_left === UED_numeric);
    assert(GED_symb === GED_left);
    assert(GED_left === UED_numeric);
  )
///

TEST ///  -- PNN function space
  d = (3,2,3);  -- layers
  r = 2; -- activation function

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
  image = flatten apply(d_2, i -> (
    f = Phi_(i,0);
    apply(mons, m -> coefficient(m, f))
  ));

  -- Get the defining ideal
  paramMap = map(R, S, image);
  I = kernel paramMap;
  F = gens I;

  assert(leftKernelUnitEDDegree(F) === numericUnitEDDegree(F));
  assert(leftKernelGenericEDDegree(F) === numericGenericEDDegree(F))
///

TEST ///  -- Attention mechanism
///

TEST ///
--load concatenate(MultiprojectiveWitnessSets#"source directory","./AEO/TST/Example1.tst.m2")
///

end

restart
loadPackage "EuclideanDistanceDegree"
installPackage "EuclideanDistanceDegree"
check "EuclideanDistanceDegree"
