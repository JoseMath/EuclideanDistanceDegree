
newPackage(
    "EuclideanDistanceDegree",
    Version => "1.1", 
    Date => "Feb 2026",
    Authors => {
   {Name => "Jose Israel Rodriguez",
       Email => "Jose@Math.wisc.edu",
       HomePage => "http://www.math.wisc.edu/~jose/"}
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
--    "vanishTally", --Needs to be document if exported
    "ReturnCriticalIdeal",
    "homotopyEDDegree",
    "symbolicWeightEDDegree",
    "determinantalUnitEuclideanDistanceDegree",
    "determinantalGenericEuclideanDistanceDegree",
    "leftKernelWeightEDDegree",
    "leftKernelUnitEDDegree",
    "leftKernelGenericEDDegree",
    "runBertiniEDDegree",
    "newNumericalComputationOptions",
    "NumericalComputationOptions",
    "numericWeightEDDegree",
    "numericEDDegree"
            }

----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------


--##########################################################################--
-- INTERNAL METHODS
--##########################################################################--
parString = (aString)->("("|toString(aString)|")");
addSlash = (aString) -> (
    if #aString =!= 0 then (
	if aString_-1 === " " then error (aString | " cannot end with whitespace.");
	if aString_-1 =!= "/" then aString = aString | "/";
	);
    aString
    );
makeJac = (system,unknowns)->(--it is a list of lists of partial derivatives of a polynomial
         for i in system list for j in unknowns list  diff(j,i))
checkZero=(aSol,eps)->if aSol/abs//min<eps then false else true
sortPointFunction = (aSol)->(if not (apply(aSol,i->{realPart i,imaginaryPart i}/abs//max)//min<1e-8) then true else false	    );

--beginDocumentation()
--load "./DOC_EDD.m2";
--##########################################################################--
-- DOCUMENTATION
--##########################################################################--

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
	2==determinantalUnitEuclideanDistanceDegree(F)
      Text
      	Using numeric computation, this code computes the (unit) ED degree of a circle. 		 
      Example
    	R=QQ[x,y];
	F={x^2+y^2-1};	c=1;
    	leftKernelUnitEDDegree(storeBM2Files,F)
	2==runBertiniEDDegree(storeBM2Files)
      Text
      	This package also computes generic ED degrees. The generic ED degree of X is always greater than or equal to the unit ED degree X.
      Example
    	R=QQ[x,y];
	F={x^2+y^2-1};	c=1
    	4==determinantalGenericEuclideanDistanceDegree(F)
    	leftKernelGenericEDDegree(storeBM2Files,F)
	4==runBertiniEDDegree(storeBM2Files)
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

///;


doc ///--symbolicWeightEDDegree
 Key
   symbolicWeightEDDegree
   (symbolicWeightEDDegree,List,List,List)
   determinantalUnitEuclideanDistanceDegree
   (determinantalUnitEuclideanDistanceDegree,List)
   determinantalGenericEuclideanDistanceDegree
   (determinantalGenericEuclideanDistanceDegree,List)
   ReturnCriticalIdeal
 Headline
   a method to determine Euclidean distance degrees using symbolic computation
 Usage
   UED = determinantalUnitEuclideanDistanceDegree(F)
   GED = determinantalGenericEuclideanDistanceDegree(F)
   GED = symbolicWeightEDDegree(F,U,W)
 Inputs
   F:List
     polynomials (system need not be square)
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
     The default is to return a degree (GED or UED) but the option ReturnCriticalIdeal=>true will change the option to ICP 
   Example
     R = QQ[x,y];     
     F = {x^2+y^2-1}
     (U,W) = ({12,23},{15,331})
     UED = determinantalUnitEuclideanDistanceDegree F 
     GED = determinantalGenericEuclideanDistanceDegree F 
     ICP = symbolicWeightEDDegree(F,U,W,ReturnCriticalIdeal=>true)
 Caveat
   none
      
///;



doc ///--NumericalComputationOptions
 Key
   NumericalComputationOptions
   newNumericalComputationOptions
   (newNumericalComputationOptions,String,Sequence)
   homotopyEDDegree   
   (homotopyEDDegree,NumericalComputationOptions,String,Boolean,Boolean)
   numericWeightEDDegree
   (numericWeightEDDegree,String,Sequence,List)
   (numericWeightEDDegree,Sequence,List)
   numericEDDegree
   (numericEDDegree,Sequence,String)
 Headline
   a method to determine Euclidean distance degrees of projective varieties (affine cones) using numerical computation
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
      
///;


doc ///--leftKernel
 Key
   leftKernelWeightEDDegree
   (leftKernelWeightEDDegree,String,List,List)
   leftKernelUnitEDDegree
   (leftKernelUnitEDDegree,String,List)
   leftKernelGenericEDDegree
   (leftKernelGenericEDDegree,String,List)
   runBertiniEDDegree
   (runBertiniEDDegree,String)
 Headline
   a method to determine Euclidean distance degrees of affine varieties that are complete intersections using numerical computation
 Usage
   GED = runBertiniEDDegree leftKernelWeightEDDegree(dir,c,F,W)
   GED = runBertiniEDDegree leftKernelGenericEDDegree(dir,c,F)
   UED = runBertiniEDDegree leftKernelUnitEDDegree(dir,c,F)
 Inputs
   F:List
     polynomials (system need not be square)
   dir:String
     a directory 
   c:ZZ
     a codimension
   W:List
     a (generic) weight vector
 Outputs
   GED:ZZ
     a generic Euclidean distance degree 
   UED:ZZ
     a unit Euclidean distance degree
 Description
   Text
     The Bertini input files are written in dir and then Bertini is ran.  
   Example
     R = QQ[x,y];     
     F = {x^2+y^2-1}
     W = {.15,.331}
     dir=temporaryFileName(); mkdir dir;
     GED = runBertiniEDDegree leftKernelWeightEDDegree(dir,F,W)
     GED = runBertiniEDDegree leftKernelGenericEDDegree(dir,F)
     UED = runBertiniEDDegree leftKernelUnitEDDegree(dir,F)
 Caveat
   none
      
///;



--##########################################################################--
-- END DOCUMENTATION
--##########################################################################--


TEST///
--load concatenate(MultiprojectiveWitnessSets#"source directory","./AEO/TST/Example1.tst.m2")
///


end
