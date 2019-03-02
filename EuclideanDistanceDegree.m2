
newPackage(
    "EuclideanDistanceDegree",
    Version => "1.0", 
    Date => "January 2019",
    Authors => {
   {Name => "Jose Israel Rodriguez",
       Email => "Jose@Math.wisc.edu",
       HomePage => "http://www.math.wisc.edu/~jose/"}
    },
    Headline => "Produces equations and computes ED degrees. ",
    DebuggingMode => true,
    AuxiliaryFiles => false,
    PackageImports => {"SimpleDoc","Bertini"},
    PackageExports => {"SimpleDoc","Bertini"},
  Configuration => { "RandomCoefficients"=>CC,
      "Continuation"=>Bertini },
  CacheExampleOutput => false
)


--path=prepend("/Users/jo/Documents/GoodGit/EuclideanDistanceDegree",path)
--loadPackage("EuclideanDistanceDegree",Reload=>true)


randomCC=()->random CC
randCC=()->random CC
randomRR=()->((-1)^(random(1,2)) *random RR)
randomZZ=()->random(1,30103)
randomValue=(kk)-> if kk===CC then randomCC() else if kk===RR then randomRR() else randomZZ() 
randomVector=method(Options=>{		})
randomVector(ZZ,Thing):= o->(n,R) ->apply(n,i->randomValue(R))--list of length n of randomValue

load"EDD_Determinantal.m2"
load"EDD_LeftKernel.m2"
load"EDD_Numerical.m2"


export { 
--load"EDD_Determinantal.m2"
    "symbolicWeightEDDegree",
    "determinantalUnitEuclideanDistanceDegree",
    "determinantalGenericEuclideanDistanceDegree",
    --load"EDD_LeftKernel.m2"
    "leftKernelWeightEDDegree",
    "leftKernelUnitEDDegree",
    "leftKernelGenericEDDegree",
    "runBertiniEDDegree",
    "writeLeftKernelProjectiveGenericEDDegree",
    "writeLeftKernelProjectiveUnitEDDegree",
    "runBertiniProjectiveEDDegree",
    "startEDDegree","runBertiniStartEDDegree",
--Options
    "Data","Weight","UseRegeneration","TargetSlice"
        }


--###################################
-- TYPE DEFINITIONS
--###################################
--RemovalMLDegree=new Type of MutableHashTable





 



----------------------------------------------------------------------------------------------------------------
checkZero=(aSol,eps)->if aSol/abs//min<eps then false else true
----------------------------------------------------------------------------------------------------------------
sortPointFunction=(aSol)->(if not (apply(aSol,i->{realPart i,imaginaryPart i}/abs//max)//min<1e-8) then true else false	    );
----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------


--##########################################################################--
-- INTERNAL METHODS
--##########################################################################--
----------------------------------------
parString=(aString)->("("|toString(aString)|")");
addSlash=(aString)->(
    if aString_-1===" " then error (aString|" cannot end with whitespace.");
    if aString_-1=!="/" then aString=aString|"/";
    return aString    )
--newHyperplanes=A->for i to (numColumns A)+1 list randomVector(numRows A)
makeJac=(system,unknowns)->(--it is a list of lists of partial derivatives of a polynomial
         for i in system list for j in unknowns list  diff(j,i))


beginDocumentation()

--load "./EuclideanDistanceDegree/DOC.m2";

TEST///
--load concatenate(MultiprojectiveWitnessSets#"source directory","./AEO/TST/Example1.tst.m2")
///


end

(mRow,nCol)=(3,3)
R=QQ[x_(1,1)..x_(mRow,nCol)]
M=transpose genericMatrix(R,mRow,nCol)

 determinantalUnitEuclideanDistanceDegree flatten entries gens minors(2,M)

win1= determinantalGenericEuclideanDistanceDegree flatten entries gens minors(3,M)

win2= determinantalGenericEuclideanDistanceDegree flatten entries gens minors(2,M)


R=QQ[x,y,z]
F={x^2*y*z + x^2*z^2 - y^3*z - y^3 }
determinantalUnitEuclideanDistanceDegree(F)--15
primaryDecomposition ideal singularLocus ideal F
CV=conormalVariety(ideal F);




R=QQ[x,y,z]
F={x^2*y*z + x^2*z^2 - y^3*z - y^3 }
determinantalUnitEuclideanDistanceDegree(F)--15
primaryDecomposition ideal singularLocus ideal F
CV=conormalVariety(ideal F);





eliminate(drop(gens ring CV,-3),CV)


R=QQ[x0,x1,x2,x3]
I=ideal sum apply(gens R,i->i^3);


R=QQ[x,y,z]
F={x^2*y*z + x^2*z^2 - y^3*z - y^3 }
determinantalUnitEuclideanDistanceDegree(F)--15
primaryDecomposition ideal singularLocus ideal F
CV=conormalVariety(ideal F);

R=QQ[y,z]
gens coefficientRing R



---ED degree Equation 3.6
R=QQ[s,t]
sectionalEDdegree=(m,n,d1,d2,j)->(    
    L:= 4*(1+t)^m*(1+s)^n*(t+s)^j*sum(apply(d1+1,i->(-2*t)^i))*sum(apply(d2+1,i->(-2*s)^i));
    L=diff(t^(m-2)*s^(n-2),L);
    L=sub(L,{t=>0,s=>0});
    L=1/((m-2)!)*1/((n-2)!)*L	  
	  )

(m,n,d1,d2,jjj)=(5,5,20,20,2)
sectionalEDdegree(m,n,d1,d2,jjj)



R=QQ[x,y,z,w]
F={det genericMatrix(R,2,2) -1}
determinantalUnitEuclideanDistanceDegree(F)
leftKernelGenericEDDegree(storeBM2Files,1,F)
runBertiniEDDegree(storeBM2Files)

leftKernelUnitEDDegree(storeBM2Files,1,F)
runBertiniEDDegree(storeBM2Files)
primaryDecomposition ideal singularLocus ideal F
CV=conormalVariety(ideal F);


---
nSize=4
R=QQ[a_(1,1)..a_(nSize,nSize)]
F={det genericMatrix(R,nSize,nSize) -1}
determinantalUnitEuclideanDistanceDegree(F)
nSize/2*2^nSize
determinantalGenericEuclideanDistanceDegree(F)

--Generic case: {8,120}
leftKernelGenericEDDegree(storeBM2Files,1,F)
runBertiniEDDegree(storeBM2Files)