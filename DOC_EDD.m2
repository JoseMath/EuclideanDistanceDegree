
doc /// --EuclideanDistanceDegree  
    Key
        EuclideanDistanceDegree 
    Headline
        a package to determine Euclidean distance degrees
    Description
      Text
        This package provides several routines for determining the Euclidean distance degree of an algebraic variety.
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
    	leftKernelUnitEDDegree(storeBM2Files,c,F)
	2==runBertiniEDDegree(storeBM2Files)
      Text
      	This package also computes generic ED degrees. The generic ED degree of X is always greater than or equal to the unit ED degree X.
      Example
    	R=QQ[x,y];
	F={x^2+y^2-1};	c=1
    	4==determinantalGenericEuclideanDistanceDegree(F)
    	leftKernelGenericEDDegree(storeBM2Files,c,F)
	4==runBertiniEDDegree(storeBM2Files)
      Text
      	The most general method for computing ED degrees with symbolic computation is symbolicWeightEDDegree
      Example
   	R=QQ[x,y];
	F={x^2+y^2-1};
	genericWeightVector={2,3}
	unitWeightVector={1,1}
	dataVector={5,7}
	4==symbolicWeightEDDegree(F,dataVector,genericWeightVector)
	2==symbolicWeightEDDegree(F,dataVector,unitWeightVector)		
      Text
      	When the variety is an affine cone, one is able to compute ED degrees using ED degree homotopies.
	The easiest case is when the variety is a hypersurface (or more generally, a complete intersection)  
      Example
        R=QQ[x1,x2,x3,x4]
	F={det genericMatrix(R,2,2)};
    	G=F--when F defines a complete intersection we may take F=G
    	P=(F,G)
	6==numericWeightEDDegree(storeBM2Files,P)
	2==numericUnitEDDegree(storeBM2Files,P)	
      Text
      	When an affine cone is not a complete intersection we use membership tests to compute ED Degrees.
	Here V(F) is an irreducible component of V(G) (a reducible variety) and #G===codim ideal F.
	These methods employ an equation by equation method called regeneration. 
      Example
        R=QQ[x1,x2,x3,x4,x5,x6]
	F=(minors(2,genericMatrix(R,3,2)))_*;
    	G=drop(F,-1);	
    	P=(F,G)
    	#G==codim ideal F;
	10==numericWeightEDDegree(storeBM2Files,P)
      Text
      	One may also determine (Unit) ED degrees using a parameter homotopy called a Weight-ED Degree Homotopy. 
      Example
      	printingPrecision=100
	R=QQ[x1,x2,x3,x4]
	F=(minors(2,genericMatrix(R,2,2)))_*;
    	G=drop(F,0);	
	L={}
	stageOne=1
    	stageTwo=2
	P=(F,G,L)
    	theDir=temporaryFileName()
	if not fileExists theDir then mkdir theDir
	NCO=newNumericalComputationOptions(theDir,P)
	TWV={1,1,1,1}
	NCO#"TargetWeight"=TWV 
	ht=(0,0,null)
	startEDDegree(NCO,ht,stageOne)
	runBertiniStartEDDegree(NCO,ht,stageOne,#F)
	readFile(NCO#"Directory","member_points",10000)
--
	filterSolutionFile(NCO,"start",#gens ring first F+1+#G+#L)     
	readFile(NCO#"Directory","start",100000)
--
	startEDDegree(NCO,ht,stageTwo)
	runBertiniStartEDDegree(NCO,ht,stageTwo,#F)
--(4,3) (3,1) (3,2)
        weightEDDegreeHomotopy(theDir,P,TWV)
        R=QQ[x1,x2,x3,x4,x5,x6]
	F=(minors(2,genericMatrix(R,3,2)))_*;
    	G=drop(F,-1);	
	L={}
	stageOne=1
    	stageTwo=2
	P=(F,G,L)
	NCO=newNumericalComputationOptions(theDir,P)
	TWV={1,1,1,1,1,1}
	NCO#"TargetWeight"=TWV
	startEDDegree(NCO,(0,0,null),stageOne)
	10==runBertiniStartEDDegree(NCO,#F)
	readFile(NCO#"Directory","member_points",10000)
	filterSolutionFile(NCO,"start",#gens ring first F+1+#G+#L)     
	readFile(NCO#"Directory","start",100000)
	startEDDegree(NCO,(0,0,null),stageTwo)
	runBertiniStartEDDegree(NCO,#F)--2



	
    	
///;
end
restart
loadPackage("EuclideanDistanceDegree",Reload=>true)
stageOne=1
stageTwo=2
P=(F,G,L)





restart
loadPackage"EuclideanDistanceDegree"
NCO#"StartWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,0,0),1)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)
NCO#"TrackSolutions"
peek NCO
   	printingPrecision=100
	F=flatten entries gens minors(2,transpose genericMatrix(R,3,2))
	G=drop(F,-1)
L={}
NCO=newNumericalComputationOptions(theDir,(F,G))
NCO#"StartWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,0,1),1)
runBertiniStartEDDegree(#F,NCO)
filterSolutionFile(theDir,"start",#gens ring first F+1+#G+#L,NCO)     
readFile("start",1000)

	F={det genericMatrix(R,2,2)};
	genericWeightVector={2,3,1,5}
	unitWeightVector={1,1,1,1}
	dataVector={5,1,11,13}
	
	4==symbolicWeightEDDegree(F,dataVector,genericWeightVector)
	2==symbolicWeightEDDegree(F,dataVector,unitWeightVector)		

restart
loadPackage("EuclideanDistanceDegree",Reload=>true)
help EuclideanDistanceDegree

doc ///--multiaffineDimension
 Key
   multiaffineDimension
   (multiaffineDimension,List,List,List,Point)--(F,E,L,pt)
   (multiaffineDimension,List,List,Point)--(F,L,pt)
   (multiaffineDimension,List,Sequence,Point) --(F,S,pt)
   (multiaffineDimension,List,Point) --(F,pt)
 Headline
   a method to determine multiaffine dimensions
 Usage
   P = multiaffineDimension(F,E,L,pt)
   P = multiaffineDimension(F,L,pt)
   P = multiaffineDimension(F,S,pt)
   P = multiaffineDimension(F,pt)
 Inputs
   F:List
     polynomials (system need not be square)
   E:List
     an index set for the affine factors 
   L:List
     entries in E and lenth equal to the number of indeterminants 
   S:Sequence
     a sequence of positive integers that sum to the number of indeterminants
   pt:Point
     a general point of an irreducible component of V(F)
 Outputs
   P:Polyhedron
     multidimension of the irreducible component of V(F) containing pt in terms of a polymatroid polytope, use latticePoints P to recover the multidimension explicitly. 
 Description
   Text
     multiaffineDimension(F,E,L,pt) is the most flexible way to use this function. 
   Example
     R = CC[x1,x2,x3,y];
     F = {x1-x2-y};
     pt = point{{2_CC,1_CC,1_CC,3_CC}};
     (E,L)=({0,1},{0,0,0,1})
     P = multiaffineDimension(F,E,L,pt)
     latticePoints P
   Text
     multiaffineDimension(F,L,pt) calls multiaffineDimension(F,E,L,pt) with  E=toList set L 
   Example
     R = CC[x1,x2,x3,y];
     F = {x1-x2-y};
     pt = point{{2_CC,1_CC,1_CC,3_CC}};
     L={0,0,0,1}
     P = multiaffineDimension(F,L,pt)
     Q = multiaffineDimension(F,{x,x,x,y},pt)
     Q==P
   Text
     multiaffineDimension(F,S,pt) calls multiaffineDimension(F,E,L,pt) with E={0,..,#S-1} and L the first S_0 elements of L are 0, the next S_1 elements are 1, and so on, and then 
   Example
     R = CC[x1,x2,x3,y];
     F = {x1-x2-y};
     pt = point{{2_CC,1_CC,1_CC,3_CC}};
     S = (3,1)
     P = multiaffineDimension(F,S,pt)
     (E,L)=({0,1},{0,0,0,1})
     Q = multiaffineDimension(F,E,L,pt)
     Q==P
   Text
     multiaffineDimension(F,pt) calls multiaffineDimension(F,E,L,pt) with E=L={0,..,# variables -1} 
   Example
     R = CC[x1,x2,x3,y];
     F = {x1-x2-y};
     pt = point{{2_CC,1_CC,1_CC,3_CC}};
     P = multiaffineDimension(F,pt)
     E=L={0,1,2,3}
     Q = multiaffineDimension(F,E,L,pt)
     Q==P
 Caveat
   pt is assumed to be a generic point and the numerical rank is computed correctly
      
///;


doc /// --restrictionCodimension
 Key
   restrictionCodimension
   (restrictionCodimension,Polyhedron,List)
   (restrictionCodimension,Polyhedron)
   (restrictionCodimension,Polyhedron,Nothing)
   RestrictionInformation
 Headline
   a method to restrict the dimension polytope
 Usage
   ri = restrictionCodimension(P,s)
   ri = restrictionCodimension(P)
   ri = restrictionCodimension(P,)
 Inputs
   P:Polyhedron
      encodes the dimension information.
   s:List
     a permutation, i.e., a length #affine factors-1 list with distinct integers from 0 to #affine factors-1
 Outputs
   ri:RestrictionInformation
     a mutable hashtable with keys "CodimensionRestriction", "DimensionReduction", and "Permutation".
     The respective values of these keys are m, Q, and s.
     These are the outputs from Algorithm ?? in arxiv: TBD.
 Description
   Text
     restrictionCodimension(P,s) is the most flexible way to use this function. 
   Example
     P=convexHull(matrix transpose {{8,1,9},{9,9,0},{2,6,10}})
     latticePoints P
     ri=restrictionCodimension(P,{0,1,2})
     ri#"CodimensionRestriction"
     Q=ri#"DimensionReduction"
     Q//latticePoints
   Text
     restrictionCodimension(P) calls restrictionCodimension(P,s) with s={0,1,...,#variables -1} 
   Example
     P=convexHull(matrix transpose {{8,1,9},{9,9,0},{2,6,10}})
     ri1=restrictionCodimension(P,{0,1,2})
     Q1=ri1#"DimensionReduction"
     ri2=restrictionCodimension(P)
     Q2=ri2#"DimensionReduction"
     Q1==Q2      
   Text
     restrictionCodimension(P,) calls restrictionCodimension(P,s) with s={#variables -1,...,1,0} 
   Example
     P=convexHull(matrix transpose {{8,1,9},{9,9,0},{2,6,10}})
     ri1=restrictionCodimension(P,{2,1,0})
     Q1=ri1#"DimensionReduction"
     ri2=restrictionCodimension(P,)
     Q2=ri2#"DimensionReduction"
     Q1==Q2      
///;


doc /// --sequentialDimension
 Key
   sequentialDimension
   (sequentialDimension,Polyhedron)
   (sequentialDimension,Polyhedron,List)
   (sequentialDimension,Polyhedron,Nothing)
 Headline
   a method to determine how to slice the variety
 Usage
   LT = sequentialDimension(P,s)
   LT = sequentialDimension(P)
   LT = sequentialDimension(P,)
 Inputs
   P:Polyhedron
      encodes the dimension information.
   s:List
     a permutation, i.e., a length #affine factors-1 list with distinct integers from 0 to #affine factors-1
 Outputs
   LT:LinearSpaceType
     every entry is a 0/1 list with length equal to the number of affine factors
 Description
   Text
     restrictionCodimension(P,s) is the most flexible way to use this function. 
   Example
     P=convexHull(matrix transpose {{2,1,0},{0,1,2}});
     m1 =sequentialDimension(P,{0,1,2})    
     m2 =sequentialDimension(P,{2,1,0})    
     m3 =sequentialDimension(P,{1,0,2})    
     m1 == sequentialDimension(P)    
     m2 == sequentialDimension(P,)    

///;


doc /// --wellPositionedLinearSpace
 Key
   wellPositionedLinearSpace
   (wellPositionedLinearSpace,Polyhedron,List,List)
   (wellPositionedLinearSpace,Polyhedron)
   (wellPositionedLinearSpace,Polyhedron,Nothing)
 Headline
   a method that returns the type of linear polynomials that define a well-positioned linear space
 Usage
   LT = wellPositionedLinearSpace(P,r,s)
   LT = wellPositionedLinearSpace(P)
   LT = wellPositionedLinearSpace(P,)
 Inputs
   P:Polyhedron
      encodes the dimension information.
   r:List
     a permutation, i.e., a length #affine factors-1 list with distinct integers from 0 to #affine factors-1
   s:List
     a permutation, i.e., a length #affine factors-1 list with distinct integers from 0 to #affine factors-1
 Outputs
   LT:LinearSpaceType
     every entry is a 0/1 list with length equal to the number of affine factors
 Description
   Text
     wellPositionedLinearSpace(P,r,s) is the most flexible way to use this function. 
   Example
     P=convexHull(matrix transpose {{1,2,1},{0,2,2}});
     r={0,1,2}
     s={2,1,0}
     LT=wellPositionedLinearSpace(P,r,s)
     ri=(restrictionCodimension(P,r))         
     LT==(ri#"RestrictionType"|sequentialDimension(ri#"DimensionReduction",s))	 
   Text
     wellPositionedLinearSpace(P) calls wellPositionedLinearSpace(P,r,s) with r=s={0,1,..#affine factors -1} 
   Example
     P=convexHull(matrix transpose {{1,2,1},{0,2,2}});
     r={0,1,2}
     s={0,1,2}
     wellPositionedLinearSpace(P,r,s)==wellPositionedLinearSpace(P)
   Text
     wellPositionedLinearSpace(P,) calls wellPositionedLinearSpace(P,r,s) with r=s={#affine factors -1,...,1,0} 
   Example
     P=convexHull(matrix transpose {{1,2,1},{0,2,2}});
     r={2,1,0}
     s={2,1,0}
     wellPositionedLinearSpace(P,r,s)==wellPositionedLinearSpace(P,)

///;


doc /// --wellPositionedIdeal
 Key
   wellPositionedIdeal
   (wellPositionedIdeal,Ring,List,List,LinearSpaceType) 
 Headline
   a method that returns an ideal of a well-positioned linear space
 Usage
   LT = wellPositionedIdeal(R,E,L,LT)
   LT = wellPositionedIdeal(R,L,LT)
   LT = wellPositionedIdeal(R,S,LT)
   LT = wellPositionedIdeal(R,LT)
 Inputs
   P:Polyhedron
      encodes the dimension information.
   E:List
     an index set for the affine factors 
   L:List
     entries in E and lenth equal to the number of indeterminants 
   S:Sequence
     a sequence of positive integers that sum to the number of indeterminants
   LT:LinearSpaceType
     every entry is a 0/1 list with length equal to the number of affine factors
 Outputs
   I:Ideal
     an ideal of a well-positioned linear space
 Description
   Text
     wellPositionedIdeal(R,E,L,LT) is the most flexible way to use this function. 
   Example
     R=QQ[x1,x2,y,z]
     E={x,y,z}
     L={x,x,y,z}
     LT=new LinearSpaceType from {{1,0,0},{1,1,1}}
     wellPositionedIdeal(R,E,L,LT)    
   Text
     wellPositionedIdeal(R,L,LT) calls wellPositionedIdeal(R,E,L,lT) with  E=toList set L 
   Example
     R=QQ[x1,x2,y,z]
     L={x,x,y,z}
     LT=new LinearSpaceType from {{1,0,0},{1,1,1}}
     wellPositionedIdeal(R,L,LT)    
   Text
     wellPositionedIdeal(R,S,LT) calls wellPositionedIdeal(R,E,L,LT) with E={0,..,#S-1} and L the first S_0 elements of L are 0, the next S_1 elements are 1, and so on, and then 
   Example
     R=QQ[x1,x2,y,z]
     S=(2,1,1)
     LT=new LinearSpaceType from {{1,0,0},{1,1,1}}
     wellPositionedIdeal(R,S,LT)    
   Text
     wellPositionedIdeal(R,LT) calls wellPositionedIdeal(R,E,L,LT) with E=L={0,..,# variables -1} 
   Example
     R=QQ[x1,x2,y,z]
     LT=new LinearSpaceType from {{1,1,0,0},{1,1,1,1}}
     wellPositionedIdeal(R,LT)    

///;


doc/// --LinearSpaceType
 Key
   LinearSpaceType
 Headline
   a new type of List that specifies the supports of linear polynomials  
 Description
   Text
     An entry L of a LinearSpaceType is a list. 
     This list has 0/1 entries and the length of L is the number of affine factors of the ambient space for the variety V(F).
     The list L also induces a linear polynomial that is a sum of general linear polynomials of CC[X_i] where i varies over the positions of L with a 1.  
   Example
     R=QQ[x1,x2,y,z]
     L={x,x,y,z}
     LT=new LinearSpaceType from {{1,0,0},{1,0,0},{1,1,1}}
     wellPositionedIdeal(R,L,LT)    
     
///;



--Examples from the paper.
--path=prepend("/Users/jo/Documents/GoodGit/MultiaffineDimension",path)
--loadPackage("MultiaffineDimension",Reload=>true)
--needsPackage("Bertini")
