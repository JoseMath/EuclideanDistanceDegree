--restart
--Projective formulation for intersections with linear spaces
rand:=randomValue
--Assume ring is a complex inexact field
--G is a subset of F. 
NumericalComputationOptions=new Type of MutableHashTable
(stageOne,stageTwo)=(1,2); 
  
    
    
-*
restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
check EuclideanDistanceDegree
R=QQ[x1,x2,x3,x4]
M=matrix{{x1,x2,x3},{x2,x3,x4}}
F=(minors(2,M))_*
G=drop(F,-1)--7
L={}
P=(F,G,L)
sf=temporaryFileName();mkdir sf
NCO=newNumericalComputationOptions(sf,P)
--determinantalUnitEuclideanDistanceDegree(F)
--determinantalGenericEuclideanDistanceDegree(F)
homotopyEDDegree(NCO,"Weight",true,true)
*-

    
parameterKeys={    "StartWeight",    "TargetWeight",
    "StartData",    "TargetData", 
    "GammaVector"}

jacKeys={    "JacobianWitnessModel","JacobianStartSubmodel","JacobianTargetSubmodel"}
modelKeys={    "Model","WitnessModel",    "StartSubmodel", "TargetSubmodel"}
degreeKeys={ "DegreeWitnessModel","DegreeSubmodel"}
--
bertiniKeys={    "BertiniStartFiberSolveConfiguration","BertiniMembershipTestConfiguration",    "BertiniSubstitute","BertiniConstants"}
coordinateKeys={    "PrimalCoordinates",    "HomogeneousVariableGroups",    "AffineVariableGroups"    }
directoryKeys={"Directory"} 
solutionKeys={"TrackSolutions"}
fixValues={"FixedData","FixedWeight","FixedSubmodel","FixedJacobianSubmodel"}
nocKeys=parameterKeys|jacKeys|modelKeys|degreeKeys|bertiniKeys|coordinateKeys|directoryKeys|solutionKeys|fixValues

defaultFixValues={"StartData","StartWeight","StartSubmodel","JacobianStartSubmodel"}



newNumericalComputationOptions=method()
newNumericalComputationOptions(String,Sequence):=(theDir,P)->(
    if #P==3 then (F,G,L):=P;
    if #P==2 then (F,G,L)=(P_0,P_1,{});
    NCO:=new NumericalComputationOptions from apply(nocKeys,i->i=>null);
    NCO#"Directory"=theDir;----directory where the files will be stored.
    NCO#"Model"=F,
    NCO#"WitnessModel"=G;
    NCO#"StartSubmodel"=L;
    NCO#"TargetSubmodel"=L;
    NCO#"DegreeSubmodel"=L/degree/first;
    NCO#"DegreeWitnessModel"=G/degree/first;
    NCO#"JacobianWitnessModel"=diff(matrix{gens ring first G}, transpose matrix{G} );
    NCO#"JacobianTargetSubmodel"=NCO#"JacobianStartSubmodel"=diff(matrix{gens ring first G}, transpose matrix {L} );
    NCO#"PrimalCoordinates"=gens ring first F; ---This is different when working with a parameterization
    numX:=#gens ring first G;
    NCO#"TargetData"=NCO#"StartData"=apply(numX,i->random CC); 
    NCO#"TargetWeight"=apply(numX,i->1);
    NCO#"StartWeight"=apply(numX,i->random CC); 
--    NCO#"GammaVector"=apply(numX-1,i->random CC); 
    scan(bertiniKeys,i->NCO#i={});
    NCO#"HomogeneousVariableGroups"={gens ring first F};
    NCO#"AffineVariableGroups"={};    
    NCO#"FixedData"="StartData";
    NCO#"FixedWeight"="StartWeight";
    NCO#"FixedSubmodel"="StartSubmodel";
    NCO#"FixedJacobianSubmodel"="JacobianStartSubmodel";
    NCO#"JacobianVars"="jv";
    NCO#"GradientVars"="gv";
    NCO#"ScaleVars"="sv";
    NCO#"DataVars"="dv";
    NCO#"WeightVars"="wv";
    NCO#"Infinity"=null;
    NCO#"PairGeneralHyperplaneList"=null;
    return NCO
    )

defaultFixValues={"StartData","StartWeight","StartSubmodel","JacobianStartSubmodel"}
-*
restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x,y,z,w]
F=G={det genericMatrix(R,2,2)}
L={}
P=(F,G,L)
sf=temporaryFileName();mkdir sf
NCO=newNumericalComputationOptions(sf,P)
homotopyEDDegree(NCO,"Weight",true,true)--2
*-

-*
restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x,y,z,w]

F=G={x^3-y*w+z+1,random({1},R)}
L={}
P=(F,G,L)
sf=temporaryFileName();mkdir sf
NCO=newNumericalComputationOptions(sf,P)
homotopyEDDegree(NCO,"Weight",true,true)
determinantalUnitEuclideanDistanceDegree(F)
*-

-*
restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x,y,z,w]
F=G={x^3-y*w+z+1,x*y-z*w}
L={}
P=(F,G,L)
sf=temporaryFileName();mkdir sf
NCO=newNumericalComputationOptions(sf,P)
determinantalUnitEuclideanDistanceDegree(F)
determinantalGenericEuclideanDistanceDegree(F)
homotopyEDDegree(NCO,"Weight",true,true)
*-

-*
restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x,y,z,w]
F=G={random({1},R),random({1},R)}
L={}
P=(F,G,L)
sf=temporaryFileName();mkdir sf
NCO=newNumericalComputationOptions(sf,P)
determinantalUnitEuclideanDistanceDegree(F)
determinantalGenericEuclideanDistanceDegree(F)
homotopyEDDegree(NCO,"Weight",true,true)
*-


-*
restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x,y,z,w]
F=G={x+y,z-w}
L={}
P=(F,G,L)
sf=temporaryFileName();mkdir sf
NCO=newNumericalComputationOptions(sf,P)
determinantalUnitEuclideanDistanceDegree(F)
determinantalGenericEuclideanDistanceDegree(F)
homotopyEDDegree(NCO,"Weight",true,true)
*-



-*
restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x1,x2,x3,x4,x5,x6,x7,x8,x9]
M=genericMatrix(R,3,3)
F=(minors(2,M))_*
G={F_0,F_1,F_4,F_5}+2*{F_1,F_2,F_3,F_7}-3*{F_2,F_1,F_0,F_5}+2*{F_2,F_1,F_6,F_8}
G={F_0,F_1,F_4,F_5}
4==codim first decompose ideal G
#G
L={}
P=(F,G,L)
sf=temporaryFileName();mkdir sf
NCO=newNumericalComputationOptions(sf,P)
--determinantalUnitEuclideanDistanceDegree(F)
--determinantalGenericEuclideanDistanceDegree(F)
homotopyEDDegree(NCO,"Weight",true,true)
*-



-*
R=QQ[x1,x2,x3,x4]
M=matrix{{x1,x2,x3},{x2,x3,x4}}
F=(minors(2,M))_*
G=drop(F,-1)
codim\decompose ideal G

S=QQ[gens R|{l0,l1,l2,HX}]
CM=matrix{apply(gens R,i->sub(i,S)-HX*random(1,1000))}||sub(matrix makeJac(G,gens R),S)
critG=ideal(matrix{{l0,l1,l2}}*CM)+sub(ideal(G),S)
critF=ideal(matrix{{l0,l1,l2}}*CM)+sub(ideal(F),S)
netList decompose critG
netList apply(decompose critG,i->ideal mingens (sub(ideal F,S)+saturate(i,HX*l0)))

help saturate
decWIN=decompose(minors(3,CM)+sub(ideal G,S))
decWIN/degree
decWIN/codim

decompose (decWIN_1+sub(ideal(F),S))
oo/codim
topRow=matrix{1431-x1}
*-



homotopyEDDegree=method()
possibleHT={"Weight","Data","Submodel"}
stageOne=1
stageTwo=2
--ht="Weight"
--isStageOne=true
--isStageTwo=true
homotopyEDDegree(NumericalComputationOptions,String,Boolean,Boolean):=(NCO,ht,isStageOne,isStageTwo)->(    
--(CODE 1) First we set the type of homotopy that will be performed.
    if not member(ht,possibleHT) then error("Argument 1 is in "|toString possibleHT);
    if ht===possibleHT_0 then     weightHomotopy:=true else weightHomotopy=false ;    
    if ht===possibleHT_1 then     dataHomotopy:=true else dataHomotopy=false ;
    if ht===possibleHT_2 then     submodelHomotopy:=true else submodelHomotopy=false ;
--(CODE 2) Extract information from NCO (NumericalComputationOptions).
--The homotopy is highly customizable, which is why we use a new type of MutableHashTable.
    jacL:=NCO#(NCO#"FixedJacobianSubmodel");
    L:=NCO#(NCO#"FixedSubmodel"); 
    startWeight:=  NCO#"StartWeight";
    targetWeight:=NCO#"TargetWeight";
    startData:=  NCO#"StartData";
    targetData:=NCO#"TargetData";
    (lagMult,numerHB,denomQ,primal,tWeight):=("lagMult","numerHB","denomQ","primal","tWeight")   ; 
    (jv,gv,sv):=(NCO#"JacobianVars",    NCO#"GradientVars",    NCO#"ScaleVars");   
    (dv,wv):=(NCO#"DataVars",    NCO#"WeightVars");   
--F is the model, V(G)\cap V(L) is a complete intersection contained in V(F)\cap V(L).
    (F,G,startL,targetL,jacG):=(NCO#"Model",NCO#"WitnessModel",NCO#"StartSubmodel",NCO#"TargetSubmodel",NCO#"JacobianWitnessModel");
--    randomGamma:=NCO#"GammaVector";
    xList:=NCO#"PrimalCoordinates";
--Set hyperplane at infinity
    if NCO#"Infinity"===null 
    then NCO#"Infinity"=makeB'Section(xList,NameB'Section=>"HX");
--(FUNCTION 0) --quicjly create variabels
    vRing:=(numVars,s,kk)->kk[apply(numVars,i->s|i)];   
--(CODE 3) Now we make rings.
    nc:=#xList;
    kk0:=QQ; 
    extrinsicRing:=kk0[flatten transpose apply(#G+#L,i->apply(nc,j->jv|i|"v"|j))];
    scan({sv,gv,dv,wv},{#G+#L+1,nc,nc,nc},(s,numVars)->extrinsicRing=extrinsicRing**vRing(numVars,s,kk0));
    symbolicJac:=genericMatrix(extrinsicRing,#G+#L,nc);
    symbolicScaleMatrix:=basis({0,1,0,0,0},extrinsicRing);
    symbolicGradient:=basis({0,0,1,0,0},extrinsicRing);
-- symbolicSystem is the system we want to solve after subsituting subfunctions.
    symbolicSystem:=symbolicScaleMatrix*(symbolicGradient||symbolicJac);
    jacZero:={};
    pairJac:={};
-- (FUNCTION 1) ---pair a row of matrix with values
    pairRowFunction:=(M1,M2,hom)->(	
	arg:=flatten entries M1;
	val:=flatten entries M2;
	scan(#arg,i->if val_i==0 
	    then (jacZero=jacZero|{arg_i=>0};
		symbolicSystem=sub(symbolicSystem,{arg_i=>0}))
	    else pairJac=pairJac|{makeB'Section({val_i},
    	    	B'NumberCoefficients=>{1},	
	    	B'Homogenization=>hom,	
	    	NameB'Section=>arg_i )}));
 --   M1=matrix{{x,y},{z,w}}
--    M2=matrix{{1,2},{5,4}}
--   print peek last pairRow(M1,M2,1,"HX")
-- (FUNCTION 2) ---pair parameters of a parameter homotopy
    pairParameterFunction:=(p0,p1,r1,r2,sym,bool)->(
	if bool then pp:=apply(#p0,i->makeB'Section({p0_i,p1_i},
    	    B'NumberCoefficients=>{r1,r2},		
	    NameB'Section=>sym_i )) else
    if not bool then  pp=apply(#p0,i->makeB'Section({p0_i},
    	    B'NumberCoefficients=>{1},		
	    NameB'Section=>sym_i
	    )));
---(CODE 4) Now set up subfunctions. This is done by pairing a symbol with a value by an option or B'Section.
    weight:=symbolicWeight := gens vRing(nc,wv,kk0);
    data:=symbolicData := gens vRing(nc,dv,kk0);
    pairData:=pairParameterFunction(startData,targetData,"(1-TData)","TData",symbolicData,dataHomotopy);	    
    pairWeight:=pairParameterFunction(startWeight,targetWeight,"(1-TWeight)","TWeight",symbolicWeight,weightHomotopy);
    print ("pairData",peek first  pairData    );
    print ("pairWeight",peek first  pairWeight    );
---(CODE 5) Pair Gradient Vector 
    kk2:=ring first startWeight;
    topS:=kk2[numerHB,denomQ,tWeight];
    (topNumerHB,topDenomQ,topTWeight):=toSequence flatten entries basis({1},topS);
    pairGradient:=apply(#xList,i->makeB'Section({xList_i},
	    ContainsPoint=>{data_i},
	    B'NumberCoefficients=>{weight_i},
	    NameB'Section=>symbolicGradient_(0,i),
	    B'Homogenization=>"HX"));
---(CODE 6) Pair Jacobian: 
--create ring to homogenize rows (indexed by polynomials) of Jacobian 
    jacLG:=jacL||jacG;
    kk3:=coefficientRing ring first F;
    jacRing:=kk3[gens ring first F|{"HX"}];
    HX:=last gens jacRing;
    homogLG:= homogenize(sub(matrix{L|G},jacRing),HX)//entries//flatten;
    homogJac:=    matrix apply(numrows jacLG,i->apply(numcols jacLG,j->diff((gens jacRing)_j,homogLG_i)));
    print homogJac;
    pairRowFunction(symbolicJac,homogJac,"HX");
---(CODE 7) Pair Scaling variables (Lagrange multipliers): 
--Determine degrees to properly homogenize cols (indexed by variables) of Jacobian
--#pairJac;
--#jacZero;
    (degSubmodel,degWitnessModel):=(NCO#"DegreeSubmodel",NCO#"DegreeWitnessModel");
    degAugJac:={1}|apply(degSubmodel|degWitnessModel,i->i-1);
    maxDegree:=degAugJac//max;
    degRescale:=degAugJac/(i->maxDegree-i);        
    bLagrangeVars:=lagList:=apply(#degRescale,i->"L"|i);
    rescaling:=new MutableList from apply(#degRescale,i->lagList_i);
--Homogenize cols by multiplying by a diagonal matrix of linear products on the left. 
----The following determines these linear products. 
    generalHyperplaneList:={};
    scan(#degRescale,i->scan(
	    degRescale_i, 
	    j->(hg:="H"|i|"v"|j;---wants to be both
		rescaling#i=(rescaling#i)|"*"|hg;
	    	generalHyperplaneList=generalHyperplaneList|{hg})
	    ));
--    print(peek rescaling);
    if NCO#"PairGeneralHyperplaneList"=!=null then 
    pairGeneralHyperplanes:=NCO#"PairGeneralHyperplaneList"  else(
	pairGeneralHyperplanes=apply(#generalHyperplaneList,i->
	    makeB'Section(xList|{"HX"},
		NameB'Section=>generalHyperplaneList_i));
	NCO#"PairGeneralHyperplaneList"=pairGeneralHyperplanes);
--    print(peek first pairGeneralHyperplanes);
    pairScale:=apply(flatten entries symbolicScaleMatrix,rescaling,(i,j)->i=>j);
-- (CODE 8)  Set up inputs for bertini. 
    bModelVars:=gens ring first F|{"HX"}   ;
    bPoly:=homogLG|flatten entries symbolicSystem;
    bConfiguration:={"UseRegeneration"=>1,
	"TrackType"=>0,
	"PrintPathProgress"=>1000}|(NCO#"BertiniStartFiberSolveConfiguration");    
    BF:=pairData|pairWeight|pairJac|pairGradient|pairGeneralHyperplanes|pairScale;
-- (FUNCTIONS 2) Functions for solving (write input)
    writeSolveInputFunction:=(stage,nif)->(
	if stage===stageOne then (PG:={"adfadfdf"}; BC:={"TData"=>0,"TWeight"=>0})
	else if stage===stageTwo
	then (BC={};
	    if dataHomotopy then PG={"TData"}
	    else if weightHomotopy then PG={"TWeight"});
    	makeB'InputFile(NCO#"Directory",
    	    NameB'InputFile=>nif,
	    HomVariableGroup=>{bLagrangeVars,bModelVars},
    	    B'Configs=>bConfiguration|{"ParameterHomotopy"=>stage},
	    B'Polynomials=>bPoly,
	    B'Constants=>jacZero|BC,
	    ParameterGroup=>PG,
    	    B'Functions=>BF
	    ));
-- (FUNCTIONS 2) Functions for solving (run file input)
--our solution file will always be member points. 
    criticalPointName:="criticalPointFile";
    runSolveInputFunction:=(stage,nif)->(
	writeSolveInputFunction(stage,nif); 
    	if stage==stageTwo then(
    	    writeParameterFile(NCO#"Directory",{0},NameParameterFile=>"start_parameters");
    	    writeParameterFile(NCO#"Directory",{1},NameParameterFile=>"final_parameters"));
	runBertini(NCO#"Directory",NameB'InputFile=>nif);
    	readFile(NCO#"Directory");
    	if stage==stageOne then(	
	    moveB'File(NCO#"Directory","nonsingular_solutions","stageOne_solutions",CopyB'File=>true);
	    moveB'File(NCO#"Directory","nonsingular_solutions","start",CopyB'File=>true);
	    moveB'File(NCO#"Directory","nonsingular_solutions","member_points",CopyB'File=>true);
	    moveB'File(NCO#"Directory","nonsingular_solutions",criticalPointName,CopyB'File=>true));
    	if stage==stageTwo then(
	    moveB'File(NCO#"Directory","nonsingular_solutions","stageTwo_solutions",CopyB'File=>true);
	    --moveB'File(NCO#"Directory","nonsingular_solutions","start",CopyB'File=>true);
	    moveB'File(NCO#"Directory","nonsingular_solutions","member_points",CopyB'File=>true);
	    moveB'File(NCO#"Directory","nonsingular_solutions",criticalPointName,CopyB'File=>true));	    
	    );
-- (FUNCTIONS 3) Functions for membership test and returning incidence matrix (IM).
    ttOne:=1;
    ttThree:=3;    
    nameFileFunction:=(stage,case,indexCase,hypersurface,theTrackType)->("input_first_MT_"|case|"_"|indexCase|"_"|theTrackType);
    writeIsMembershipFunction:=(stage,case,indexCase,hypersurface,theTrackType)->(
	nif:=nameFileFunction(stage,case,indexCase,hypersurface,theTrackType);
    	if stage===stageOne then BC:={"TData"=>0,"TWeight"=>0};
    	if stage===stageTwo then BC={"TData"=>1,"TWeight"=>1};
    	if not member(stage,{1,2}) then error"stage is in {1,2}";
    	makeB'InputFile(NCO#"Directory",
    	    NameB'InputFile=>nif,
	    AffVariableGroup=>flatten flatten {bLagrangeVars,bModelVars},
    	    B'Configs=>bConfiguration|{"TrackType"=>theTrackType},
	    B'Polynomials=>{hypersurface},
	    B'Constants=>jacZero|BC,
--	    ParameterGroup=>PG,
    	    B'Functions=>BF
	    ));
--    isMembershipFunction(stageOne,"TT",0,"x1*x2-x3*x4")--Test!
    isMembershipFunction:=(stage,case,indexCase,hypersurface)->(
	--Pos dim solve TrackType=>1
	writeIsMembershipFunction(stage,case,indexCase,hypersurface,ttOne);
	nif:=nameFileFunction(stage,case,indexCase,hypersurface,ttOne);
	runBertini(NCO#"Directory",NameB'InputFile=>nif);
    	moveB'File(NCO#"Directory","bertini_session.log","bertini_session_"|nif|".log",CopyB'File => false);
--    	print nif;
	--MT TrackType=>3
	writeIsMembershipFunction(stage,case,indexCase,hypersurface,ttThree);
	nif=nameFileFunction(stage,case,indexCase,hypersurface,ttThree);
	runBertini(NCO#"Directory",NameB'InputFile=>nif);
    	moveB'File(NCO#"Directory","bertini_session.log","bertini_session_"|nif|".log",CopyB'File => false);
       	outIM:=importIncidenceMatrix(NCO#"Directory");
    	print nif;	
	print outIM;
	return outIM	)  ;  
-- (FUNCTIONS 4) Functions for filtering based off of incidence matrix
    filterSolutionFunction:=(nsf,kp,ns,numCoords)->(     
    	print("RUN FILTER",kp=>numCoords);    
    	firstLine := true;
    	countSol  := 0;
    	countLine := 0;
    	groupSize := 1+numCoords;
    	isSelected:= null;
    	sf:=openOut(NCO#"Directory"|"/"|nsf);
    	scanLineSolutionFunction := (ell)->(
      	    if firstLine 
      	    then (firstLine=false; sf<< toString(#kp)<<endl)
      	    else if countSol < ns
      	    then (
    	  	if countLine==0 then isSelected=member(countSol,kp);
	  	countLine=countLine+1;
    	  	if isSelected then sf <<ell<<endl;
      	  	if countLine==groupSize 
      	  	then (
    	    	    print(countSol=>isSelected);    	    	
	      	    --print (countLine,groupSize,"grp");
	      	    countLine=0; 
	      	    countSol=countSol+1;
	      	    )));
      scanLines(scanLineSolutionFunction,(NCO#"Directory")|"/"|"member_points");      
      close sf;
      return (nsf));
  --filterSolutionFunction("T1",{1,2,3,4,5,6,7},8)
--      saturateFunction=positionFunction=positionFilterFunction;
    positionFilterFunction:=(stage,case,indexCase,hypersurface,bin)->(--(stage,case,indexCase,hypersurface)
	isMembershipFunction(stage,case,indexCase,hypersurface);      
--    	(kp,ns):=positionMembershipFunction(stage,case,indexCase,hypersurface);
    	if bin==="typeA" 
	then isOffHypersurface:=(m->(m==={}))
    	else if bin==="typeB" 
	then isOffHypersurface=(m->(m=!={}))
	else error"last argument is typeA or typeB";
	imMT:=isMembershipFunction(stage,case,indexCase,hypersurface);
    	kp:={};
    	scan(#imMT,i->if isOffHypersurface(imMT_i) then kp=kp|{i});
    	ns:= #imMT;
    	print("kp",kp,"num kp",#kp,"num sols",ns);
	(nsf,nc):=("filterFile",#flatten {bLagrangeVars,bModelVars});
	filterSolutionFunction("filterFile",kp,ns,nc);
    	moveB'File(NCO#"Directory","filterFile","member_points",CopyB'File=>true);
	return #kp
	);
    stageEDDegBound:=new MutableList from {"empty",null,null};   
-- (FUNCTIONS 5) Functions to iterate filtering
    runSaturateUnionFunction:=(polyList,stage)->(
    	(case,bin):=("SaturateH","typeA");
    	scan(#polyList,i->(
		stageEDDegBound#stage=positionFilterFunction(stageOne,case,i,polyList_i,bin);
	    	print(peek stageEDDegBound,case,bin,polyList_i)));	 
    	print(peek stageEDDegBound));	        
-- (FUNCTIONS 6) Functions to restrict to the variety 
    runRestrictIntersectionFunction:=(polyList,stage)->(
    	(case,bin):=("IntersectF","typeB");
    	scan(#polyList,i->(
		print(peek stageEDDegBound,case,bin,polyList_i);
		stageEDDegBound#stage=positionFilterFunction(stageOne,case,i,polyList_i,bin)));	 
    	print(peek stageEDDegBound));	        
-- (Function 7) 
    runComputationStage:=(stage,offPolyList,onPolyList)->(
	if stage==stageOne then
	runSolveInputFunction(stageOne,"input_first_solve") else 
	runSolveInputFunction(stageTwo,"input_second_solve");
	print("offPolyList",offPolyList);
	runSaturateUnionFunction(offPolyList,stage);
    	print("WIN","SATURATE");
--    	moveB'File(NCO#"Directory","member_points","filterFile",CopyB'File=>true);
	print("onPolyList",onPolyList);
	runRestrictIntersectionFunction(onPolyList,stage);
	print("WIN","RESTRICT");
	print("WIN",stage)	);    
    offPolyList:={HX,"L0"}|((pairGeneralHyperplanes/(i->i#NameB'Section)));
    onPolyList :=F/(i->homogenize(sub(i,jacRing),HX));
    if isStageOne then runComputationStage(stageOne,offPolyList,onPolyList);
    if isStageTwo then runComputationStage(stageTwo,offPolyList,onPolyList);
    if isStageTwo then return stageEDDegBound#2 else if isStageOne then return stageEDDegBound#1
      )





numericWeightEDDegree=method()

numericWeightEDDegree(String,Sequence,List):=(dir,P,wv)->(    
    NCO:=newNumericalComputationOptions(dir,P);
--    WV:=apply(#gens ring first first P,i->random CC);
    NCO#"StartWeight"=wv;
    ht:="Weight";
    isStageOne:=true;
    isStageTwo:=false;
    homotopyEDDegree(NCO,ht,isStageOne,isStageTwo)
    )

numericWeightEDDegree(Sequence,List):=(P,wv)->(
    dir:=temporaryFileName();    
    if not fileExists dir then mkdir dir;
    numericWeightEDDegree(dir,P,wv)    )
--numericWeightEDDegree(P,{1,2,3,4})

edTypes={"Generic","Unit"}
numericEDDegree=method()
numericEDDegree(Sequence,String):=(P,typeED)->(    
    if typeED ==="Generic" then  wv:=apply(#gens ring first first P,i->random CC)
    else if typeED==="Unit" then wv=apply(#gens ring first first P,i->1)
    else error("last argument needs to be in "|toString edTypes);    
    numericWeightEDDegree(P,wv));    
--numericEDDegree(P,"Generic")
--numericEDDegree(P,"Unit")








































-*
peek first pairWeight
peek first pairData
    stage=stageTwo;
    if stage>=stageOne then (
    	writeSolveInputFile(stageOne); 
    	runBertini(NCO#"Directory",NameB'InputFile=>);
    	readFile(NCO#"Directory");
	moveB'File(NCO#"Directory","nonsingular_solutions","stageOne_solutions",CopyB'File=>true);
	moveB'File(NCO#"Directory","nonsingular_solutions","start",CopyB'File=>true);
	moveB'File(NCO#"Directory","nonsingular_solutions","member_points",CopyB'File=>true));
    readFile(NCO#"Directory");
    if stage==stageTwo then(
    	writeSolveInputFile(stageTwo) 	;
    	writeParameterFile(NCO#"Directory",{0},NameParameterFile=>"start_parameters");
    	writeParameterFile(NCO#"Directory",{1},NameParameterFile=>"final_parameters");
    	runBertini(NCO#"Directory",NameB'InputFile=>"input_second_solve");
	moveB'File(NCO#"Directory","nonsingular_solutions","stageTwo_solutions",CopyB'File=>true);
	moveB'File(NCO#"Directory","nonsingular_solutions","start",CopyB'File=>true);
	moveB'File(NCO#"Directory","nonsingular_solutions","member_points",CopyB'File=>true));
    readFile(NCO#"Directory");
    print(pairJac);
    print (peek first pairGradient);
    print(pairScale);
    print(peek first pairData);
    print(peek first pairWeight);
--hyperplane
    if stage==1 or stage==0 or stage==2 then hyperplaneAtInfinity:=makeB'Section(
    	sum apply(nc,i->{data_i*(primalList_i)}),
	B'NumberCoefficients=>{1},
	NameB'Section=>numerHB);
----Input file for stage 1 and stage 2.
    win:=L|G|randomizedCritConditions;
--Gradient conditions (RHS)  (STAGE ONE)
    if stage==0 or stage==1 then    RHS:=apply(nc,i->makeB'Section(
	{topDenomQ*data_i,-topNumerHB*primalList_i*startWeight_i},
	NameB'Section=>"RHS"|i,
	B'NumberCoefficients=>{topL0*topNumerHB^(first degRescale),toString(topL0*topNumerHB^(first degRescale))}
	));
    if stage==0 or stage==1 then isoQ:=makeB'Section(
    	sum apply(nc,i->{startWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ);
--Jacobian ring
    theR:=ring first F;
    numX:=#gens theR;
    kk1:=coefficientRing ring first F;
    jacS:=theR**kk1[apply(#L+#G,i->lagMult|i+1)]**kk1[{numerHB,denomQ}|{tWeight}];
    jacLamList:=flatten entries basis({0,1,0},jacS);
    (jacNumerHB,jacDenomQ,jacTWeight):=toSequence flatten entries basis({0,0,1},jacS);
    jacLV:=apply(jacLamList,drop(degRescale,1),(lam,j)->if j==0 
	then lam else if j>0 then lam*jacNumerHB^j 
	else print "Error: Homogenized incorrectly.");
    --print LV; 
    print 7;
    startJacStartCondition:=flatten entries (matrix{jacLV}*sub((jacL||jacG),jacS));    
    LHS:=apply(nc,i->"LHS"|i=>startJacStartCondition_i);
--Jacobian definition conditions (LHS) 
---SUPPORT FUNCTIONS
--Writing functions
    print 8;
    runWrite:=(stage,PG,nif,bfs)->(
	makeB'InputFile(NCO#"Directory",NameB'InputFile=>nif,
	    HomVariableGroup=>(NCO#"HomogeneousVariableGroups")|{{topL0}|jacLamList},
    	    AffVariableGroup=>NCO#"AffineVariableGroups",
	    B'Configs=>{
	    	"UseRegeneration"=>1,
	    	"TrackType"=>0,
	    	"ParameterHomotopy"=>stage,
	    	"PrintPathProgress"=>1000}|(NCO#"BertiniStartFiberSolveConfiguration"),
	    B'Polynomials=>win,
    	    ParameterGroup=>PG,
	    B'Functions=>bfs,
	    B'Constants=>NCO#"BertiniConstants"));    
    writeMTFile:=(s,k,bp,bfs)->makeB'InputFile(theDir,NameB'InputFile=>(defaultMTName|s|toString k),
	--Which variables come first?
	AffVariableGroup=>flatten flatten{NCO#"HomogeneousVariableGroups",{{topL0}|jacLamList},NCO#"AffineVariableGroups"},
	B'Configs=>{"TrackType"=>k,"PrintPathProgress"=>1000}|(NCO#"BertiniMembershipTestConfiguration"),
	B'Polynomials=>bp,
	B'Functions=>bfs,
	B'Constants=>NCO#"BertiniConstants"
	);
    print 9;
    writeManyMT:=(stage,bfs)->(
	writeMTFile("Residual"|stage|"_",1,{first last critConditions},bfs);
    	writeMTFile("Residual"|stage|"_",3,{first last critConditions},bfs);
    	writeMTFile("Degenerate"|stage|"_"|0|"_",1,{lagMult|"0"},bfs);
  	writeMTFile("Degenerate"|stage|"_"|1|"_",1,{numerHB},bfs);
  	writeMTFile("Degenerate"|stage|"_"|2|"_",1,{denomQ},bfs);
    	writeMTFile("Degenerate"|stage|"_"|0|"_",3,{lagMult|"0"},bfs);
  	writeMTFile("Degenerate"|stage|"_"|1|"_",3,{numerHB},bfs);
  	writeMTFile("Degenerate"|stage|"_"|2|"_",3,{denomQ},bfs);
    	scan(#F,i->(
    	    writeMTFile("Hypersurface"|stage|"_"|i|"_",1,{F_i},bfs);
    	    writeMTFile("Hypersurface"|stage|"_"|i|"_",3,{F_i},bfs);        ))  );
---Run membershipTest
    print 10;
    runMT:=(s)->(
   	runBertini(theDir,NameB'InputFile=>defaultMTName|s|"1");
	print (defaultMTName|s|"1");
    	moveB'File(theDir,"bertini_session.log","bertini_session_MT_"|s|"1.log",CopyB'File => true);
    	runBertini(theDir,NameB'InputFile=>defaultMTName|s|"3");
	print (defaultMTName|s|"3");
    	moveB'File(theDir,"bertini_session.log","bertini_session_MT_"|s|"3.log",CopyB'File => true);
    	outIM:=importIncidenceMatrix(theDir);
    	print outIM;
    	return outIM);                                                                           
    --Run Sort
    print 11;
    runSort:=(stage)->(
--	moveB'File(theDir,"nonsingular_solutions","member_points");
	print "isResidual";
	imResidual:=runMT("Residual"|stage|"_");
	--Degenerate
	print "isDegenerate";
	imAllDegenerate:=apply(3,i->runMT("Degenerate"|stage|"_"|i|"_"));
	--Hypersurfaces
	print "isMemberOfVariety";
	imAllHypersurfaces:=apply(#F,i->runMT("Hypersurface"|stage|"_"|i|"_"));
	EDDeg:=0;
	isMemberEveryHypersurface:=(i)->(
	    output:=true;
	    scan(imAllHypersurfaces,j->if {}==j_i then (output=false;break));
	    return output);
	isMemberAnyDegenerate:=(i)->(
	    output:=false;
	    scan(imAllDegenerate,j->if {}=!=j_i then (output=true;break));
	    return output);
    	ts:={};
	scan(#imResidual,i->print( imResidual_i=!={}, not isMemberAnyDegenerate(i), isMemberEveryHypersurface(i)));
	scan(#imResidual,i->if imResidual_i=!={} and not isMemberAnyDegenerate(i) and isMemberEveryHypersurface(i) 
	    then (EDDeg=EDDeg+1; ts=ts|{1}) else ts=ts|{0});
	NCO#"TrackSolutions"=ts;
	print("EDDeg",EDDeg,"Stage",stage);
	return EDDeg);
    --RUN  more writing functions
    print 11;    
--END Support functions
    print 12;
--STAGE ONE
    if stage==1 or stage==0
    then (
	nif:="inputCriticalPointSuperSet";
    	bfs:=NCO#"BertiniSubstitute"|primalSub|{hyperplaneAtInfinity,isoQ}|RHS|LHS|critConditions;
    	print(13,"A");
    	runWrite(stageOne,{"asdfadagds"},nif,bfs);
    	runBertini(theDir,NameB'InputFile=>nif);
	print nif; 
	moveB'File(theDir,"bertini_session.log","bertini_session_"|nif|".log",CopyB'File => true);
	moveB'File(theDir,"nonsingular_solutions","member_points",CopyB'File => true);
	moveB'File(theDir,"nonsingular_solutions","dummy");
    	print(13,"B");
	writeManyMT(stageOne,bfs);
    	print(13,"C");
    	GED:=runSort(stageOne));
    print 9;
    filterSolutionFile(NCO,"start",nc);
    print 10;
--Gradient conditions (RHS)  (STAGE TWO)
    if stage==2  or stage==0
    then (writeParameterFile(theDir,{0},NameParameterFile=>"start_parameters");
    	writeParameterFile(theDir,{1},NameParameterFile=>"final_parameters"));
    print("D",1);
    if stage==0 or stage==2 then    RHS=apply(nc,i->makeB'Section(
	{topDenomQ*data_i,-topNumerHB*primalList_i*startWeight_i,-topNumerHB*primalList_i*targetWeight_i},
	NameB'Section=>"RHS"|i,
	B'NumberCoefficients=>{topL0*topNumerHB^(first degRescale),
		    "(1-"|tWeight|")*"|toString(topL0*topNumerHB^(first degRescale)),
		    tWeight|"*"|toString(topL0*topNumerHB^(first degRescale))}
	));
    if stage==2 then print("RHSMT",peek first RHS);
    if stage==0 or stage==2 then    RHSMT:=apply(nc,i->makeB'Section(
	{topDenomQ*data_i,-topNumerHB*primalList_i*targetWeight_i},
	NameB'Section=>"RHS"|i,
	B'NumberCoefficients=>{topL0*topNumerHB^(first degRescale),		   
		    toString(topL0*topNumerHB^(first degRescale))}
	));
    if stage==2 then     print("RHSMT",peek first RHSMT);
    if stage==0 or stage==2 then isoQ=makeB'Section(
    	sum apply(nc,i->{startWeight_i*(primalList_i)^2,targetWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{"(1-"|tWeight|")",tWeight},
	NameB'Section=>denomQ);
    if stage==2 then     print ("ISOQ",peek isoQ);
    if stage==0 or stage==2 then isoQMT:=makeB'Section(
    	sum apply(nc,i->{targetWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ);
    if stage==2 then     print ("ISOQMT",peek isoQMT);
    if stage==2 or stage==0
    then (
	nif="inputParameterHomotopySuperSet";
    	bfs=NCO#"BertiniSubstitute"|primalSub|{hyperplaneAtInfinity,isoQ}|RHS|LHS|critConditions;
    	print(8,"B");
    	runWrite(stageTwo,{tWeight},nif,bfs);
    	runBertini(theDir,NameB'InputFile=>nif);
	print nif; 
	moveB'File(theDir,"bertini_session.log","bertini_session_"|nif|".log",CopyB'File => true);
	moveB'File(theDir,"nonsingular_solutions","member_points",CopyB'File => true);
	moveB'File(theDir,"nonsingular_solutions","dummy2");
    	print(8,"B");
	RHS=RHSMT;
	isoQ=isoQMT;
    	bfs2:=NCO#"BertiniSubstitute"|primalSub|{hyperplaneAtInfinity,isoQ}|RHS|LHS|critConditions;
	writeManyMT(stageTwo,bfs2);
    	print(8,"C");
    	UED:=runSort(stageTwo));
    return(GED=>UED));
*-
-*   
runBertiniStartEDDegree=method()
runBertiniStartEDDegree(NumericalComputationOptions):=(NCO)->runBertiniStartEDDegree(NCO,(0,0,0),1,#NCO#"Model")
runBertiniStartEDDegree(NumericalComputationOptions,Sequence,ZZ):=(NCO,ht,stage)->runBertiniStartEDDegree(NCO,ht,stage,#NCO#"Model")
runBertiniStartEDDegree(NumericalComputationOptions,Sequence,ZZ,ZZ):=(NCO,ht,stage,n)->(
    theDir:=NCO#"Directory";
    if stage===1 then nif:="inputCriticalPointSuperSet" 
    else if stage===2 then nif="inputParameterHomotopySuperSet";        
    (pvM,pvD,pvW):=ht;--weight is last
    pgStart:=ht/(i->if i===null then 0 else i)//toList;
    pgTarget:=ht/(i->if i===null then 1 else i)//toList;
    if stage==2 
    then (writeParameterFile(theDir,pgStart,NameParameterFile=>"start_parameters");
    	writeParameterFile(theDir,pgTarget,NameParameterFile=>"final_parameters"));
    runBertini(theDir,NameB'InputFile=>nif);
    print nif; 
    moveB'File(theDir,"bertini_session.log","bertini_session_"|nif|".log",CopyB'File => true);
    moveB'File(theDir,"nonsingular_solutions","member_points");
--Residual
    print "isResidual";
    imResidual:=runMT("Residual"|stage|"_");
--Degenerate
    print "isDegenerate";
    imDegenerate:=runMT("Degenerate"|stage|"_");
--Hypersurfaces
    print "isMemberOfVariety";
    imAllHypersurfaces:=apply(n,i->runMT("Hypersurface"|stage|"_"|i|"_"));
    EDDeg:=0;
    memberEveryHypersurface:=(i)->(
	output:=true;
	scan(imAllHypersurfaces,j->if {}==j_i then (output=false;break));
	return output);
    ts:={};
    scan(#imResidual,i->if imResidual_i=!={} and imDegenerate_i=={} and memberEveryHypersurface(i) 
	then (EDDeg=EDDeg+1; ts=ts|{1}) else ts=ts|{0});
    NCO#"TrackSolutions"=ts;
--    moveB'File(theDir,"bertini_session_CriticalPointSuperSet.log","bertini_session.log");
--    return(imResidual,imDegenerate)
    return(EDDeg)
     )
*- 













-*
weightEDDegreeHomotopy=method()
weightEDDegreeHomotopy(String,Sequence,List):=(theDir,P,TWV)->(
    print ("A",0);
    NCO:=newNumericalComputationOptions(theDir,P);
    NCO#"TargetWeight"=TWV;
    nc:=#gens ring first (NCO#"Model")+1+#NCO#"WitnessModel"+#NCO#"StartSubmodel";
    print("nc",nc);
    numF:=#(NCO#"Model");
    homotopyType:=(0,0,null);
    print ("A",1);
    startEDDegree(NCO,homotopyType,stageOne);
    print ("A",2);
    GED:=runBertiniStartEDDegree(NCO,homotopyType,stageOne);    
    print (GED,"GED");
    filterSolutionFile(NCO,"start",nc) ;
    readFile(NCO#"Directory","start",100000);
    print ("A",3);
    startEDDegree(NCO,homotopyType,stageTwo);
    UED:=runBertiniStartEDDegree(NCO,homotopyType,stageTwo);
    print ("A",4);
    return  (GED=>UED)
    )
*-    
end

 
---EXAMPLES
restart
--path=prepend("/Users/jo/Documents/GoodGit/EuclideanDistanceDegree",path)
--loadPackage("EuclideanDistanceDegree",Reload=>true)

--[EX 1] 
R=QQ[x0,x1,x2,y0,y1,y2];
F=flatten entries gens minors(2,genericMatrix(R,3,2))
G=drop(F,-(#F-2))
L={}
startEDDegree(storeBM2Files,F,G,L,Weight=>"Generic")
runBertiniStartEDDegree(storeBM2Files,#F)--10

startEDDegree(storeBM2Files,F,G,{},Weight=>"Unit")
runBertiniStartEDDegree(storeBM2Files,#F)--2

--[EX 2] 
----Key words: 3 by 3 matrices rank two matrices, singular variety
--storeBM2Files="/Users/jo/Desktop/BertiniOutputFiles/****"
R=QQ[x11,x12,x13,
    x21,x22,x23,
    x31,x32,x33];
F=flatten entries gens minors(3,genericMatrix(R,3,3))
G=F
L={}

startEDDegree(storeBM2Files,F,G,{},Weight=>"Unit")
runBertiniStartEDDegree(storeBM2Files,#F)--3

startEDDegree(storeBM2Files,F,G,L,Weight=>"Generic")
runBertiniStartEDDegree(storeBM2Files,#F)--39

--These results agree with those on page 5 of https://www.researchgate.net/publication/258374161_Exact_Solutions_in_Structured_Low-Rank_Approximation/download
L={x12-x21,x31-x22,x31-x13,x32-x23}

startEDDegree(storeBM2Files,F,G,L,Weight=>"Unit")
runBertiniStartEDDegree(storeBM2Files,#F)--9
startEDDegree(storeBM2Files,F,G,L,Weight=>"Generic")
runBertiniStartEDDegree(storeBM2Files,#F)--13


---
xVars=flatten apply(5,i->apply(5,j->value("x"|i|j)))
R=CC[xVars];
M=genericMatrix(R,5,5)
F=flatten entries gens minors(3,M)
--G=apply(5-2,i->apply(5-2,j->({0,1,2}+{i,i,i},{0,1,2}+{j,j,j})))
G=flatten apply(5-2,i->apply(5-2,j->det submatrix(M,{0,1,2}+{i,i,i},{0,1,2}+{j,j,j})))
--L={random({1},R)}
L={}
startEDDegree(storeBM2Files,F,G,{},Weight=>"Unit")
runBertiniStartEDDegree(storeBM2Files,#F)--5 choose 2

startEDDegree(storeBM2Files,F,G,L,Weight=>"Generic")
runBertiniStartEDDegree(storeBM2Files,#F)--


R=QQ[t]






---
--[EX 1] 
R=QQ[x,y];
F={x^4+y^4-1}
G=F
L={}
apropos "EDD"
determinantalUnitEuclideanDistanceDegree(F)
determinantalGenericEuclideanDistanceDegree(F)
leftKernelUnitEDDegree(storeBM2Files,1,F)
runBertiniEDDegree(storeBM2Files)

symbolicWeightEDDegree
startEDDegree(storeBM2Files,F,G,L,Weight=>"Generic")
runBertiniStartEDDegree(storeBM2Files,#F)--10

stageWeightEDDegreeHomotopy=method()
stageWeightEDDegreeHomotopy(NumericalComputationOptions,ZZ):=(NCO,stage)->(    
--stage is either 0, 1 or 2.
--Now we extract information from NCO.
    theDir:=NCO#"Directory"; 
    data:=NCO#(NCO#"FixedData");
    jacL:=NCO#(NCO#"FixedJacobianSubmodel");
    L:=NCO#(NCO#"FixedSubmodel");
--Other strings that are used
    (lagMult,numerHB,denomQ,primal,tWeight):=("lagMult","numerHB","denomQ","primal","tWeight")   ; 
    startWeight:=NCO#"StartWeight";
    targetWeight:=NCO#"TargetWeight";
--First the model and submodel constraints
------F is the model, V(G)\cap V(L) is a complete intersection contained in V(F)\cap V(L).
    (F,G,startL,targetL,jacG):=(NCO#"Model",NCO#"WitnessModel",NCO#"StartSubmodel",NCO#"TargetSubmodel",NCO#"JacobianWitnessModel");
    randomGamma:=NCO#"GammaVector";
    xList:=NCO#"PrimalCoordinates";
    nc:=#xList;
-- --Homogenize appropriately
    (degSubmodel,degWitnessModel):=(NCO#"DegreeSubmodel",NCO#"DegreeWitnessModel");
    maxDegree:=(degSubmodel|degWitnessModel|{3})//max;
    degRescale:=({3}|degSubmodel|degWitnessModel)/(i->maxDegree-i);        
    print 1;
--We have three rings. One ring for manipulating the Jacobian, one ring for manipulating the gradient, and one for manipulating both
--Gradient ring
    print 2;
    kk2:=ring first startWeight;
    topS:=kk2[apply(nc,i->primal|i)]**kk2[{lagMult|"0",numerHB,denomQ,tWeight}];
    (topL0,topNumerHB,topDenomQ,topTWeight):=toSequence flatten entries basis({0,1},topS);
    primalList:=flatten entries basis({1,0},topS);
    primalSub:=transpose{primalList,xList};
---Critical ring (LHS-RHS)
    kk3:=ring first randomGamma    ;
    critRing:=kk3[apply(nc,i->"critCondition"|i)];
    critConditions:=apply(nc,i->(gens critRing)_i=>"LHS"|i|"-RHS"|i);
    print critConditions;
    randomizedCritConditions:=apply(drop(gens critRing,-1),randomGamma,(i,j)->i+j*last gens critRing);
--hyperplane
    if stage==1 or stage==0 or stage==2 then hyperplaneAtInfinity:=makeB'Section(
    	sum apply(nc,i->{data_i*(primalList_i)}),
	B'NumberCoefficients=>{1},
	NameB'Section=>numerHB);
----Input file for stage 1 and stage 2.
    win:=L|G|randomizedCritConditions;
--Gradient conditions (RHS)  (STAGE ONE)
    if stage==0 or stage==1 then    RHS:=apply(nc,i->makeB'Section(
	{topDenomQ*data_i,-topNumerHB*primalList_i*startWeight_i},
	NameB'Section=>"RHS"|i,
	B'NumberCoefficients=>{topL0*topNumerHB^(first degRescale),toString(topL0*topNumerHB^(first degRescale))}
	));
    if stage==0 or stage==1 then isoQ:=makeB'Section(
    	sum apply(nc,i->{startWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ);
--Jacobian ring
    theR:=ring first F;
    numX:=#gens theR;
    kk1:=coefficientRing ring first F;
    jacS:=theR**kk1[apply(#L+#G,i->lagMult|i+1)]**kk1[{numerHB,denomQ}|{tWeight}];
    jacLamList:=flatten entries basis({0,1,0},jacS);
    (jacNumerHB,jacDenomQ,jacTWeight):=toSequence flatten entries basis({0,0,1},jacS);
    jacLV:=apply(jacLamList,drop(degRescale,1),(lam,j)->if j==0 
	then lam else if j>0 then lam*jacNumerHB^j 
	else print "Error: Homogenized incorrectly.");
    --print LV; 
    print 7;
    startJacStartCondition:=flatten entries (matrix{jacLV}*sub((jacL||jacG),jacS));    
    LHS:=apply(nc,i->"LHS"|i=>startJacStartCondition_i);
--Jacobian definition conditions (LHS) 
---SUPPORT FUNCTIONS
--Writing functions
    print 8;
    runWrite:=(stage,PG,nif,bfs)->(
	makeB'InputFile(theDir,NameB'InputFile=>nif,
	    HomVariableGroup=>(NCO#"HomogeneousVariableGroups")|{{topL0}|jacLamList},
    	    AffVariableGroup=>NCO#"AffineVariableGroups",
	    B'Configs=>{
	    	"UseRegeneration"=>1,
	    	"TrackType"=>0,
	    	"ParameterHomotopy"=>stage,
	    	"PrintPathProgress"=>1000}|(NCO#"BertiniStartFiberSolveConfiguration"),
	    B'Polynomials=>win,
    	    ParameterGroup=>PG,
	    B'Functions=>bfs,
	    B'Constants=>NCO#"BertiniConstants"));    
    writeMTFile:=(s,k,bp,bfs)->makeB'InputFile(theDir,NameB'InputFile=>(defaultMTName|s|toString k),
	--Which variables come first?
	AffVariableGroup=>flatten flatten{NCO#"HomogeneousVariableGroups",{{topL0}|jacLamList},NCO#"AffineVariableGroups"},
	B'Configs=>{"TrackType"=>k,"PrintPathProgress"=>1000}|(NCO#"BertiniMembershipTestConfiguration"),
	B'Polynomials=>bp,
	B'Functions=>bfs,
	B'Constants=>NCO#"BertiniConstants"
	);
    print 9;
    writeManyMT:=(stage,bfs)->(
	writeMTFile("Residual"|stage|"_",1,{first last critConditions},bfs);
    	writeMTFile("Residual"|stage|"_",3,{first last critConditions},bfs);
    	writeMTFile("Degenerate"|stage|"_"|0|"_",1,{lagMult|"0"},bfs);
  	writeMTFile("Degenerate"|stage|"_"|1|"_",1,{numerHB},bfs);
  	writeMTFile("Degenerate"|stage|"_"|2|"_",1,{denomQ},bfs);
    	writeMTFile("Degenerate"|stage|"_"|0|"_",3,{lagMult|"0"},bfs);
  	writeMTFile("Degenerate"|stage|"_"|1|"_",3,{numerHB},bfs);
  	writeMTFile("Degenerate"|stage|"_"|2|"_",3,{denomQ},bfs);
    	scan(#F,i->(
    	    writeMTFile("Hypersurface"|stage|"_"|i|"_",1,{F_i},bfs);
    	    writeMTFile("Hypersurface"|stage|"_"|i|"_",3,{F_i},bfs);        ))  );
---Run membershipTest
    print 10;
    runMT:=(s)->(
   	runBertini(theDir,NameB'InputFile=>defaultMTName|s|"1");
	print (defaultMTName|s|"1");
    	moveB'File(theDir,"bertini_session.log","bertini_session_MT_"|s|"1.log",CopyB'File => true);
    	runBertini(theDir,NameB'InputFile=>defaultMTName|s|"3");
	print (defaultMTName|s|"3");
    	moveB'File(theDir,"bertini_session.log","bertini_session_MT_"|s|"3.log",CopyB'File => true);
    	outIM:=importIncidenceMatrix(theDir);
    	print outIM;
    	return outIM);                                                                           
    --Run Sort
    print 11;
    runSort:=(stage)->(
--	moveB'File(theDir,"nonsingular_solutions","member_points");
	print "isResidual";
	imResidual:=runMT("Residual"|stage|"_");
	--Degenerate
	print "isDegenerate";
	imAllDegenerate:=apply(3,i->runMT("Degenerate"|stage|"_"|i|"_"));
	--Hypersurfaces
	print "isMemberOfVariety";
	imAllHypersurfaces:=apply(#F,i->runMT("Hypersurface"|stage|"_"|i|"_"));
	EDDeg:=0;
	isMemberEveryHypersurface:=(i)->(
	    output:=true;
	    scan(imAllHypersurfaces,j->if {}==j_i then (output=false;break));
	    return output);
	isMemberAnyDegenerate:=(i)->(
	    output:=false;
	    scan(imAllDegenerate,j->if {}=!=j_i then (output=true;break));
	    return output);
    	ts:={};
	scan(#imResidual,i->print( imResidual_i=!={}, not isMemberAnyDegenerate(i), isMemberEveryHypersurface(i)));
	scan(#imResidual,i->if imResidual_i=!={} and not isMemberAnyDegenerate(i) and isMemberEveryHypersurface(i) 
	    then (EDDeg=EDDeg+1; ts=ts|{1}) else ts=ts|{0});
	NCO#"TrackSolutions"=ts;
	print("EDDeg",EDDeg,"Stage",stage);
	return EDDeg);
    --RUN  more writing functions
    print 11;    
--END Support functions
    print 12;
--STAGE ONE
    if stage==1 or stage==0
    then (
	nif:="inputCriticalPointSuperSet";
    	bfs:=NCO#"BertiniSubstitute"|primalSub|{hyperplaneAtInfinity,isoQ}|RHS|LHS|critConditions;
    	print(13,"A");
    	runWrite(stageOne,{"asdfadagds"},nif,bfs);
    	runBertini(theDir,NameB'InputFile=>nif);
	print nif; 
	moveB'File(theDir,"bertini_session.log","bertini_session_"|nif|".log",CopyB'File => true);
	moveB'File(theDir,"nonsingular_solutions","member_points",CopyB'File => true);
	moveB'File(theDir,"nonsingular_solutions","dummy");
    	print(13,"B");
	writeManyMT(stageOne,bfs);
    	print(13,"C");
    	GED:=runSort(stageOne));
    print 9;
    filterSolutionFile(NCO,"start",nc);
    print 10;
--Gradient conditions (RHS)  (STAGE TWO)
    if stage==2  or stage==0
    then (writeParameterFile(theDir,{0},NameParameterFile=>"start_parameters");
    	writeParameterFile(theDir,{1},NameParameterFile=>"final_parameters"));
    print("D",1);
    if stage==0 or stage==2 then    RHS=apply(nc,i->makeB'Section(
	{topDenomQ*data_i,-topNumerHB*primalList_i*startWeight_i,-topNumerHB*primalList_i*targetWeight_i},
	NameB'Section=>"RHS"|i,
	B'NumberCoefficients=>{topL0*topNumerHB^(first degRescale),
		    "(1-"|tWeight|")*"|toString(topL0*topNumerHB^(first degRescale)),
		    tWeight|"*"|toString(topL0*topNumerHB^(first degRescale))}
	));
    if stage==2 then print("RHSMT",peek first RHS);
    if stage==0 or stage==2 then    RHSMT:=apply(nc,i->makeB'Section(
	{topDenomQ*data_i,-topNumerHB*primalList_i*targetWeight_i},
	NameB'Section=>"RHS"|i,
	B'NumberCoefficients=>{topL0*topNumerHB^(first degRescale),		   
		    toString(topL0*topNumerHB^(first degRescale))}
	));
    if stage==2 then     print("RHSMT",peek first RHSMT);
    if stage==0 or stage==2 then isoQ=makeB'Section(
    	sum apply(nc,i->{startWeight_i*(primalList_i)^2,targetWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{"(1-"|tWeight|")",tWeight},
	NameB'Section=>denomQ);
    if stage==2 then     print ("ISOQ",peek isoQ);
    if stage==0 or stage==2 then isoQMT:=makeB'Section(
    	sum apply(nc,i->{targetWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ);
    if stage==2 then     print ("ISOQMT",peek isoQMT);
    if stage==2 or stage==0
    then (
	nif="inputParameterHomotopySuperSet";
    	bfs=NCO#"BertiniSubstitute"|primalSub|{hyperplaneAtInfinity,isoQ}|RHS|LHS|critConditions;
    	print(8,"B");
    	runWrite(stageTwo,{tWeight},nif,bfs);
    	runBertini(theDir,NameB'InputFile=>nif);
	print nif; 
	moveB'File(theDir,"bertini_session.log","bertini_session_"|nif|".log",CopyB'File => true);
	moveB'File(theDir,"nonsingular_solutions","member_points",CopyB'File => true);
	moveB'File(theDir,"nonsingular_solutions","dummy2");
    	print(8,"B");
	RHS=RHSMT;
	isoQ=isoQMT;
    	bfs2:=NCO#"BertiniSubstitute"|primalSub|{hyperplaneAtInfinity,isoQ}|RHS|LHS|critConditions;
	writeManyMT(stageTwo,bfs2);
    	print(8,"C");
    	UED:=runSort(stageTwo));
    return(GED=>UED));
