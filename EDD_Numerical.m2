--restart
--Projective formulation for intersections with linear spaces
rand:=randomValue
--Assume ring is a complex inexact field
--G is a subset of F. 
NumericalComputationOptions=new Type of MutableHashTable
(stageOne,stageTwo)=(1,2); 
  
    
    
-*
installPackage"EuclideanDistanceDegree"
check EuclideanDistanceDegree

restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
check EuclideanDistanceDegree
help EuclideanDistanceDegree
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
    NCO#"IsProjective"=false;
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
    print("homogenized jacobian",homogJac);
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
    	print("Membership tests",nif);	
	print outIM;
	return outIM	)  ;  
-- (FUNCTIONS 4) Functions for filtering based off of incidence matrix
    filterSolutionFunction:=(nsf,kp,ns,numCoords)->(     
--    	print("RUN FILTER",kp=>numCoords);    
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
--	isMembershipFunction(stage,case,indexCase,hypersurface);      
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
	(nsf,nc):=("filterFile",#flatten {bLagrangeVars,bModelVars});
    	print("Filter",kp,"num kp",#kp,"num sols",ns,"num coordinates",nc,bin);
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
--	print("onPolyList",onPolyList);
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


--restart
--loadPackage"Bertini"
-*
parameterizeConormalMatrixRankConstraint=(r,m,n,L,kk)->(
    matrixVars:=(a,b,s)->for i to a-1 list for j to b-1 list (s|i|"v"|j);    
    getMatrix=(aRing,a,b)->transpose genericMatrix(aRing,b,a);
    getId=(m,v)->diagonalMatrix apply(m,i->v);
--
    ringPrimalSpan:=kk[(m-r,r,"PS")//matrixVars//flatten|{"PSH"}];
    ringPrimalCenter:=kk[(r,r,"PC")//matrixVars//flatten|{"PCH"}];
    ringDualSpan:=kk[(n-r,r,"DS")//matrixVars//flatten|{"DSH"}];
    ringDualCenter:=kk[(m-r,n-r,"DC")//matrixVars//flatten|{"DCH"}];
    ringData:=kk[flatten matrixVars(m,n,"u")];
    ringWeight:=kk[flatten matrixVars(m,n,"w")];
    allRings:={ringPrimalSpan,ringPrimalCenter,ringDualSpan,ringDualCenter,ringData,ringWeight};
    sizeMatrix={(m-r,r),(r,r),(n-r,r), (m-r,n-r), (m,n),  (m,n)}; 
    allMatrices= apply(allRings,sizeMatrix,(i,j)->getMatrix(i,j_0,j_1));
    print allMatrices; 
    lagRing=kk[apply(#L+1,i->"lag"|i)];
    bigRing=ringData**ringWeight**lagRing;
    scan(reverse drop(allRings,-2),i->bigRing=i**bigRing);
    allHomogenizers= apply(allRings,i->last gens i);
    print allHomogenizers;
    (hps,hpc,hds,hdc,hu,hw)= (allHomogenizers/(i->sub(i,bigRing)))//toSequence;
    print(hps,hpc,hds,hdc,hu,hw);
    U=sub(transpose genericMatrix(ringData,n,m),bigRing);
    print U;
    W=transpose sub(genericMatrix(ringWeight,n,m),bigRing);
    print W;
    print 2;
    mp0=sub(getId(last sizeMatrix_0,allHomogenizers_0)||allMatrices_0,bigRing);
    print 3;
    mp1=sub(allMatrices_1,bigRing);
    print W;
    print 4;
    mp2=sub(getId(last sizeMatrix_0,allHomogenizers_2)|transpose allMatrices_2,bigRing);
    print (mp0,mp1,mp2);
    md0=sub((transpose allMatrices_0)||getId(first sizeMatrix_3,-allHomogenizers_0),bigRing);
    md1=sub(allMatrices_3,bigRing);
    print md1;    
    md2=transpose sub(transpose allMatrices_2||getId(last sizeMatrix_3,-allHomogenizers_2),bigRing);
    print md2;
    print (md0,md1,md2);
    P=hdc*(product\\{mp0,mp1,mp2});
    print P;
    Q=hpc*(product\\{md0,md1,md2});
    print Q;
    scan(#L,i->Q=Q+hpc*hps*hds*sub((gens lagRing)_(i+1),bigRing )*sub(L_i,bigRing));
    linearSpace:=apply(#L,i->makeB'Section(flatten entries P));    
    print linearSpace;
    dataSlices={};
--    (hps,hpc,hds,hdc,hu,hw)
    scan(m,i->scan(n,
	    j->dataSlices=dataSlices|{
		makeB'Section({P_(i,j),-Q_(i,j),hpc*hdc*hps*hds*U_(i,j)},
		    B'NumberCoefficients=>{1,W_(i,j),-1})})
    ); --(x-u=w y) (x-u)
--    return((mp0,mp1,mp2),(md0,md1,md2),dataSlices,linearSpace)
    dir=temporaryFileName();
    print dir;
    if not fileExists dir then mkdir dir;
--    HVG={gens ringPrimalSpan, gens ringPrimalCenter, 
--	gens ringDualSpan, (gens ringDualCenter)|drop(gens lagRing,1)};    
    HVG=flatten{drop(gens ringPrimalSpan,-1), drop(gens ringPrimalCenter,-1), 
	drop(gens ringDualSpan,-1), (gens ringDualCenter)|drop(gens lagRing,1)};    
    makeB'InputFile(dir,
    	HomVariableGroup=>HVG,
    	B'Polynomials=>linearSpace|dataSlices,
    	B'Functions=>{PSH=>DCH,PCH=>DCH,DSH=>DCH},
    	B'Configs=>{"UseRegeneration"=>1},
    	RandomComplex=>(flatten entries U),
	B'Constants=>apply(flatten entries W,i->i=>1 )	
	);
    return dir
    )
-*

-*
rm=()->matrix for i to 5-1 list for j to 5-1 list random CC
rm()
dir =parameterizeConormalMatrixRankConstraint(2,5,5,{rm(),rm(),rm()},CC)        
runBertini(dir)
R=QQ[p1,p2,p3,p4]

pm=genericMatrix(R,2,2)
l=p1-p2
diff(pm,l)
*-

-*
GEd(X cap L^s)-UED(X cap L^2)=GED(Z cap L^s)
where Z=sing(X cap isoQ)
R=QQ[x1,x2,x3,x4]
Q=ideal sum (gens R/(i->i^2))
X=minors(2,genericMatrix(R,2,2))
decompose ideal singularLocus(Q+X)--4 points

loadPackage"Bertini"
R=QQ[x1,x2,x3,x4,x5,x6,x7,x8,x9]
Q=ideal sum (gens R/(i->i^2))
X=minors(2,genericMatrix(R,3,3))
X=minors(2,genericMatrix(R,2,4))--2 components of degree 2 in dim3
--bertiniPosDimSolve flatten entries gens ideal singularLocus(Q+X)
Z=flatten entries gens radical ideal singularLocus(Q+X)
bertiniPosDimSolve Z--This has two components of degree 2. 
printGens  
ideal Z==radical ideal Z
conjecture--
loadPackage"Bertini"
(r,m,n)=(1,2,4)

R=QQ[apply(m*n,i->"x"|i)]**QQ[apply((m-r)*(n-r),i->"y"|i)];
Q=basis({1,0})//entries//flatten/(i->i^2)//sum//ideal
X=genericMatrix(R,m,n)


bertiniPosDimSolve flatten entries gens ideal singularLocus(Q+X)


*-



-*
experimentDualityDifference=method()
dualVariety=(F,codF)->(
    R1:=ring first F;
    kk:=coefficientRing R1;
    R2:=kk[apply(#gens R1,i->"y"|i)];
    R:=R1**R2;
    M:=sub(matrix {gens R2},R)||sub(matrix makeJac(F,gens R1),R);
    win:=sub(ideal(F),R)+minors(codF+1,M);        
    E:=sub(eliminate(flatten entries basis({1,0},R),first decompose win),R2);
    print"computed dual variety";
    return E
    )
experimentDualityDifference(List,ZZ,ZZ):=(G,s,codF)->(
    R1:=ring first G;
    HS:=apply(s,i->random({1},R1));
    F:=G|HS;    
    adeg:=determinantalGenericEuclideanDistanceDegree F;
    print("generic X cap Hs",adeg);
    bdeg:=determinantalUnitEuclideanDistanceDegree F;
    print("unit X cap Hs",bdeg);
    R2:=QQ[apply(#gens R1,i->"y"|i)];
    dF:=flatten entries gens dualVariety(F,codF+s);
    cdeg:=determinantalGenericEuclideanDistanceDegree dF;
    print("generic X^* cap Hs",cdeg);
    ddeg:=determinantalUnitEuclideanDistanceDegree dF;
    print("unit X* cap Hs",ddeg);
    print(adeg,bdeg,cdeg,ddeg);
    adeg-bdeg==cdeg-ddeg )
*-


-*
restart
loadPackage"EuclideanDistanceDegree"
experimentDualityDifference
R=QQ[x0,x1,x2,x3,x4,x5,x6]
F={sum apply(gens R,i->i*i^2)}
IQ=(L)->sum  apply(L,i->i^2)
F={x0+x1,x2-x3,IQ({x1,x2,x3,x4,x5,x6})}
(s,codF)=(1,#F)
experimentDualityDifference(F,s,codF)

R=QQ[x0,x1,x2,x3,x4,x5,x6]
F={sum apply(gens R,i->i*i^2)}
IQ=(L)->sum  apply(L,i->i^2)
F={x0+x1,x2-x3,IQ {x1,x2,x3}, IQ{x2,x4,x5}}
(s,codF)=(2,#F)
experimentDualityDifference(F,s,codF)

*-



-*
restart
loadPackage"EuclideanDistanceDegree"
R=QQ[x,y,z,w]
--Calyx
F={x^2+y^2*z^3-z^4}/(i->homogenize(i,w))
experimentDualityDifference(F,1,1)

determinantalUnitEuclideanDistanceDegree(F)--21
primaryDecomposition ideal singularLocus ideal F

--Calypso
F={x^2+y^2*z-z^2}/(i->homogenize(i,w))
experimentDualityDifference(F,1,1)

determinantalUnitEuclideanDistanceDegree(F)--10
primaryDecomposition ideal singularLocus ideal F

--Crixxi
F={(y^2+z^2-1)^2 +(x^2+y^2-1)^3}/(i->homogenize(i,w))
experimentDualityDifference(F,1,1)

determinantalUnitEuclideanDistanceDegree(F)--22
primaryDecomposition ideal singularLocus ideal F

--Cube
R=QQ[x,y,z]
F={x^6+y^6+z^6-1}/(i->homogenize(i,w))
experimentDualityDifference(F,1,1)

determinantalUnitEuclideanDistanceDegree(F)--180
primaryDecomposition ideal singularLocus ideal F

--Geisha
F={x^2*y*z + x^2*z^2 - y^3*z - y^3 }/(i->homogenize(i,w))
experimentDualityDifference(F,1,1)

determinantalUnitEuclideanDistanceDegree(F)--15
primaryDecomposition ideal singularLocus ideal F

--Helix
F={6*x^2-2*x^4-y^2*z^2}/(i->homogenize(i,w))
experimentDualityDifference(F,1,1)
determinantalUnitEuclideanDistanceDegree(F)--20
primaryDecomposition ideal singularLocus ideal F

--Himmel und Hölle --RECUCIBLE
F={x^2-y^2*z^2}/(i->homogenize(i,w))
experimentDualityDifference(F,1,1)
determinantalUnitEuclideanDistanceDegree(F)--10
primaryDecomposition ideal singularLocus ideal F

--Kolibri
F={x^3 + x^2*z^2 - y^2}/(i->homogenize(i,w))
experimentDualityDifference(F,1,1)
determinantalUnitEuclideanDistanceDegree(F)--15
primaryDecomposition ideal singularLocus ideal F




*-



































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
--restart
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
P=(F,F,L)
NCO=newNumericalComputationOptions(storeBM2Files,P)

methods homotopyEDDegree
homotopyEDDegree(NCO,"Weight",true,true)

(storeBM2Files,F,G,{},Weight=>"Unit")
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

restart
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x0,x1,x2]
--f=x0^2*x1-x2*x3^2
f=x0*x1-x2*(x1-x0)
F={f}
Z=radical  (ideal(sum apply(gens R,i->i^2))+ideal F)
bertiniPosDimSolve flatten entries gens Z-- two lines
determinantalGenericEuclideanDistanceDegree F
determinantalUnitEuclideanDistanceDegree F


determinantalGenericEuclideanDistanceDegree Z
determinantalUnitEuclideanDistanceDegree Z




restart
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x0,x1,x2,t,s,lam]
xList={x0,x1,x2,t,s}
f=x0*x1-x2*(x1-x0)
--f=x0^2+x1^2-x2^2
--F={f,t*(x0^2+x1^2+x2^2)}

xList={x0,x1,x2}
F={f}
Z=radical  (ideal(sum apply(xList,i->i^2))+ideal F)
bertiniPosDimSolve flatten entries gens Z-- two lines
determinantalGenericEuclideanDistanceDegree F
determinantalUnitEuclideanDistanceDegree F
win=ideal F+ minors(#F+2,(
    matrix transpose apply(xList,v->{v*(s+t*random(1,1000)),random(1,1000)})  )||matrix makeJac(F,xList)
    )
decWin= decompose win;
netList decWin
decWin/degree
I= sub(first decWin,{t=>0,s=>1})
decompose  I
factor ( eliminate({x2},I))_0
coefficients( ( eliminate({x2},first decWin))_0,Variables=>xList)
decompose I




restart
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x0,x1,x2,x3,t,s]
f=x0*x1-x2*x3
xList={x0,x1,x2,x3}
F={f}
Z=radical  (ideal(sum apply(xList,i->i^2))+ideal F)
bertiniPosDimSolve flatten entries gens Z-- two lines
determinantalGenericEuclideanDistanceDegree F
determinantalUnitEuclideanDistanceDegree F
win=ideal F+ minors(#F+2,(
    matrix transpose apply(xList,v->{v*(s+t*random(1,1000)),random(1,1000)})  )||matrix makeJac(F,xList)
    )
decWin= decompose win;

netList decWin
decWin/degree
I= sub(first decWin,{t=>0,s=>1})
decI=decompose I
decI/(i->i==i+Z)
oo/degree



I= sub(first decWin,{})

coefficients( ( eliminate({x1,x2},I))_0,Variables=>xList)


decompose  I
oo/degree
coefficients( ( eliminate({x2},I))_0,Variables=>xList)

coefficients( ( eliminate({x2},first decWin))_0,Variables=>xList)
decompose I

---
restart
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x0,x1,x2,x3,y,t,s,lam]
xList={x0,x1,x2,x3,y}
f=x0*x1-x2*x3
g=(-y^2+sum apply(xList,i->i^2))-y
F={f,g}
Z=radical  (ideal(sum apply(xList,i->i^2))+ideal F)
bertiniPosDimSolve flatten entries gens Z-- two lines
determinantalGenericEuclideanDistanceDegree F
determinantalUnitEuclideanDistanceDegree F
win=ideal F+ minors(#F+2,(
    matrix transpose apply(xList,v->{v*(random(1,1000)),random(1,1000)})  )||matrix makeJac(F,xList)
    )

win=ideal F+ minors(#F+2,(
    matrix transpose apply(xList,v->{v*(s+t*random(1,1000)),random(1,1000)})  )||matrix makeJac(F,xList)
    )

printGens win
decWin= decompose win;

netList decWin
decWin/degree
I= sub(first decWin,{t=>0,s=>1})
decompose I

I= sub(first decWin,{})

coefficients( ( eliminate({x1,x2},I))_0,Variables=>xList)


restart
loadPackage"EuclideanDistanceDegree"
R=QQ[x1,x2,x3,x4]
f=det matrix{{x1,x2,x4},
    {x4,0,x3},
    {0,x3,x4}}
factor f
F={f}
determinantalGenericEuclideanDistanceDegree F
determinantalUnitEuclideanDistanceDegree F


restart
loadPackage"EuclideanDistanceDegree"
R=QQ[x1,x2,x3,x4,x5,x6,x7,x8,x9]
M=transpose genericMatrix(R,3,3)
F={x2-x4,x8-x6,x2-x8,det M}--(unit,generic)=(15,21)
--bertiniPosDimSolve F--codim 4 and degree is 3. 
#F==4
P=(F,F,{})
theDir
NCO=newNumericalComputationOptions(theDir,P)
help EuclideanDistanceDegree
homotopyEDDegree(NCO,"Weight",true,true)                    



restart
loadPackage"EuclideanDistanceDegree"
R=QQ[x1,x2,x3,x4,x5,x6,x7,x8,x9]
M=transpose genericMatrix(R,3,3)
F={x2-x4,x8-x6,det M}--(generic,unit)=(39,15)
bertiniPosDimSolve F--codim 4 and degree is 3. 
Q=ideal sum apply(gens R,i->i^2)
--sl= ideal F+ideal mingens ideal singularLocus(Q+ideal F)
--decSL=decompose sl;
 netList decSL;
decSL/codim---{7,7,7}
decSL/degree --{2, 2, 4}  -->ED degrees are (2,2,10)
---each of the first two lines intersects each of the second two lines (like a square)
--2 of these points are also on the quartic. 
--Intersection lattice: {line1a,line1b,line2a,line2b,quartic}--{2,2,2,2,2}
primaryDecomposition ideal mingens ideal singularLocus sl ---4points
radical (decSL_0+decSL_2)
degree radical sum decSL

#F==3
P=(F,F,{})
theDir
NCO=newNumericalComputationOptions(theDir,P)
help EuclideanDistanceDegree
 NCO#"StartData"=apply(#gens R,i->random RR)
 NCO#"StartData"=(.525992249534688+.888630754215339*ii)*{.474323, .269152, .859334, .370283, .785498, .263867, .468432, .743158, .514817}
homotopyEDDegree(NCO,"Weight",true,true)                    
cp=importSolutionsFile(theDir,NameSolutionsFile=>"criticalPointFile");
first cp

#cp
S01=decompose radical(decSL_0+decSL_1)
S02=decompose radical(decSL_0+decSL_2)
S12=decompose radical(decSL_1+decSL_2)
netList apply(S01|S02|S12,i->apply(S01|S02|S12,j->i==j))
S01--two different components
radical (sum S01)
decompose radical sum decSL


radical sum (S01|S02|S12)

--singular locus is four lines. 
quartic=sub(ideal(x3-x7,x1+x5+x9,x6-x8,x2-x4,x4^2+x5^2+x5*x9,x8^2-x5*x9,x4*x7+x5*x8+x8*x9,x5*x7-x4*x8,x7*x8-x4*x9,x7^2+x5*x9+x9^2),R)
singularLocus quartic
bertiniPosDimSolve flatten entries gens quartic
codim quartic
numgens quartic
printGens quartic
G={(x3-x7),
(x1+x5+x9),
(x6-x8),
(x2-x4),
(x4^2+x5^2+x5*x9),
(x8^2-x5*x9),
--(x4*x7+x5*x8+x8*x9),
--(x5*x7-x4*x8),
--(x7*x8-x4*x9),
(x7^2+x5*x9+x9^2)}/(i->sub(i,R))
#G
#G
codim quartic
determinantalUnitEuclideanDistanceDegree
bertiniPosDimSolve G
P2=(flatten entries gens quartic,G,{})
theDir2=storeBM2Files
NCO2=newNumericalComputationOptions(theDir2,P2)

--GED degree of quartic is 10. 
help EuclideanDistanceDegree
homotopyEDDegree(NCO2,"Weight",true,true)                    







restart
loadPackage"EuclideanDistanceDegree"
R=QQ[m1,m2,m3,u1,u2,u3,x,y,z,w1,w2,w3]
xList={x,y,z}
mList={m1,m2,m3}
mList={u1,u2,u3}
f=x^2+y^2-4*z^2
F={f}
Z=radical  (ideal(sum apply(xList,i->i^2))+ideal F)
bertiniPosDimSolve flatten entries gens Z-- two lines
determinantalGenericEuclideanDistanceDegree F
determinantalUnitEuclideanDistanceDegree F
--GED--generic m. 
genericM=matrix{{(1+m1)*(x-u1),(1+m2)*(y-u2),(1+m3)*(z-u3)}}||matrix makeJac(F,xList)
fixGEDProblem={u1=>1,u2=>1,u3=>1}|{m1=>12,m2=>124,m3=>99}

gedI=saturate((ideal F +minors(2,sub(genericM,fixGEDProblem))),ideal(x,y,z))
4==degree gedI
gedE=(eliminate({y,z},gedI))_*//first
--UED
fixUEDProblem={m1=>0,m2=>0,m3=>0}
uedI=saturate((ideal F +minors(2,sub(genericM,fixUEDProblem))),ideal(x,y,z))
2==degree sub(uedI,{u1=>134,u2=>3414,u3=>13444})
uedE=(eliminate({y,z},uedI))_*//first


wM=matrix{{w1*(x-u1),w2*(y-u2),w3*(z-u3)}}||matrix makeJac(feg Z,xList)
zedI=saturate((Z+minors(codim Z+1,wM)),ideal(x,y,z))

degree sub( saturate((Z+minors(codim Z+1,wM)),ideal(x,y,z)),fixUEDProblem|{w1=>134,w2=>2355,w3=>455})
zedE=(eliminate({y,z},zedI))_*//first

doubleProblem=ideal mingens ideal last coefficients(zedE*uedE-gedE,Variables=>{x})
codim 
printGens ideal mingens sub(doubleProblem,{u2=>124,u3=>35,w1=>342,w2=>242})


bigI=uedI+gedI+zedI
decompose bigI

coefficients((eliminate({z,y},saturate(uedI,radical ideal singularLocus ideal F)))_0,Variables=>{x})
coefficients((eliminate({z,y},saturate(gedI,radical ideal singularLocus ideal F)))_0,Variables=>{x})


g1=(eliminate({z,y},saturate(uedI,radical ideal singularLocus ideal F)))_0
g2=(eliminate({z,y},saturate(gedI,radical ideal singularLocus ideal F)))_0

S=frac (QQ[m1,m2,m3,u1,u2,u3])[x,y,z]
coefficientRing S
h1=sub(g1,S)
h2=sub(g2,S)
win=h2%ideal h1
toString (last coefficients win)_(0,0)


T=QQ[m1,m2,m3,u1,u2,u3]
umI=gcd(sub((last coefficients win)_(0,0),T),sub((last coefficients win)_(1,0),T))


--decompose umI

apropos"rder"
help MonomialOrderin
leadCoefficient(g1)
g2% g1



matrix transpose apply(xList,v->{v*(1+mLis),random(1,1000)})  )||matrix makeJac(F,xList)
gedProblem=
ideal F

minors(#F+2,(
    
    )

win=ideal F+ 


----Let's consider 2x2 rank one matrices
restart
loadPackage"EuclideanDistanceDegree"
R=QQ[x1,x2,x3,x4]
M=transpose genericMatrix(R,2,2)
F={det M}--(generic,unit)=(6,2)
bertiniPosDimSolve F--projective dim 2 and degree is 2
Q=ideal sum apply(gens R,i->i^2)
sl= ideal mingens ideal singularLocus(Q+ideal F)
decSL=decompose sl;
 netList decSL
decSL/codim
decSL/degree 
--singular locus is four projective points 
bertiniPosDimSolve flatten entries gens sl--4 projective points. 
#F==1
P=(F,F,{})
theDir
NCO=newNumericalComputationOptions(theDir,P)
help EuclideanDistanceDegree
homotopyEDDegree(NCO,"Weight",true,true)   --(6,2)                 

cp=importSolutionsFile(theDir,NameSolutionsFile=>"criticalPointFile")
netList cp
---
R=QQ[m0,m1,m2,m3,m4,u1,u2,u3,u4,x1,x2,x3,x4,w1,w2,w3,w4,s,t]
xList={x1,x2,x3,x4}
mList={m1,m2,m3,m4}
uList={u1,u2,u3,u4}
wList={w1,w2,w3,w4}
Q=ideal sum apply(xList,i->i^2)
f=x1*x2-x3*x4
F={f}
Z=radical  (ideal(sum apply(xList,i->i^2))+ideal F)
bertiniPosDimSolve flatten entries gens Z-- two lines
determinantalGenericEuclideanDistanceDegree F
determinantalUnitEuclideanDistanceDegree F
--GED--generic m. 

M=matrix{uList}||matrix{wList}*diagonalMatrix xList||matrix makeJac(F,xList)
--B=matrix{}||matrix{mList}||matrix{wList}
wSub=apply(mList,uList/(i->1),(a,b)->a*s+m0*b*t);wSub=apply(wList,wSub,(i,j)->i=>j)
wSub
I=minors(codim ideal F+2,sub(M,wSub))
I=saturate(I,ideal(t,s))
I=saturate(I+ideal F,ideal(t,s))
I=saturate(I,ideal xList);
I=saturate(I,ideal {m0,m1,m2,m3,m4});
I=saturate(I,ideal {m0});

--uSub=uList/(i->i=>random(1,100))
uSub={u1 => 32, u2 => 28, u3 => 13, u4 => 29}
E=eliminate({x1,x2},sub(I,uSub))
printGens sub(E,{s=>0})
(t)^4*(x3-x4)^2*(x3+x4)^2*(m0)^4*(99161*x3^2-226802*x3*x4+99161*x4^2)*(16)


I2= sub(I,{s=>0,t=>1})
A=saturate(I2,m0)
decA=decompose A
apply(decA,i->i+Q==i)--two of the three components are in the isotrpica quadric
twoPointsA= sub(first decA,uSub)
3==codim twoPointsA
2==degree twoPointsA
eliminate({x1,x2},twoPointsA)


twoPointsStandard=first decompose (minors(3,sub(M,uSub|(wList/(i->i=>1))))+ideal F)
codim twoPointsStandard
2===degree twoPointsStandard
twoPointsStandard==twoPointsA







---Let's try 2x3 matrices. Want positive dimensional Z.
---
R=QQ[m0,m1,m2,m3,m4,m5,m6,u1,u2,u3,u4,u5,u6,x1,x2,x3,x4,x5,x6,w1,w2,w3,w4,w5,w6,s,t]
xList={x1,x2,x3,x4,x5,x6}
mList={m1,m2,m3,m4,m5,m6}
uList={u1,u2,u3,u4,u5,u6}
wList={w1,w2,w3,w4,w5,w6}
Q=ideal sum apply(xList,i->i^2)
F=flatten entries gens minors(2,matrix{{x1,x2,x3},{x4,x5,x6}})
Z=radical  (ideal(sum apply(xList,i->i^2))+ideal F)

bertiniPosDimSolve flatten entries gens Z-- #gens R-24==3--codimension is 3 and projective dimension of Z is 5-3
(decompose Z)/gens/entries/flatten/bertiniPosDimSolve
--We have 2 lines and a quartic
determinantalGenericEuclideanDistanceDegree F--10
determinantalUnitEuclideanDistanceDegree F--2
--GED--generic m. 

M=matrix{uList}||matrix{wList}*diagonalMatrix xList||matrix makeJac(F,xList)
--B=matrix{}||matrix{mList}||matrix{wList}
wSub=apply(mList,uList/(i->1),(a,b)->a*s+m0*b*t);wSub=apply(wList,wSub,(i,j)->i=>j)
wSub
uSub={u1 => 32, u2 => 28, u3 => 13, u4 => 29}
I=sub(minors(codim ideal F+2,sub(M,wSub)),uSub);
I=ideal mingens ideal(gens I % ideal F)
I=sub(I,{m0=>1})
printGens I

I=saturate(I,ideal(t,s));
I=saturate(I+ideal F,ideal(t,s));
I=saturate(I,ideal xList);
I=saturate(I,ideal {m0,m1,m2,m3,m4});
I=saturate(I,ideal {m0});

--uSub=uList/(i->i=>random(1,100))

E=eliminate({x1,x2},sub(I,uSub))
printGens sub(E,{s=>0})
(t)^4*(x3-x4)^2*(x3+x4)^2*(m0)^4*(99161*x3^2-226802*x3*x4+99161*x4^2)*(16)


I2= sub(I,{s=>0,t=>1})
A=saturate(I2,m0)
decA=decompose A
apply(decA,i->i+Q==i)--two of the three components are in the isotrpica quadric
twoPointsA= sub(first decA,uSub)
3==codim twoPointsA
2==degree twoPointsA
eliminate({x1,x2},twoPointsA)


twoPointsStandard=first decompose (minors(3,sub(M,uSub|(wList/(i->i=>1))))+ideal F)
codim twoPointsStandard
2===degree twoPointsStandard
twoPointsStandard==twoPointsA



----Max's prediction. 
restart
printingPrecision =100
loadPackage"EuclideanDistanceDegree"
help EuclideanDistanceDegree
R=QQ[x1,x2,x3,x4,x5,x6]
M=transpose genericMatrix(R,3,2)
F=flatten entries gens minors(2,M)
G=drop(F,1)
P=(F,G,{})
theDir=temporaryFileName();mkdir theDir
NCO=newNumericalComputationOptions(theDir,P)
fixM={.935612+.780809*ii, .71813+.874296*ii, .056595+.850869*ii, .980549+.247701*ii, .595609+.191267*ii, .02298+.610646*ii}
fixW=apply(fixM,i->i+1)
fixV={.72466+.412908*ii, .20736+.413345*ii, .161095+.942996*ii, .653861+.917835*ii, .621858+.739543*ii, .22654+.870631*ii}
fixU=apply(#fixW,(i)->fixV_i/fixW_i*fixM_i)
--(1+m)*u=m*u
NCO#"TargetWeight"=fixW
NCO#"StartData"=fixU

homotopyEDDegree(NCO,"Weight",true, true)
elevenPoints=importSolutionsFile(theDir,NameSolutionsFile=>"criticalPointFile")
#elevenPoints==11

Z=radical ideal mingens ideal singularLocus (ideal F+ideal sum apply (gens R,i->i^2))
printGens Z
--bertiniPosDimSolve flatten entries gens Z
--      dim 2:  (dim=2,deg=2) (dim=2,deg=2)
--0==determinantalUnitEuclideanDistanceDegree flatten entries gens Z
--8==determinantalGenericEuclideanDistanceDegree flatten entries gens Z

4==codim Z
printGens Z
zg=Z_*
GZ=ideal {zg_0,zg_2,zg_4,zg_8}
codim Z==codim GZ
decompose GZ
PZ=(Z_*,GZ_*,{})
theDir=temporaryFileName();mkdir theDir
ZNCO=newNumericalComputationOptions(theDir,PZ)
ZNCO#"StartWeight"=fixM
ZNCO#"StartData"=fixV
homotopyEDDegree(ZNCO,"Weight",true, false)

eightPoints=importSolutionsFile(theDir,NameSolutionsFile=>"filterFile");
#eightPoints==8
pickOneOfEight=(i,eightPoints)->(
    onePoint=eightPoints_i;
    onePoint=drop(onePoint,#(GZ_*)+1);
    onePoint=(1/first onePoint)*onePoint    )
onePoint=pickOneOfEight(1,eightPoints)

--drop 5==#(GZ_*)+1 coordinates 
first\ apply(elevenPoints,i->((1/i_(#G+1))*drop(i,#G+1)))
netList\\apply(elevenPoints,i->((1/i_(#G+1))*drop(i,#G+2)))
onePoint

---Let's try symbolic computation.
R=QQ[s,t,x1,x2,x3,x4,x5,x6,m0]
xList={x1,x2,x3,x4,x5,x6}
M=matrix{{x1,x2,x3},{x4,x5,x6}}
F=flatten entries gens minors(2,M)
wFix=apply(xList,i->random(1,100))
uFix=apply(xList,i->random(1,100))
mFix={53, 67, 83, 1, 64, 96}

critM=(matrix{uFix}||
matrix{apply(#xList,i->(mFix_i*s+t)*xList_i)}||
matrix makeJac(F,xList))

Z=radical ideal mingens ideal singularLocus (ideal F+ideal sum apply (xList,i->i^2))
zg=Z_*
GZ=ideal {zg_0,zg_2,zg_4,zg_8}

critZ=Z+minors(2+codim Z,
    matrix{uFix}||(matrix{apply(#xList,i->(mFix_i*s+t*m0)*xList_i)}
	)||matrix makeJac(zg,xList));
decCZ=decompose ideal mingens critZ;
decCZ/degree
last decCZ
printGens last decCZ
printGens 
--degree first decCZ
decCZ/codim

intersectZF=ideal mingens(minors(codim ideal F+2,critM)+ideal F    +last decCZ)
decompose radical intersectZF
degree first oo
