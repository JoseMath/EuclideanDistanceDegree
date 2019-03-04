--restart
--Projective formulation for intersections with linear spaces
rand:=randomValue
--Assume ring is a complex inexact field
--G is a subset of F. 
NumericalComputationOptions=new Type of MutableHashTable

  
    
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
nocKeys=parameterKeys|jacKeys|modelKeys|degreeKeys|bertiniKeys|coordinateKeys|directoryKeys|solutionKeys

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
    NCO#"TargetWeight"=NCO#"StartWeight"=apply(numX,i->random CC); 
    NCO#"GammaVector"=apply(numX-1,i->random CC); 
    scan(bertiniKeys,i->NCO#i={});
    NCO#"HomogeneousVariableGroups"={gens ring first F};
    NCO#"AffineVariableGroups"={};    
    return NCO
   )



defaultMTName="input_MT_"
startEDDegree=method()
startEDDegree(NumericalComputationOptions,Sequence,ZZ):=(NCO,TP,stage)->(    
--TP=type-parameters is a length three sequence with entries in {0,1,null}
--stage is either 1 or 2.
--Organize parameters for the homotopy. 
    theDir:=NCO#"Directory";
    (pvM,pvD,pvW):=TP;--weight is last
    pgStart:=TP/(i->if i===null then 0 else i)//toList;
    pgTarget:=TP/(i->if i===null then 1 else i)//toList;
    if stage==2 
    then (writeParameterFile(theDir,pgStart,NameParameterFile=>"start_parameters");
    	writeParameterFile(theDir,pgTarget,NameParameterFile=>"final_parameters"));
    PG:={"TModel","TData","TWeight"};
    if stage==1 then pgMT:=pgStart
    else if stage==2 then   pgMT=pgTarget
    else error"stage is either 1 or 2";
    print pgMT;
--Other strings that are used
    (lagMult,numerHB,denomQ,primal):=("lagMult","numerHB","denomQ","primal")   ; 
--Now we extract information from NCO.
--First the model and submodel constraints
------F is the model, V(G)\cap V(L) is a complete intersection contained in V(F)\cap V(L).
    (F,G,startL,targetL):=(NCO#"Model",NCO#"WitnessModel",NCO#"StartSubmodel",NCO#"TargetSubmodel");
    (startData,startWeight,targetData,targetWeight):=(NCO#"StartData",NCO#"StartWeight",NCO#"TargetData",NCO#"TargetWeight");
    (jacStartL,jacTargetL,jacG):=(NCO#"JacobianStartSubmodel",NCO#"JacobianTargetSubmodel",NCO#"JacobianWitnessModel");
    nc:=numcols jacStartL;
    randomGamma:=NCO#"GammaVector";
    xList:=NCO#"PrimalCoordinates";
-- --Homogenize appropriately
    (degSubmodel,degWitnessModel):=(NCO#"DegreeSubmodel",NCO#"DegreeWitnessModel");
    maxDegree:=(degSubmodel|degWitnessModel|{3})//max;
    degRescale:=({3}|degSubmodel|degWitnessModel)/(i->maxDegree-i);        
    print 1;
--Set L
    if pvM===null
    then L:=apply(startL,targetL,(i,j)->makeB'Section({i,j},B'NumberCoefficients=>{"1-"|PG_0,PG_0}))
    else if pvM===0 then L=startL
    else if pvM===1 then L=targetL;
    print 7;
--We have three rings. One ring for manipulating the Jacobian, one ring for manipulating the gradient, and one for manipulating both
--Gradient ring
    print 2;
    kk2:=ring first startData;
    topS:=kk2[apply(nc,i->primal|i)]**kk2[{lagMult|"0",numerHB,denomQ}|PG];
    primalList:=flatten entries basis({1,0},topS);
    primalSub:=transpose{primalList,xList};
    (topL0,topNumerHB,topDenomQ,topTModel,topTData,topTWeight):=toSequence flatten entries basis({0,1},topS);
--Tracking  subs
    pgSub:={};
    if pvD=!=null and stage==2 then pgSub=pgSub|{topTData   => pgMT_1};
    if pvW=!=null and stage==2 then pgSub=pgSub|{topTWeight => pgMT_2};
    if stage==1 then pgSub={topTData=>pgMT_1,topTWeight=>pgMT_2};    
    print pgSub;
--MT  subs
    pgMTSub:={topTData=>pgMT_1,topTWeight=>pgMT_2};
    print pgMTSub;
--Gradient conditions (RHS)
----
    if pvW===null and stage==2
    
    then isoQ:=makeB'Section(
    	sum apply(nc,i->{startWeight_i*(primalList_i)^2*,targetWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{"(1-"|PG_2|")*"|lagMult|"0*"|numerHB|"^"|first degRescale,"("|PG_2|")*"|lagMult|"0*"|numerHB|"^"|first degRescale},
	NameB'Section=>denomQ)
    else if pvW===0 or (pvW===null and stage==1) then isoQ=makeB'Section(
    	sum apply(nc,i->{startWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ)
    else if pvW===1  then isoQ=makeB'Section(
    	sum apply(nc,i->{targetWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ)
    else error"isoQ wrong";	
    if (pvW===null and stage==1) or pvW===0 
    then isoQMT:=makeB'Section(
    	sum apply(nc,i->{startWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ)
    else if (pvW===null and stage==2) or pvW===1 then isoQMT=makeB'Section(
    	sum apply(nc,i->{targetWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ)
    else error"isoQMT is wrong.";

    RHS:=apply(nc,i->"RHS"|i=>
	topL0*topNumerHB^(first degRescale)*sub((((1-topTData)*topDenomQ*startData_i-(1-topTWeight)*topNumerHB*primalList_i*startWeight_i+
	(topTData)*topDenomQ*targetData_i-(topTWeight)*topNumerHB*primalList_i*targetWeight_i)),pgSub));
    print RHS;
-*
    hyperplaneAtInfinity:=sub(
	sum apply(nc,i->(1-topTData)*startData_i*primalList_i+
	    (topTData)*targetData_i*primalList_i),pgSub);
    hyperplaneAtInfinityMT:=sub(hyperplaneAtInfinity,pgMTSub);
    print hyperplaneAtInfinity;
    print hyperplaneAtInfinityMT;
    print pgSub;
*-
--    
    if pvD===null and stage==2
    then hyperplaneAtInfinity:=makeB'Section(
    	sum apply(nc,i->{startData_i*(primalList_i),targetData_i*(primalList_i)}),
	B'NumberCoefficients=>{1-topTData,topTData},
	NameB'Section=>numerHB)
    else if pvD===0 or (pvW===null and stage==1) then hyperplaneAtInfinity=makeB'Section(
    	sum apply(nc,i->{startData_i*(primalList_i)}),
	B'NumberCoefficients=>{1},
	NameB'Section=>numerHB)
    else if pvD===1 then hyperplaneAtInfinity=makeB'Section(
    	sum apply(nc,i->{targetData_i*(primalList_i)}),
	B'NumberCoefficients=>{1},
	NameB'Section=>numerHB)
    else error"hyperplaneAtInfinity wrong";	
    if (pvD===null and stage===1) or pvD===0 then hyperplaneAtInfinityMT:=makeB'Section(
    	sum apply(nc,i->{startData_i*(primalList_i)}),
	B'NumberCoefficients=>{1},
	NameB'Section=>numerHB)
    else if (pvD===null and stage===2) or pvD==1 then hyperplaneAtInfinityMT=makeB'Section(
    	sum apply(nc,i->{targetData_i*(primalList_i)}),
	B'NumberCoefficients=>{1},
	NameB'Section=>numerHB)
    else error"hyperplaneAtInfinityMT is wrong";
----ISOTROPIC Q
    if pvW===null and stage==2
    then isoQ:=makeB'Section(
    	sum apply(nc,i->{startWeight_i*(primalList_i)^2,targetWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1-topTWeight,topTWeight},
	NameB'Section=>denomQ)
    else if pvW===0 or (pvW===null and stage==1) then isoQ=makeB'Section(
    	sum apply(nc,i->{startWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ)
    else if pvW===1  then isoQ=makeB'Section(
    	sum apply(nc,i->{targetWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ)
    else error"isoQ wrong";	
    if (pvW===null and stage==1) or pvW===0 
    then isoQMT:=makeB'Section(
    	sum apply(nc,i->{startWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ)
    else if (pvW===null and stage==2) or pvW===1 then isoQMT=makeB'Section(
    	sum apply(nc,i->{targetWeight_i*(primalList_i)^2}),
	B'NumberCoefficients=>{1},
	NameB'Section=>denomQ)
    else error"isoQMT is wrong.";
--Jacobian ring
    theR:=ring first F;
    numX:=#gens theR;
    kk1:=coefficientRing ring first F;
    jacS:=theR**kk1[apply(#startL+#G,i->lagMult|i+1)]**kk1[{numerHB,denomQ}|PG];
    jacLamList:=flatten entries basis({0,1,0},jacS);
    (jacNumerHB,jacDenomQ,jacTModel,jacTData,jacTWeight):=toSequence flatten entries basis({0,0,1},jacS);
    jacLV:=apply(jacLamList,drop(degRescale,1),(lam,j)->if j==0 
	then lam else if j>0 then lam*jacNumerHB^j 
	else print "Error: Homogenized incorrectly.");
    --print LV; 
    print 7;
    startJacStartCondition:=flatten entries (matrix{jacLV}*sub((jacStartL||jacG),jacS));    
    startLHS:=apply(nc,i->"SLHS"|i=>startJacStartCondition_i);
    targetJacCondition:=flatten entries (matrix{jacLV}*sub((jacTargetL||jacG),jacS));    
    targetLHS:=apply(nc,i->"TLHS"|i=>targetJacCondition_i);
--Jacobian definition conditions (LHS) 
    if  #startL>0 and pvM===null    
    then LHS:=apply(nc,i->"LHS"|i=>"SLHS"|i|"*(1-TModel)+TLHS"|i|"*(TModel)")
    else if pvM===0 then LHS=apply(nc,i->"LHS"|i=>"SLHS"|i)
    else if pvM===1 then LHS=apply(nc,i->"LHS"|i=>"TLHS"|i)
    else error"last argument is a sequence of three elements in {0,1,null}^3";
    print 8;
---Critical ring (LHS-RHS)
    kk3:=ring first randomGamma    ;
    critRing:=kk3[apply(nc,i->"critCondition"|i)];
    critConditions:=apply(nc,i->(gens critRing)_i=>"LHS"|i|"-RHS"|i);
    print critConditions;
    randomizedCritConditions:=apply(drop(gens critRing,-1),randomGamma,(i,j)->i+j*last gens critRing);
----Input file for stage 1 and stage 2.
    win:=L|G|randomizedCritConditions;
    bfs:=NCO#"BertiniSubstitute"|primalSub|{hyperplaneAtInfinity,isoQ}|startLHS|targetLHS|RHS|LHS|critConditions;
    print bfs;
    print 12;
    if stage==1 then PG={"dsadadf"};
    if stage==1 
    then nif:="inputCriticalPointSuperSet" 
    else if stage==2 
    then nif="inputParameterHomotopySuperSet";
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
	B'Constants=>NCO#"BertiniConstants"
	);
    print 14;
-----Membership test files function (s,k,bp)=(string to name the file, tracktype, list of polynomials)
    bfs2:=NCO#"BertiniSubstitute"|primalSub|{hyperplaneAtInfinityMT,isoQMT}|startLHS|targetLHS|RHS|LHS|critConditions;    
    print bfs2;
    imt:=(s,k,bp)->makeB'InputFile(theDir,NameB'InputFile=>(defaultMTName|s|toString k),
	--Which variables come first?
	AffVariableGroup=>flatten flatten{NCO#"HomogeneousVariableGroups",{{topL0}|jacLamList},NCO#"AffineVariableGroups"},
	B'Configs=>{"TrackType"=>k,"PrintPathProgress"=>1000}|(NCO#"BertiniMembershipTestConfiguration"),
	B'Polynomials=>bp,
	B'Functions=>bfs2,
	B'Constants=>pgMTSub|NCO#"BertiniConstants"
	);
--Filter Residuals
    imt("Residual"|stage|"_",1,{first last critConditions});
    imt("Residual"|stage|"_",3,{first last critConditions});
    print numerHB;
--UBeta
    imt("Degenerate"|stage|"_",1,{numerHB|"*"|denomQ|"*"|lagMult|"0"});
    imt("Degenerate"|stage|"_",3,{numerHB|"*"|denomQ|"*"|lagMult|"0"});
    print hyperplaneAtInfinity;
    print hyperplaneAtInfinityMT;
--Filer component    
    scan(#F,i->(
    	    imt("Hypersurface"|stage|"_"|i|"_",1,{F_i});
    	    imt("Hypersurface"|stage|"_"|i|"_",3,{F_i});
        ))  )


   
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
---Run membershipTest
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
 
filterSolutionFile=method()
filterSolutionFile(NumericalComputationOptions,String,ZZ):=(NCO,newFileName,numCoords)->(     
    theDir:=NCO#"Directory";
    if fileExists(addSlash(theDir)|newFileName)
    then removeFile(addSlash(theDir)|newFileName);    
    selectPoints:=NCO#"TrackSolutions";
    numSols   := #selectPoints;
    firstLine := true;
    countSol  := -1;
    countLine := 0;
    groupSize := 1+numCoords;
    isSelected:= null;
------
    sf:=openOutAppend(addSlash(theDir)|newFileName);
    scanLineSolutionFunction := (ell)->(
--      print (countLine,countSol,"every");
      if firstLine 
      then (firstLine=false; sf<< toString(sum selectPoints)<<endl)
      else if countSol < numSols 
      then (
    	  if countLine==0 then isSelected=selectPoints_countSol;
	  countLine=countLine+1;
    	  if isSelected==1 then sf <<ell<<endl;
      	  if countLine==groupSize 
      	  then (
	      --print (countLine,groupSize,"grp");
	      countLine=0; 
	      countSol=countSol+1;
	      )));
      scanLines(scanLineSolutionFunction,addSlash(theDir)|"member_points");      
      close sf;
      return (theDir,newFileName))


(stageOne,stageTwo)=(1,2);
numericWeightEDDegree=method()
numericWeightEDDegree(String,Sequence):=(theDir,P)->(    
    WV:=apply(#gens ring first first P,i->random CC);
    NCO:=newNumericalComputationOptions(theDir,P);
    numericWeightEDDegree(NCO,WV)
    )

numericWeightEDDegree(String,Sequence,List):=(theDir,P,WV)->(    
    NCO:=newNumericalComputationOptions(theDir,P);
    numericWeightEDDegree(NCO,WV)
    )
numericWeightEDDegree(NumericalComputationOptions,List):=(NCO,WV)->(
    NCO#"StartWeight"=WV;
    homotopyType:=(0,0,0);
    startEDDegree(NCO,homotopyType,stageOne);
    return runBertiniStartEDDegree(NCO,homotopyType,stageOne)    
    )


numericUnitEDDegree=method()
numericUnitEDDegree(String,Sequence):=(theDir,P)->(
    NCO:=newNumericalComputationOptions(theDir,P);
    numericUnitEDDegree(NCO))
numericUnitEDDegree(NumericalComputationOptions):=(NCO)->(
    WV:=apply(#gens ring first (NCO#"Model"),i->1);
    numericWeightEDDegree(NCO,WV))


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
