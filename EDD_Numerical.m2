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
    if #P==1 then (F,G,L)=(P_0,P_0,{});
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
startEDDegree(NumericalComputationOptions,Sequence,ZZ):=(NCO,parameters,stage)->(    
    (pvM,pvD,pvW):=parameters;
    pgStart:=toList apply(parameters,i->if i===null then 0 else i);
    pgTarget:=toList apply(parameters,i->if i===null then 1 else i);
    --
    writeParameterFile(storeBM2Files,pgStart,NameParameterFile=>"start_parameters");
    writeParameterFile(storeBM2Files,pgTarget,NameParameterFile=>"target_parameters");
--F is the model, V(G)\cap V(L) is a complete intersection contained in V(F)\cap V(L)
    (theDir,F,G,L):=(NCO#"Directory",NCO#"Model",NCO#"WitnessModel",NCO#"StartSubmodel");
    (data,weight,randomGamma):=(NCO#"StartData",NCO#"StartWeight",NCO#"GammaVector");
    (targetData,targetWeight):=(NCO#"TargetData",NCO#"TargetWeight");
    theR:=ring first F;
    numX:=#gens theR;
    kk1:=coefficientRing ring first F;
    print 1;
    jacS:=theR**kk1[apply(#L+#G,i->"L"|i+1)]**kk1["numerHB","denomQ","TData","TWeight"];
    jacLamList:=flatten entries basis({0,1,0},jacS);
    (jacStartL,jacTargetL,jacG):=(NCO#"JacobianStartSubmodel",NCO#"JacobianTargetSubmodel",NCO#"JacobianWitnessModel");
    nc:=numcols jacStartL;
    print 3;
    (jacNumerHB,jacDenomQ,jacTData,jacTWeight):=toSequence flatten entries basis({0,0,1},jacS);
    print 3;
--
    kk2:=ring first data;
    topS:=kk2[apply(nc,i->"primal"|i)]**kk2["L0","numerHB","denomQ","TData","TWeight"];
    primalList:=flatten entries basis({1,0},topS);
    primalSub:=transpose{primalList,NCO#"PrimalCoordinates"};
    print 4;
    (topL0,topNumerHB,topDenomQ,topTData,topTWeight):=toSequence flatten entries basis({0,1},topS);
    print 5;
--homogenize jacobian. 
    degRescale:={3}|NCO#"DegreeSubmodel"|NCO#"DegreeWitnessModel";
    maxDeg:=(max degRescale);
    --print degRescale; 
    degRescale=apply(degRescale,i->maxDeg-i);
    --print degRescale;
    jacLV:=apply(jacLamList,drop(degRescale,1),(lam,j)->if j==0 
	then lam else if j>0 then lam*jacNumerHB^j 
	else print "Error: Homogenized incorrectly.");
    --print LV; 
    print 6;
    startJacStartCondition:=flatten entries (matrix{jacLV}*sub((jacStartL||jacG),jacS));    
    startLHS:=apply(nc,i->"SLHS"|i=>startJacStartCondition_i);
    targetJacCondition:=flatten entries (matrix{jacLV}*sub((jacTargetL||jacG),jacS));    
    targetLHS:=apply(nc,i->"TLHS"|i=>targetJacCondition_i);
    print 4;
--Model conditions. 
    if #L>0 or  pvM===null    
    then LHS:=apply(nc,i->"LHS"|i=>"SLHS"|i|"*(1-TModel)+TLHS"|i|"*(TModel)")
    else if pvM===0 then LHS=apply(nc,i->"LHS"|i=>"SLHS"|i)
    else if pvM===1 then LHS=apply(nc,i->"LHS"|i=>"TLHS"|i)
    else error"last argument is a sequence of three elements in {0,1,null}^3";
    print 5;
--TopRow conditios
    pSub:={};
    print pSub;
    print (topTData=>pvD);
    if pvD=!=null then pSub=pSub|{topTData=>pvD};
    print 5;
    if pvW=!=null then pSub=pSub|{topTWeight=>pvD};
    RHS:=apply(nc,i->"RHS"|i=>topL0*topNumerHB^(first degRescale)*sub((((1-topTData)*topDenomQ*data_i-(1-topTWeight)*topNumerHB*primalList_i*weight_i+
	(topTData)*topDenomQ*targetData_i-(topTWeight)*topNumerHB*primalList_i*targetWeight_i)),pSub));
    print 7;
    print RHS;
    hyperplaneAtInfinity:=sub(sum apply(nc,i->(1-topTData)*data_i*primalList_i+(topTData)*targetData_i*primalList_i),pSub);
    isoQ:=sub(sum apply(nc,i->(1-topTWeight)*weight_i*(primalList_i)^2+(topTWeight)*targetWeight_i*(primalList_i)^2),pSub);
    kk3:=ring first randomGamma    ;
    print 8;
    critRing:=kk3[apply(nc,i->"critCondition"|i)];
    print 9;
    critConditions:=apply(nc,i->(gens critRing)_i=>"LHS"|i|"-RHS"|i);
    print 10;
    randomizedCritConditions:=apply(drop(gens critRing,-1),randomGamma,(i,j)->i+j*last gens critRing);
    win:=L|G|randomizedCritConditions;
----Input file 
    bfs:=NCO#"BertiniSubstitute"|primalSub|{topNumerHB=>hyperplaneAtInfinity,topDenomQ=>isoQ}|startLHS|targetLHS|RHS|LHS|critConditions;
    print 12;
    pgv:={"TModel","TData","TWeight"};
    makeB'InputFile(theDir,NameB'InputFile=>"inputCriticalPointSuperSet",
	HomVariableGroup=>(NCO#"HomogeneousVariableGroups")|{{topL0}|jacLamList},
    	AffVariableGroup=>NCO#"AffineVariableGroups",
	B'Configs=>{
	    "UseRegeneration"=>1,
	    "TrackType"=>0,
	    "ParameterHomotopy"=>stage,
	    "PrintPathProgress"=>1000}|(NCO#"BertiniStartFiberSolveConfiguration"),
	B'Polynomials=>win,
    	ParameterGroup=>pgv,
	B'Functions=>bfs,
	B'Constants=>NCO#"BertiniConstants"
	);
    print 14;
-----Membership test files function (s,k,bp)=(string to name the file, tracktype, list of polynomials)
    bfs2:=NCO#"BertiniSubstitute"|primalSub|{topNumerHB=>hyperplaneAtInfinity,topDenomQ=>isoQ}|startLHS|targetLHS|RHS|LHS|critConditions;    
    if stage==1 then pgStageConstants:=transpose{pgv,pgStart}
    else if stage==2 then   pgStageConstants=transpose{pgv,pgTarget}
    else error"stage is either 1 or 2";
    imt:=(s,k,bp)->makeB'InputFile(theDir,NameB'InputFile=>(defaultMTName|s|toString k),
	--Which variables come first?
	AffVariableGroup=>flatten flatten{NCO#"HomogeneousVariableGroups",{{topL0}|jacLamList},NCO#"AffineVariableGroups"},
	B'Configs=>{"TrackType"=>k,"PrintPathProgress"=>1000}|(NCO#"BertiniMembershipTestConfiguration"),
	B'Polynomials=>bp,
	B'Functions=>bfs2,
	B'Constants=>(pgStageConstants)|NCO#"BertiniConstants"
	);
--Filter Residuals
    imt("Residual",1,{first last critConditions});
    imt("Residual",3,{first last critConditions});
--UBeta
    imt("Degenerate",1,{"numerHB*denomQ*L0"});
    imt("Degenerate",3,{"numerHB*denomQ*L0"});
--Filer component    
    scan(#F,i->(
    	    imt("Hypersurface_"|i|"_",1,{F_i});
    	    imt("Hypersurface_"|i|"_",3,{F_i});
        ))  )


   
-*
restart--
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x1,x2,x3,x4]
F=flatten entries gens minors(2,transpose genericMatrix(R,2,2))
G=F
NCO=newNumericalComputationOptions(storeBM2Files,(F,G))
startEDDegree(NCO,(0,0,0),1)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)
NCO#"StartWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,0,0),1)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)
NCO#"TrackSolutions"
peek NCO

*-

runBertiniStartEDDegree=method(Options=>{})
runBertiniStartEDDegree(String,ZZ,NumericalComputationOptions):= o->(storeBM2Files,n,NCO)->(
--    writeParameterFile(storeBM2Files,{0,0,0},NameParameterFile=>"start_parameters");
--    writeParameterFile(storeBM2Files,{1,1,1},NameParameterFile=>"final_parameters");
    runBertini(storeBM2Files,NameB'InputFile=>"inputCriticalPointSuperSet");
    print "inputCriticalPointSuperSet";
    moveB'File(storeBM2Files,"bertini_session.log","bertini_session_CriticalPointSuperSet.log",CopyB'File => true);
    moveB'File(storeBM2Files,"nonsingular_solutions","member_points");
---Run membershipTest
    runMT:=(s)->(
   	runBertini(storeBM2Files,NameB'InputFile=>defaultMTName|s|"1");
	print (defaultMTName|s|"1");
    	moveB'File(storeBM2Files,"bertini_session.log","bertini_session_MT_"|s|"1.log",CopyB'File => true);
    	runBertini(storeBM2Files,NameB'InputFile=>defaultMTName|s|"3");
	print (defaultMTName|s|"3");
    	moveB'File(storeBM2Files,"bertini_session.log","bertini_session_MT_"|s|"3.log",CopyB'File => true);
    	return importIncidenceMatrix(storeBM2Files));                                                                           
--Residual
    imResidual:=runMT("Residual");
--Degenerate
    imDegenerate:=runMT("Degenerate");
--Hypersurfaces
    imAllHypersurfaces:=apply(n,i->runMT("Hypersurface_"|i|"_"));
    EDDeg:=0;
    memberEveryHypersurface:=(i)->(
	output:=true;
	scan(imAllHypersurfaces,j->if {}==j_i then (output=false;break));
	return output);
    ts:={};
    scan(#imResidual,i->if imResidual_i=!={} and imDegenerate_i=={} and memberEveryHypersurface(i) 
	then (EDDeg=EDDeg+1; ts=ts|{1}) else ts=ts|{0});
    NCO#"TrackSolutions"=ts;
    moveB'File(storeBM2Files,"bertini_session_CriticalPointSuperSet.log","bertini_session.log");
--    return(imResidual,imDegenerate)
    return(EDDeg)
     )
 
filterSolutionFile=method()
filterSolutionFile(String,String,ZZ,NumericalComputationOptions):=(storeBM2Files,newFileName,numCoords,NCO)->(     
    if fileExists(addSlash(storeBM2Files)|newFileName)
    then removeFile(addSlash(storeBM2Files)|newFileName);    
    selectPoints:=NCO#"TrackSolutions";
    numSols   := #selectPoints;
    firstLine := true;
    countSol  := -1;
    countLine := 0;
    groupSize := 1+numCoords;
    isSelected:= null;
------
    sf:=openOutAppend(addSlash(storeBM2Files)|newFileName);
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
      scanLines(scanLineSolutionFunction,addSlash(storeBM2Files)|"member_points");      
      close sf;
      return (storeBM2Files,newFileName))

-*
restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x1,x2,x3,y1,y2,y3]
F=flatten entries gens minors(2,transpose genericMatrix(R,3,2))
G=drop(F,-1)
L={}
NCO=newNumericalComputationOptions(storeBM2Files,(F,G))
NCO#"StartWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,0,1),1)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)
filterSolutionFile(storeBM2Files,"start",#gens ring first F+1+#G+#L,NCO)     
readFile(storeBM2Files,"start",1000)
*-

-*---Homotopy
restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x1,x2,x3,x4]
F=flatten entries gens minors(2,transpose genericMatrix(R,2,2))
G=drop(F,0)
L={}

--HOMOTOPY
NCO=newNumericalComputationOptions(storeBM2Files,(F,G))
--STAGE ONE--10
stage=1
startEDDegree(NCO,(0,0,0),stage)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)
--FILTER
filterSolutionFile(storeBM2Files,"start",#gens ring first F+1+#G+#L,NCO)     
readFile(storeBM2Files,"start",1000)
-- STAGE TWO
stage=2
NCO#"TargetWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,0,null),stage)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)
*-

    
end

runBertiniWeightTargetEDDegree=method(Options=>{
	})
runBertiniWeightTargetEDDegree(String,List):= o->(storeBM2Files,weight)->(
    moveB'File(storeBM2Files,"weight_start_parameters","start_parameters",CopyB'File=>true);
    writeParameterFile(storeBM2Files,weight|{0},NameParameterFile=>"final_parameters");
    runBertini(storeBM2Files,NameB'InputFile=>"input_PH_weight");
--Residual
    imResidual:=runMT("Residual");
--Degenerate
    imDegenerate:=runMT("Degenerate");
--Hypersurfaces
    imAllHypersurfaces:=apply(n,i->runMT("Hypersurface_"|i|"_"));
    EDDeg:=0;
    memberEveryHypersurface:=(i)->(
	output:=true;
	scan(imAllHypersurfaces,j->if {}==j_i then (output=false;break));
	return output);
    scan(#imResidual,i->if imResidual_i=!={} and imDegenerate_i=={} and memberEveryHypersurface(i) then EDDeg=EDDeg+1);
    moveB'File(storeBM2Files,"bertini_session_CriticalPointSuperSet.log","bertini_session.log");
--    return(imResidual,imDegenerate,imComponent)
    return(EDDeg)
     )

end


 
---EXAMPLES

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
