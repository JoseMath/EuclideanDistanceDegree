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
startEDDegree(NumericalComputationOptions,Sequence,ZZ):=(NCO,TP,stage)->(    
--TP=type-parameters is a length three sequence with entries in {0,1,null}
--stage is either 1 or 2.
--Organize parameters for the homotopy. 
    theDir:=NCO#"Directory";
    (pvM,pvD,pvW):=TP;--weight is last
    pgStart:=TP/(i->if i===null then 0 else i)//toList;
    pgTarget:=TP/(i->if i===null then 1 else i)//toList;
    writeParameterFile(theDir,pgStart,NameParameterFile=>"start_parameters");
    writeParameterFile(theDir,pgTarget,NameParameterFile=>"final_parameters");
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
    RHS:=apply(nc,i->"RHS"|i=>
	topL0*topNumerHB^(first degRescale)*sub((((1-topTData)*topDenomQ*startData_i-(1-topTWeight)*topNumerHB*primalList_i*startWeight_i+
	(topTData)*topDenomQ*targetData_i-(topTWeight)*topNumerHB*primalList_i*targetWeight_i)),pgSub));
    print RHS;
    hyperplaneAtInfinity:=sub(
	sum apply(nc,i->(1-topTData)*startData_i*primalList_i+
	    (topTData)*targetData_i*primalList_i),pgSub);
    hyperplaneAtInfinityMT:=sub(hyperplaneAtInfinity,pgMTSub);
    print hyperplaneAtInfinity;
    print hyperplaneAtInfinityMT;
    print pgSub;
    isoQ:=sub(sum apply(nc,i->
	    (1-topTWeight)*startWeight_i*(primalList_i)^2+
	    (topTWeight)*targetWeight_i*(primalList_i)^2),pgSub);
    isoQMT:=sub(isoQ,pgMTSub);
    print isoQ;
    print isoQMT;
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
----Input file 
    win:=L|G|randomizedCritConditions;
    bfs:=NCO#"BertiniSubstitute"|primalSub|{topNumerHB=>hyperplaneAtInfinity,topDenomQ=>isoQ}|startLHS|targetLHS|RHS|LHS|critConditions;
    print 12;
    makeB'InputFile(theDir,NameB'InputFile=>"inputCriticalPointSuperSet",
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
    bfs2:=NCO#"BertiniSubstitute"|primalSub|{topNumerHB=>hyperplaneAtInfinityMT,topDenomQ=>isoQMT}|startLHS|targetLHS|RHS|LHS|critConditions;    
    imt:=(s,k,bp)->makeB'InputFile(theDir,NameB'InputFile=>(defaultMTName|s|toString k),
	--Which variables come first?
	AffVariableGroup=>flatten flatten{NCO#"HomogeneousVariableGroups",{{topL0}|jacLamList},NCO#"AffineVariableGroups"},
	B'Configs=>{"TrackType"=>k,"PrintPathProgress"=>1000}|(NCO#"BertiniMembershipTestConfiguration"),
	B'Polynomials=>bp,
	B'Functions=>bfs2,
	B'Constants=>pgMTSub|NCO#"BertiniConstants"
	);
--Filter Residuals
    imt("Residual",1,{first last critConditions});
    imt("Residual",3,{first last critConditions});
    print numerHB;
--UBeta
    imt("Degenerate",1,{numerHB|"*"|denomQ|"*"|lagMult|"0"});
    imt("Degenerate",3,{numerHB|"*"|denomQ|"*"|lagMult|"0"});
    print hyperplaneAtInfinity;
    print hyperplaneAtInfinityMT;
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
runBertiniStartEDDegree(ZZ,NumericalComputationOptions):= o->(n,NCO)->(
    theDir:=NCO#"Directory";
--    writeParameterFile(storeBM2Files,{0,0,0},NameParameterFile=>"start_parameters");
--    writeParameterFile(storeBM2Files,{1,1,1},NameParameterFile=>"final_parameters");
    runBertini(theDir,NameB'InputFile=>"inputCriticalPointSuperSet");
    print "inputCriticalPointSuperSet"; 
    moveB'File(theDir,"bertini_session.log","bertini_session_CriticalPointSuperSet.log",CopyB'File => true);
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
    imResidual:=runMT("Residual");
--Degenerate
    print "isDegenerate";
    imDegenerate:=runMT("Degenerate");
--Hypersurfaces
    print "isMemberOfVariety";
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
    moveB'File(theDir,"bertini_session_CriticalPointSuperSet.log","bertini_session.log");
--    return(imResidual,imDegenerate)
    return(EDDeg)
     )
 
filterSolutionFile=method()
filterSolutionFile(String,ZZ,NumericalComputationOptions):=(newFileName,numCoords,NCO)->(     
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

-*
restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x1,x2,x3,y1,y2,y3]
F=flatten entries gens minors(2,transpose genericMatrix(R,3,2))
G=drop(F,-1)
L={}
NCO=newNumericalComputationOptions(theDir,(F,G))
NCO#"StartWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,0,1),1)
runBertiniStartEDDegree(#F,NCO)
filterSolutionFile(theDir,"start",#gens ring first F+1+#G+#L,NCO)     
readFile("start",1000)
*-

-*
restart
printingPrecision=100
loadPackage("EuclideanDistanceDegree",Reload=>true)
R=QQ[x1,x2,x3,x4]
F=flatten entries gens minors(2,transpose genericMatrix(R,2,2))
G=drop(F,0)
L={}


--STAGE ONE TESTS
stageOne=1
stageTwo=2
NCO=newNumericalComputationOptions(storeBM2Files,(F,G,L))
NCO#"TargetWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,0,0),stageOne)
runBertiniStartEDDegree(#F,NCO)--6
readFile(storeBM2Files,"member_points",10000)
10                                               

      -6.289977362595892e-01 1.912844762422852e-01
      3.229330587679995e-01 6.032119326084742e-01
      -1.874164479402360e-01 9.424130040076552e-01
      1.000000000000000e+00 0.000000000000000e+00
      1.000000000000000e+00 0.000000000000000e+00
      1.649271782698879e-01 5.456541500970949e-01

      5.023820446423847e-01 4.548167910299095e-01
      -2.769246618401497e-01 1.103170119655528e-02
      1.000000000000000e+00 0.000000000000000e+00
      -2.920103070944934e-01 2.863217217871001e-01
      1.000000000000000e+00 0.000000000000000e+00
      -6.555591290157402e-01 -6.693750148030577e-01

      1.000000000000000e+00 0.000000000000000e+00
      1.686452067024666e-02 8.243968675810902e-01
      1.809342382215986e-01 -8.298517551782524e-01
      6.871785567260653e-01 1.351665671496124e-01
      1.000000000000000e+00 0.000000000000000e+00
      -3.983540147857419e-01 -7.286234417087349e-01

      1.000000000000000e+00 0.000000000000000e+00
      7.136091195660894e-02 4.911638365853592e-01
      -3.076656085394328e-01 6.399390901710897e-01
      -3.362702371124366e-01 -1.054475836043053e-01
      1.000000000000000e+00 0.000000000000000e+00
      9.889773344040313e-02 6.186440695556037e-01

      -8.727865247599269e-01 -4.238955238272207e-01
      2.493561439182775e-01 5.112914443215998e-01
      1.000000000000000e+00 0.000000000000000e+00
      -4.613858711962510e-01 -3.617288189127669e-01
      -9.167465833888927e-01 -6.267791280300083e-03
      1.000000000000000e+00 0.000000000000000e+00

      -6.666891893131216e-01 5.074842990005828e-01
      1.000000000000000e+00 0.000000000000000e+00
      1.651006729939873e-01 -4.238158851386274e-01
      -4.631679365721978e-01 2.831385787464482e-01
      -2.796203070225711e-01 2.400879156300438e-01
      1.000000000000000e+00 0.000000000000000e+00

      1.000000000000000e+00 0.000000000000000e+00
      -3.495186762173846e-01 -4.238479836156782e-01
      -2.195951495118892e-01 -8.881869498496547e-02
      3.910698118440107e-02 1.241187540268990e-01
      -5.570299013891629e-01 -3.123798700592145e-01
      1.000000000000000e+00 0.000000000000000e+00

      1.000000000000000e+00 0.000000000000000e+00
      7.423524831263540e-01 -5.635439165128736e-01
      6.860077553104628e-01 7.058976849524824e-01
      9.070639066341516e-01 1.374294020717560e-01
      1.000000000000000e+00 0.000000000000000e+00
      4.017558780880642e-01 -1.046104868497161e-01

      1.000000000000000e+00 0.000000000000000e+00
      8.350017415760391e-01 -3.035440128760801e-01
      3.351517066290223e-01 8.286065487115769e-01
      5.313708156187256e-01 5.901546173030379e-01
      1.000000000000000e+00 0.000000000000000e+00
      2.737797565592244e-01 -5.837061663165312e-02

      1.000000000000000e+00 0.000000000000000e+00
      -2.556076727434309e-02 -5.863529384777062e-01
      -2.526213176462230e-01 3.482485287048019e-01
      2.106537428354917e-01 1.392237523281170e-01
      1.000000000000000e+00 0.000000000000000e+00
      -7.992284323716657e-01 -4.053661776456201e-01




startEDDegree(NCO,(0,0,1),stageOne)
runBertiniStartEDDegree(#F,NCO)--6 residuals and 4 degenerates
readFile(storeBM2Files,"member_points",10000)
o18 = 10                                               

      1.000000000000000e+00 0.000000000000000e+00
      -4.443599600198331e-01 -4.798357879691401e-01
      -9.735263078992706e-02 -7.381498972988272e-02
      7.840537364664571e-03 7.951370219118840e-02
      1.770367044586341e-01 -2.295279826172113e-01
      1.000000000000000e+00 0.000000000000000e+00

      -7.308401557577290e-01 3.088500185581869e-01
      1.000000000000000e+00 0.000000000000000e+00
      3.614504800180809e-01 -2.747932099566459e-01
      -5.544459515121313e-01 1.416897617444813e-01
      -1.001649892761790e-01 -1.942439837792194e-01
      1.000000000000000e+00 0.000000000000000e+00

      -1.221245327087672e-15 -9.999999999999994e-01
      -9.999999999999987e-01 8.326672684688674e-16
      1.000000000000000e+00 0.000000000000000e+00
      2.997602166487923e-15 -9.999999999999982e-01
      1.000000000000000e+00 0.000000000000000e+00
      -9.999999999999984e-01 1.665334536937735e-16

      2.835724772392306e-01 3.577730674559088e-01
      -1.282521953934703e-01 -6.969204071628685e-03
      1.000000000000000e+00 0.000000000000000e+00
      -1.864654733776547e-01 2.106802497017793e-01
      1.000000000000000e+00 0.000000000000000e+00
      -9.947413110043575e-01 -4.157598393624301e-02

      5.544459515121326e-01 -1.416897617444829e-01
      3.614504800180817e-01 -2.747932099566473e-01
      1.000000000000000e+00 0.000000000000000e+00
      7.308401557577300e-01 -3.088500185581867e-01
      1.000000000000000e+00 0.000000000000000e+00
      -1.001649892761794e-01 -1.942439837792202e-01

      1.000000000000000e+00 0.000000000000000e+00
      -3.080868893334809e-15 9.999999999999982e-01
      1.221245327087672e-15 -9.999999999999982e-01
      9.999999999999956e-01 -2.581268532253489e-15
      1.000000000000000e+00 0.000000000000000e+00
      -9.999999999999992e-01 8.326672684688674e-16

      4.267535494419798e-01 8.307453247541303e-01
      1.000000000000000e+00 0.000000000000000e+00
      -3.862995205918740e-01 3.312296226509898e-01
      1.264682845279503e-01 5.299702529608813e-01
      -9.952688224517491e-01 -3.930198333806562e-02
      1.000000000000000e+00 0.000000000000000e+00

      1.387778780781446e-15 -9.999999999999978e-01
      1.000000000000000e+00 0.000000000000000e+00
      9.999999999999997e-01 2.164934898019055e-15
      -7.216449660063518e-16 9.999999999999994e-01
      9.999999999999984e-01 -1.110223024625157e-16
      1.000000000000000e+00 0.000000000000000e+00

      4.853903375103449e-01 -6.099742192880407e-01
      4.748412972865161e-01 -2.808974300402036e-01
      1.000000000000000e+00 0.000000000000000e+00
      6.612503970713429e-01 2.522676188721750e-01
      1.000000000000000e+00 0.000000000000000e+00
      1.596751771451175e-01 2.017837893940159e-01

      -9.999999999999993e-01 -3.885780586188048e-16
      7.771561172376096e-16 9.999999999999996e-01
      1.054711873393899e-15 9.999999999999984e-01
      1.000000000000000e+00 0.000000000000000e+00
      9.999999999999996e-01 2.164934898019055e-15
      1.000000000000000e+00 0.000000000000000e+00




stageOne=1
stageTwo=2

--NCO=newNumericalComputationOptions(theDir,(F,G,L))
--NCO#"TargetWeight"=apply(#gens ring first F,i->1)
NCO#"Directory"=theDir
startEDDegree(NCO,(0,0,null),stageOne)
runBertiniStartEDDegree(#F,NCO)
readFile(theDir,"member_points",10000)

     10                                               

       1.388577090353700e-01 7.353765432279112e-01
       1.000000000000000e+00 0.000000000000000e+00
       -2.038897708738092e-01 -1.840956672770518e-01
       -2.922745165763581e-01 2.220701796875007e-01
       1.000000000000000e+00 0.000000000000000e+00
       -1.444744393010470e-01 -4.406337290250787e-01

       2.016669117159307e-01 -2.675398781345260e-01
       7.967009180643085e-01 9.797873118166754e-02
       2.086734127503752e-01 -3.614724519796661e-01
       1.000000000000000e+00 0.000000000000000e+00
       1.000000000000000e+00 0.000000000000000e+00
       1.897727047356341e-01 -2.985357526268447e-01

       1.000000000000000e+00 0.000000000000000e+00
       -1.470704686225854e-01 9.051415573284101e-01
       -3.920016039633773e-01 7.243506115210331e-01
o25 = 10                                               

      1.000000000000000e+00 0.000000000000000e+00
      8.350017415760396e-01 -3.035440128760798e-01
      3.351517066290213e-01 8.286065487115772e-01
      5.313708156187256e-01 5.901546173030381e-01
      1.000000000000000e+00 0.000000000000000e+00
      2.737797565592243e-01 -5.837061663165330e-02

      1.000000000000000e+00 0.000000000000000e+00
      -2.556076727434309e-02 -5.863529384777060e-01
      -2.526213176462230e-01 3.482485287048018e-01
      2.106537428354914e-01 1.392237523281171e-01
      1.000000000000000e+00 0.000000000000000e+00
      -7.992284323716661e-01 -4.053661776456202e-01

      5.023820446423849e-01 4.548167910299097e-01
      -2.769246618401499e-01 1.103170119655542e-02
      1.000000000000000e+00 0.000000000000000e+00
      -2.920103070944935e-01 2.863217217871003e-01
      1.000000000000000e+00 0.000000000000000e+00
      -6.555591290157397e-01 -6.693750148030579e-01

      -6.666891893131216e-01 5.074842990005829e-01
      1.000000000000000e+00 0.000000000000000e+00
      1.651006729939872e-01 -4.238158851386270e-01
      -4.631679365721976e-01 2.831385787464482e-01
      -2.796203070225711e-01 2.400879156300439e-01
      1.000000000000000e+00 0.000000000000000e+00

      1.000000000000000e+00 0.000000000000000e+00
      7.423524831263542e-01 -5.635439165128751e-01
      6.860077553104640e-01 7.058976849524824e-01
      9.070639066341517e-01 1.374294020717549e-01
      1.000000000000000e+00 0.000000000000000e+00
      4.017558780880642e-01 -1.046104868497162e-01

      1.000000000000000e+00 0.000000000000000e+00
      -3.495186762173847e-01 -4.238479836156784e-01
      -2.195951495118888e-01 -8.881869498496529e-02
      3.910698118440092e-02 1.241187540268984e-01
      -5.570299013891631e-01 -3.123798700592146e-01
      1.000000000000000e+00 0.000000000000000e+00

      1.000000000000000e+00 0.000000000000000e+00
      7.136091195660893e-02 4.911638365853598e-01
      -3.076656085394331e-01 6.399390901710904e-01
      -3.362702371124361e-01 -1.054475836043058e-01
      1.000000000000000e+00 0.000000000000000e+00
      9.889773344040326e-02 6.186440695556038e-01

      -6.289977362595908e-01 1.912844762422850e-01
      3.229330587679998e-01 6.032119326084745e-01
      -1.874164479402353e-01 9.424130040076557e-01
      1.000000000000000e+00 0.000000000000000e+00
      1.000000000000000e+00 0.000000000000000e+00
      1.649271782698877e-01 5.456541500970947e-01

      -8.727865247599269e-01 -4.238955238272206e-01
      2.493561439182778e-01 5.112914443215999e-01
      1.000000000000000e+00 0.000000000000000e+00
      -4.613858711962508e-01 -3.617288189127671e-01
      -9.167465833888921e-01 -6.267791280300139e-03
      1.000000000000000e+00 0.000000000000000e+00

      1.000000000000000e+00 0.000000000000000e+00
      1.686452067024735e-02 8.243968675810909e-01
      1.809342382215980e-01 -8.298517551782532e-01
      6.871785567260658e-01 1.351665671496123e-01
      1.000000000000000e+00 0.000000000000000e+00
      -3.983540147857424e-01 -7.286234417087362e-01


filterSolutionFile("start",#gens ring first F+1+#G+#L,NCO)     
readFile(theDir,"start",100000)

oo28 = 6

       1.000000000000000e+00 0.000000000000000e+00
       8.350017415760396e-01 -3.035440128760798e-01
       3.351517066290213e-01 8.286065487115772e-01
       5.313708156187256e-01 5.901546173030381e-01
       1.000000000000000e+00 0.000000000000000e+00
       2.737797565592243e-01 -5.837061663165330e-02

       1.000000000000000e+00 0.000000000000000e+00
       7.423524831263542e-01 -5.635439165128751e-01
       6.860077553104640e-01 7.058976849524824e-01
       9.070639066341517e-01 1.374294020717549e-01
       1.000000000000000e+00 0.000000000000000e+00
       4.017558780880642e-01 -1.046104868497162e-01

       1.000000000000000e+00 0.000000000000000e+00
       -3.495186762173847e-01 -4.238479836156784e-01
       -2.195951495118888e-01 -8.881869498496529e-02
       3.910698118440092e-02 1.241187540268984e-01
       -5.570299013891631e-01 -3.123798700592146e-01
       1.000000000000000e+00 0.000000000000000e+00

       -6.289977362595908e-01 1.912844762422850e-01
       3.229330587679998e-01 6.032119326084745e-01
       -1.874164479402353e-01 9.424130040076557e-01
       1.000000000000000e+00 0.000000000000000e+00
       1.000000000000000e+00 0.000000000000000e+00
       1.649271782698877e-01 5.456541500970947e-01

       -8.727865247599269e-01 -4.238955238272206e-01
       2.493561439182778e-01 5.112914443215999e-01
       1.000000000000000e+00 0.000000000000000e+00
       -4.613858711962508e-01 -3.617288189127671e-01
       -9.167465833888921e-01 -6.267791280300139e-03
       1.000000000000000e+00 0.000000000000000e+00

       1.000000000000000e+00 0.000000000000000e+00
       1.686452067024735e-02 8.243968675810909e-01
       1.809342382215980e-01 -8.298517551782532e-01
       6.871785567260658e-01 1.351665671496123e-01
       1.000000000000000e+00 0.000000000000000e+00
       -3.983540147857424e-01 -7.286234417087362e-01

startEDDegree(NCO,(0,0,null),stageTwo)
runBertiniStartEDDegree(#F,NCO)--2
readFile(theDir,"member_points",100000)
oo32 = 6                                                

       4.853903375103445e-01 -6.099742192880406e-01
       4.748412972865166e-01 -2.808974300402035e-01
       1.000000000000000e+00 0.000000000000000e+00
       6.612503970713429e-01 2.522676188721749e-01
       1.000000000000000e+00 0.000000000000000e+00
       1.596751771451180e-01 2.017837893940159e-01

       5.544459515121329e-01 -1.416897617444826e-01
       3.614504800180816e-01 -2.747932099566470e-01
       1.000000000000000e+00 0.000000000000000e+00
       7.308401557577303e-01 -3.088500185581864e-01
       1.000000000000000e+00 0.000000000000000e+00
       -1.001649892761792e-01 -1.942439837792199e-01

       1.000000000000000e+00 0.000000000000000e+00
       -4.443599600198334e-01 -4.798357879691401e-01
       -9.735263078992658e-02 -7.381498972988323e-02
       7.840537364664869e-03 7.951370219118850e-02
       1.770367044586340e-01 -2.295279826172110e-01
       1.000000000000000e+00 0.000000000000000e+00

       1.000000000000000e+00 0.000000000000000e+00
       1.318389841742373e-16 -9.999999999999999e-01
       2.151057110211241e-16 -9.999999999999999e-01
       -9.999999999999999e-01 -3.538835890992686e-16
       1.000000000000000e+00 0.000000000000000e+00
       1.000000000000000e+00 0.000000000000000e+00

       -7.308401557577299e-01 3.088500185581869e-01
       1.000000000000000e+00 0.000000000000000e+00
       3.614504800180812e-01 -2.747932099566471e-01
       -5.544459515121326e-01 1.416897617444827e-01
       -1.001649892761796e-01 -1.942439837792201e-01
       1.000000000000000e+00 0.000000000000000e+00

filterSolutionFile("start",#gens ring first F+1+#G+#L,NCO)     
readFile(theDir,"start",10000)


oo52 = 2

       1.000000000000000e+00 0.000000000000000e+00
       -4.443599600198332e-01 -4.798357879691401e-01
       -9.735263078992645e-02 -7.381498972988332e-02
       7.840537364664765e-03 7.951370219118849e-02
       1.770367044586340e-01 -2.295279826172110e-01
       1.000000000000000e+00 0.000000000000000e+00

       1.387778780781446e-17 9.999999999999994e-01
       -9.999999999999994e-01 -1.249000902703301e-16
       1.000000000000000e+00 0.000000000000000e+00
       -1.526556658859590e-16 9.999999999999998e-01
       1.000000000000000e+00 0.000000000000000e+00
       -9.999999999999998e-01 -1.110223024625157e-16



NCO=newNumericalComputationOptions(storeBM2Files,(F,G,L))
NCO#"TargetWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,0,null),stageOne)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)--6
filterSolutionFile(storeBM2Files,"start",#gens ring first F+1+#G+#L,NCO)     
startEDDegree(NCO,(0,0,null),stageTwo)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)--FAILURE


NCO=newNumericalComputationOptions(storeBM2Files,(F,G,L))
NCO#"TargetWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,null,1),stageOne)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)--2






--STAGE TWO TESTS

NCO=newNumericalComputationOptions(storeBM2Files,(F,G,L))
NCO#"TargetWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,0,0),stage)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)--6

NCO=newNumericalComputationOptions(storeBM2Files,(F,G,L))
NCO#"TargetWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,null,0),stage)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)--6

NCO=newNumericalComputationOptions(storeBM2Files,(F,G,L))
NCO#"TargetWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,0,null),stage)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)--6


NCO=newNumericalComputationOptions(storeBM2Files,(F,G,L))
NCO#"TargetWeight"=apply(#gens ring first F,i->1)
startEDDegree(NCO,(0,null,1),stage)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)--2

--HOMOTOPY

--FILTER
filterSolutionFile(storeBM2Files,"start",#gens ring first F+1+#G+#L,NCO)     
readFile(storeBM2Files,"start",1000)
-- STAGE TWO
stage=1
startEDDegree(NCO,(0,0,1),stage)
runBertiniStartEDDegree(storeBM2Files,#F,NCO)
readFile(storeBM2Files,"inputCriticalPointSuperSet",10000)
*-

(stageOne,stageTwo)=(1,2);
numericWeightEDDegree=method()
numericWeightEDDegree(List,List,List,List):=(F,G,L,WV)->(
    NCO:=newNumericalComputationOptions(storeBM2Files,(F,G,L));
    NCO#"StartWeight"=WV;
    homotopyType:=(0,0,0);
    startEDDegree(NCO,homotopyType,stageOne);
    return runBertiniStartEDDegree(#F,NCO)    
    )

numericWeightEDDegree(List,List,List):=(F,G,WV)->numericWeightEDDegree(F,G,{},WV)
numericWeightEDDegree(List,List):=(F,G)->(
    WV:=apply(#gens ring first F,i->random CC);
    numericWeightEDDegree(F,G,WV))

numericUnitEDDegree=method()
numericUnitEDDegree(List,List,List):=(F,G,L)->(
    WV:=apply(#gens ring first F,i->1);
    numericWeightEDDegree(F,G,L,WV))
numericUnitEDDegree(List,List):=(F,G)->numericUnitEDDegree(F,G,{})




    
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
