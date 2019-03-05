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
    NCO#"TargetWeight"=NCO#"StartWeight"=apply(numX,i->random CC); 
    NCO#"GammaVector"=apply(numX-1,i->random CC); 
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
    NCO#"GeneralHyperplaneList"=null;
    return NCO
    )

defaultFixValues={"StartData","StartWeight","StartSubmodel","JacobianStartSubmodel"}


defaultMTName="input_MT_"

stageOne=1
stageTwo=2
stageWeightEDDegreeHomotopy=method()
stageWeightEDDegreeHomotopy(NumericalComputationOptions,ZZ):=(NCO,stage)->(    
--stage is either 0, 1 or 2.
--Now we extract information from NCO.
    theDir:=NCO#"Directory"; 
--    data:=NCO#(NCO#"FixedData");
    jacL:=NCO#(NCO#"FixedJacobianSubmodel");
    L:=NCO#(NCO#"FixedSubmodel");
--Other strings that are used
    (lagMult,numerHB,denomQ,primal,tWeight):=("lagMult","numerHB","denomQ","primal","tWeight")   ; 
    startWeight:=NCO#"StartWeight";
    targetWeight:=NCO#"TargetWeight";
    startData:=NCO#"StartData";
    targetData:=NCO#"StartData";
--    weight:=apply(startWeight,targetWeight,(i,j)->makeB'Section({i,j}))
    (jv,gv,sv):=(NCO#"JacobianVars",    NCO#"GradientVars",    NCO#"ScaleVars");   
    (dv,wv):=(NCO#"DataVars",    NCO#"WeightVars");   
--First the model and submodel constraints
------F is the model, V(G)\cap V(L) is a complete intersection contained in V(F)\cap V(L).
    (F,G,startL,targetL,jacG):=(NCO#"Model",NCO#"WitnessModel",NCO#"StartSubmodel",NCO#"TargetSubmodel",NCO#"JacobianWitnessModel");
    randomGamma:=NCO#"GammaVector";
    xList:=NCO#"PrimalCoordinates";
    nc:=#xList;
    kk0:=QQ;
    vRing:=(numVars,s,kk)->kk[apply(numVars,i->s|i)];   
    extrinsicRing:=kk0[flatten transpose apply(#G+#L,i->apply(nc,j->jv|i|"v"|j))];
    scan({sv,gv,dv,wv},{#G+#L+1,nc,nc,nc},(s,numVars)->extrinsicRing=extrinsicRing**vRing(numVars,s,kk0))
--    **kk0[flatten apply(#G+#L+1,i->sv|i)]**kk0[flatten apply(nc,j->gv|j)]**kk[exVars(nc,dv)];
    symbolicJac:=genericMatrix(extrinsicRing,#G+#L,nc);
    symbolicScaleMatrix:=basis({0,1,0,0,0},extrinsicRing);
    symbolicGradient:=basis({0,0,1,0,0},extrinsicRing);
    symbolicSystem:=symbolicScaleMatrix*(symbolicGradient||symbolicJac);
--    
    data:=symbolicData:=gens vRing(nc,dv,kk0);
    dataHomotopy:=false
    if dataHomotopy 
    then pairData:=apply(nc,i->makeB'Section({startData_i,targetData_i},
    	    B'NumberCoefficients=>{"(1-TData)","TData"},		
	    NameB'Section=>symbolicData_i )) else
    if not dataHomotopy then  pairData=apply(nc,i->makeB'Section({startData_i},
    	    B'NumberCoefficients=>{1},		
	    NameB'Section=>symbolicData_i
	    ));
    peek first pairData
--
    weight:=symbolicWeight:=gens vRing(nc,wv,kk0);
--This should be a function and not repeaded twice.
    weightHomotopy:=true
    if weightHomotopy 
    then pairWeight:=apply(nc,i->makeB'Section({startWeight_i,targetWeight_i},
    	    B'NumberCoefficients=>{"(1-TWeight)","TWeight"},		
	    NameB'Section=>symbolicWeight_i )) else
    if not weightHomotopy then  pairWeight=apply(nc,i->makeB'Section({startWeight_i},
    	    B'NumberCoefficients=>{1},		
	    NameB'Section=>symbolicWeight_i
	    ));        
--Gradient Vector
    print 2;
    kk2:=ring first startWeight;
    topS:=kk2[numerHB,denomQ,tWeight];
    (topNumerHB,topDenomQ,topTWeight):=toSequence flatten entries basis({1},topS);
    pairGradient:=apply(#xList,i->makeB'Section({xList_i},
	    ContainsPoint=>{data_i},
	    B'NumberCoefficients=>{weight_i},
	    NameB'Section=>symbolicGradient_(0,i),
	    B'Homogenization=>"HX"));
--Jacobian matrix
    jacLG:=jacL||jacG;
    kk3:=coefficientRing ring first F;
    jacRing:=kk3[gens ring first F|{"HX"}];
    jac=    matrix apply(numrows jacLG,i->apply(numcols jacLG,j->sub(jacLG_(i,j),jacRing)
	    ));
    homogJac:=transpose homogenize(transpose (jac),HX)    ;
    pairJac:=apply(numrows jacLG,i->apply(numcols jacLG,j->symbolicJac_(i,j)=>sub(jacLG_(i,j),jacRing)));
    (degSubmodel,degWitnessModel):=(NCO#"DegreeSubmodel",NCO#"DegreeWitnessModel");
    maxDegree:=(degSubmodel|degWitnessModel|{2})//max;
    degRescale:=({2}|degSubmodel|degWitnessModel)/(i->maxDegree-i);        
--    degRescale={1,2};
    bLagrangeVars:=lagList:=apply(#degRescale,i->"L"|i);
    rescaling:=new MutableList from apply(#degRescale,i->lagList_i);
    generalHyperplaneList:={};
    scan(#degRescale,i->scan(
	    degRescale_i,
	    j->(hg:="*H"|i|"v"|j;
		rescaling#i=(rescaling#i)|hg;
	    	generalHyperplaneList=generalHyperplaneList|{hg})
	    ))
    rescaling#0;
    rescaling#1;
--    NCO#"PairGeneralHyperplaneList"=null
    if NCO#"PairGeneralHyperplaneList"=!=null then
         pairGeneralHyperplanes:=NCO#"PairGeneralHyperplaneList"    else
         pairGeneralHyperplanes=apply(#generalHyperplaneList,i->makeB'Section(xList|{"HX"},NameB'Section=>generalHyperplaneList_i));
    pairScale:=apply(flatten entries symbolicScaleMatrix,rescaling,(i,j)->i=>j);
--    
    bModelVars:=gens ring first F|{"HX"}   ;
    bPoly:=L|G|flatten entries symbolicSystem;
    bConfiguration:={"UseRegeneration"=>1,
	"TrackType"=>0,
	"PrintPathProgress"=>1000}|(NCO#"BertiniStartFiberSolveConfiguration");    
    BF:=pairData|pairWeight|pairJac//flatten|pairGradient|pairScale|pairGeneralHyperplanes;
    writeSolveInputFile:=(stage,nif)->(
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
	    B'Constants=>BC,
	    ParameterGroup=>PG,
    	    B'Functions=>BF
	    ));
    runSolveInputFile:=(stage,nif)->(
	writeSolveInputFile(stage,nif); 
    	if stage==stageTwo then(
    	    writeParameterFile(NCO#"Directory",{0},NameParameterFile=>"start_parameters");
    	    writeParameterFile(NCO#"Directory",{1},NameParameterFile=>"target_parameters"));
	runBertini(NCO#"Directory",NameB'InputFile=>nif);
    	readFile(NCO#"Directory");
    	if stage==stageOne then(	
	    moveB'File(NCO#"Directory","nonsingular_solutions","stageOne_solutions",CopyB'File=>true);
	    moveB'File(NCO#"Directory","nonsingular_solutions","start",CopyB'File=>true);
	    moveB'File(NCO#"Directory","nonsingular_solutions","member_points",CopyB'File=>true));
    	if stage==stageTwo then(
	    moveB'File(NCO#"Directory","nonsingular_solutions","stageTwo_solutions",CopyB'File=>true);
	    --moveB'File(NCO#"Directory","nonsingular_solutions","start",CopyB'File=>true);
	    moveB'File(NCO#"Directory","nonsingular_solutions","member_points",CopyB'File=>true));	    
	    )
    ttOne=1;
    ttThree=3;    
    nameFile:=(stage,case,indexCase,hypersurface,theTrackType)->("input_first_MT_"|case|"_"|indexCase|"_"|theTrackType)
    writeIsMembershipHypersurface:=(stage,case,indexCase,hypersurface,theTrackType)->(
	nif=nameFile(stage,case,indexCase,hypersurface,theTrackType);
    	if stage===stageOne then BC={"TData"=>0,"TWeight"=>0};
    	if stage===stageTwo then BC={"TData"=>1,"TWeight"=>1};
    	if not member(stage,{1,2}) then error"stage is in {1,2}";
    	makeB'InputFile(NCO#"Directory",
    	    NameB'InputFile=>nif,
	    AffVariableGroup=>flatten flatten {bLagrangeVars,bModelVars},
    	    B'Configs=>bConfiguration|{"TrackType"=>theTrackType},
	    B'Polynomials=>{hypersurface},
	    B'Constants=>BC,
--	    ParameterGroup=>PG,
    	    B'Functions=>BF
	    ));
    runIsMembershipHypersurface:=(stage,case,indexCase,hypersurface,theTrackType)->(
	nif=nameFile(stage,case,indexCase,hypersurface,theTrackType);
	runBertini(NCO#"Directory",NameB'InputFile=>nif));
    isMembershipHypersurface:=(stage,case,indexCase,hypersurface)->(
	--Pos dim solve TrackType=>1
	writeIsMembershipHypersurface(stage,case,indexCase,hypersurface,ttOne);
	runIsMembershipHypersurface(stage,case,indexCase,hypersurface,ttOne);
	nif=nameFile(stage,case,indexCase,hypersurface,ttOne);
    	moveB'File(NCO#"Directory","bertini_session.log","bertini_session_"|nif|".log",CopyB'File => false);
    	print nif;
	--MT TrackType=>3
	writeIsMembershipHypersurface(stage,case,indexCase,hypersurface,ttThree);
	runIsMembershipHypersurface(stage,case,indexCase,hypersurface,ttThree);
	nif=nameFile(stage,case,indexCase,hypersurface,ttThree);
    	moveB'File(NCO#"Directory","bertini_session.log","bertini_session_"|nif|".log",CopyB'File => false);
    	print nif;	
       	outIM:=importIncidenceMatrix(NCO#"Directory");
	print outIM;
	return outIM
	)  ;  
    runSolveInputFile(stageOne,"input_first_solve");
    isMembershipHypersurface(stageOne,"HX",0,HX)      
    
---BACK TO COMPUTATIONS
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
    	writeParameterFile(NCO#"Directory",{1},NameParameterFile=>"target_parameters");
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
