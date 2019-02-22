--restart
--Projective formulation for intersections with linear spaces
rand:=randomValue
--Assume ring is a complex inexact field
--G is a subset of F. 

defaultMTName="input_MT_"
startEDDegree=method(Options=>{
	TargetSlice=>null,
	UseRegeneration=>"1",
	Data=>"Generic",
	Weight=>"Unit"	
	})
startEDDegree(String,List,List,List):= o->(theDir,F,G,L)->(--F is the model, V(G)\cap V(L) is a complete intersection contained in V(F)\cap V(L)
    theR:=ring first F;
    numX:=#gens theR;
--Make data    
    if o#Data=="Generic" 
    then data:=apply(numX,i->random CC) else data=o#Data;
--Make weight
    if o#Weight=="Generic" 
    then weight:=apply(numX,i->random CC) else
    if o#Weight=="Unit"
    then weight=apply(numX,i->1_CC) else weight=o#Weight;
--
    kk:=(coefficientRing theR);
    extraVars:= kk;   
    while class kk===PolynomialRing do (
	kk=(coefficientRing kk);
	extraVars=extraVars|gens kk);
    if class kk=!=ComplexField then "Error: coefficient ring needs to be a ComplexField. ";
    S:=theR**kk[apply(#L+#G+1,i->"L"|i)]**kk[apply(numX,i->"u"|i)]**kk[apply(numX,i->"w"|i)]**kk["numerHB","denomQ"]**kk[apply(numX-1,i->"gam"|i)]**kk["TL"];
    xList:=flatten entries basis(  {1,0,0,0,0,0,0},S);
    lamList:=flatten entries basis({0,1,0,0,0,0,0},S);
    uList:=flatten entries basis(  {0,0,1,0,0,0,0},S);
    wList:=flatten entries basis(  {0,0,0,1,0,0,0},S);
    (numerHB,denomQ):=toSequence flatten entries basis({0,0,0,0,1,0,0},S);
    gamList:=flatten entries basis({0,0,0,0,0,1,0},S);
    tList:=flatten entries basis(  {0,0,0,0,0,0,1},S);
    jac:=sub(matrix makeJac(apply(L|G,i->sub(i,S)),xList),S);
    topRow:=apply(#weight,i->sub(denomQ,S)*uList_i-sub(numerHB,S)*xList_i*wList_i);
    M:=matrix{topRow}||jac;
    degRescale:={3}|apply(L,i->1)|(G/degree/first);
    maxDeg:=(max degRescale);
    --print degRescale;
    degRescale=apply(degRescale,i->maxDeg-i);
    --print degRescale;
    LV:=matrix{apply(lamList,degRescale,(lam,j)->if j==0 then sub(lam,S) else if j>0 then sub(lam,S)*(sub(numerHB,S))^j else print "Error: Homogenized incorrectly.")};
    --print LV; 
    critEq:=flatten entries((LV*sub(M,S)));
    win:=L|G|apply(#critEq-1,i->critEq_i+gamList_i*last critEq);
    randomGamma:=apply(gamList,i->random CC);
    theConstants:=(transpose{uList,data})|(transpose{wList,weight})|(transpose{gamList,randomGamma});
----Input file 
    makeB'InputFile(theDir,NameB'InputFile=>"inputCriticalPointSuperSet",
	HomVariableGroup=>{xList,lamList},
	B'Configs=>{"UseRegeneration"=>1,"TrackType"=>0,"PrintPathProgress"=>1000},
	B'Polynomials=>win,
	B'Functions=>{numerHB=>sum apply(uList,xList,(u,x)->u*x), denomQ=>sum apply(wList,xList,(w,x)->w*x^2)},
	B'Constants=>theConstants
	);
----Input file Parameter Homotopy
    theT:=first tList;
    if o.TargetSlice=!=null 
    then deformL:=(theT)*L+(1-theT)*(o.TargetSlice)
    else deformL=L;
    targetJac:=sub(matrix makeJac(apply(deformL|G,i->sub(i,S)),xList),S);
    targetM:=matrix{topRow}||targetJac;
    targetCritEq:=flatten entries((LV*sub(targetM,S)));
    targetWin:=deformL|G|apply(#targetCritEq-1,i->targetCritEq_i+gamList_i*last targetCritEq);
--
    makeB'InputFile(theDir,NameB'InputFile=>"input_PH_weight",
	HomVariableGroup=>{xList,lamList},
	B'Configs=>{"ParameterHomotopy"=>2,"PrintPathProgress"=>1000},
	B'Polynomials=>targetWin,
	B'Functions=>{numerHB=>sum apply(uList,xList,(u,x)->u*x),denomQ=>sum apply(wList,xList,(w,x)->w*x^2)},
    	ParameterGroup=>{wList,tList},
	B'Constants=>(transpose{uList,data})|(transpose{gamList,randomGamma})
	);
    writeParameterFile(storeBM2Files,weight|{1},NameParameterFile=>"weight_start_parameters");
--
    makeB'InputFile(theDir,NameB'InputFile=>"input_PH_data",
	HomVariableGroup=>{xList,lamList},
	B'Configs=>{"ParameterHomotopy"=>2,"PrintPathProgress"=>1000},
	B'Polynomials=>targetWin,
	B'Functions=>{numerHB=>sum apply(uList,xList,(u,x)->u*x),denomQ=>sum apply(wList,xList,(w,x)->w*x^2)},
    	ParameterGroup=>{wList,tList},
	B'Constants=>(transpose{wList,weight})|(transpose{gamList,randomGamma})
	);
    writeParameterFile(storeBM2Files,data|{1},NameParameterFile=>"data_start_parameters");
--
-----Membership test files function (s,k,bp)=(string to name the file, tracktype, list of polynomials)
    imt:=(s,k,bp)->makeB'InputFile(theDir,NameB'InputFile=>(defaultMTName|s|toString k),
	AffVariableGroup=>flatten{xList,lamList},
	B'Configs=>{"UseRegeneration"=>o#UseRegeneration,"TrackType"=>k,"PrintPathProgress"=>1000},
	B'Polynomials=>bp,
	B'Functions=>{numerHB=>sum apply(uList,xList,(u,x)->u*x),denomQ=>sum apply(wList,xList,(w,x)->w*x^2)},
	B'Constants=>theConstants
	);
--Filter Residuals
    imt("Residual",1,{last critEq});
    imt("Residual",3,{last critEq});
--UBeta
    imt("Degenerate",1,{"numerHB*denomQ*L0"});
    imt("Degenerate",3,{"numerHB*denomQ*L0"});
--Filer component    
    scan(#F,i->(
    	    imt("Hypersurface_"|i|"_",1,{F_i});
    	    imt("Hypersurface_"|i|"_",3,{F_i});
        )))  

runBertiniStartEDDegree=method(Options=>{})
runBertiniStartEDDegree(String,ZZ):= o->(storeBM2Files,n)->(
    runBertini(storeBM2Files,NameB'InputFile=>"inputCriticalPointSuperSet");
    moveB'File(storeBM2Files,"bertini_session.log","bertini_session_CriticalPointSuperSet.log",CopyB'File => true);
    moveB'File(storeBM2Files,"nonsingular_solutions","member_points");
---Run membershipTest
    runMT:=(s)->(
   	runBertini(storeBM2Files,NameB'InputFile=>defaultMTName|s|"1");
    	moveB'File(storeBM2Files,"bertini_session.log","bertini_session_MT_"|s|"1.log",CopyB'File => true);
    	runBertini(storeBM2Files,NameB'InputFile=>defaultMTName|s|"3");
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
    scan(#imResidual,i->if imResidual_i=!={} and imDegenerate_i=={} and memberEveryHypersurface(i) then EDDeg=EDDeg+1);
    moveB'File(storeBM2Files,"bertini_session_CriticalPointSuperSet.log","bertini_session.log");
--    return(imResidual,imDegenerate,imComponent)
    return(EDDeg)
     )


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








