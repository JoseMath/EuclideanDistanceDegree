--restart
needsPackage"Bertini"


--rand:=()->random(CC)
--makeJac:=(F,x)->apply(F,f->apply(x,j->diff(j,f)))
leftKernelWeightEDDegree=method(Options=>{		})
leftKernelWeightEDDegree(String,ZZ,List,List):= o->(theDir,cod,F,weight)->(
    theR:=ring first F;
    numX:=#gens theR;
    data:=apply(numX,i->randCC());
    S:=theR**(coefficientRing theR)[apply(#F+1,i->"L"|i)]**(coefficientRing theR)[apply(numX,i->"u"|i)]**(coefficientRing theR)[apply(numX,i->"w"|i)];
    xList:=flatten entries basis({1,0,0,0},S);
    lamList:=flatten entries basis({0,1,0,0},S);
    uList:=flatten entries basis({0,0,1,0},S);
--    print uList;
    wList:=flatten entries basis({0,0,0,1},S);
--    print wList;
    c:=#lamList-1;
    jac:=sub(matrix makeJac(apply(F,i->sub(i,S)),xList),S);
    topRow:=apply(#weight,i->2*wList_i*(xList_i-uList_i));
    M:=matrix{topRow}||jac;
    critEq:=flatten entries((matrix{lamList}*sub(M,S)));
    restrictLam:=apply(#lamList-1-cod,i->makeB'Section(drop(lamList,1)));
    win:=restrictLam|F|critEq;
    theConstants:=(transpose{uList,data})|(transpose{wList,weight});
    unitQ:=sum apply(xList,i->i^2);
--    sl=radical ideal singularLocus I;
--    win=saturate(win,sl);
    makeB'InputFile(theDir,
	AffVariableGroup=>{xList},
	HomVariableGroup=>{lamList},
	B'Configs=>{"UseRegeneration"=>1,"TrackType"=>0,"PrintPathProgress"=>1000},
	B'Polynomials=>win,
	B'Constants=>theConstants
	)  )  

leftKernelUnitEDDegree=method(Options=>{		})
leftKernelUnitEDDegree(String,ZZ,List):= o->(theDir,cod,F)->leftKernelWeightEDDegree(theDir,    cod,    F,    apply(#gens ring first F,i->1_(ring first F)))

leftKernelGenericEDDegree=method(Options=>{		})
leftKernelGenericEDDegree(String,ZZ,List):= o->(theDir,cod,F)->leftKernelWeightEDDegree(theDir,    cod,    F,    apply(#gens ring first F,i->randCC()))


runBertiniEDDegree=method(Options=>{		})
runBertiniEDDegree(String):= o->(storeBM2Files)->(
    runBertini(storeBM2Files);
--    readFile(storeBM2Files);
    if storeBM2Files_-1===" " then error (storeBM2Files|" cannot end with a whitespace.");
    if storeBM2Files_-1=!="/" then aString:=storeBM2Files|"/";
    numSols:=null;
    scanLines(ell->(numSols=value ell; break),aString|"nonsingular_solutions");
    return numSols)                                                                            

--Projective formulation
writeLeftKernelProjectiveWeightEDDegree=method(Options=>{		})
writeLeftKernelProjectiveWeightEDDegree(String,List,List):= o->(theDir,F,weight)->(
    theR:=ring first F; 
    numX:=#gens theR;
    data:=apply(numX,i->randCC());
    kk:=(coefficientRing theR);
    S:=theR**kk[apply(#F+1,i->"L"|i)]**kk[apply(numX,i->"u"|i)]**kk[apply(numX,i->"w"|i)]**kk["numerHB","denomQ"]**kk[apply(numX-1,i->"gam"|i)];
    xList:=flatten entries basis({1,0,0,0,0,0},S);
    lamList:=flatten entries basis({0,1,0,0,0,0},S);
    uList:=flatten entries basis({0,0,1,0,0,0},S);
    wList:=flatten entries basis({0,0,0,1,0,0},S);
    gamList:=flatten entries basis({0,0,0,0,0,1},S);
    c:=#lamList-1;
    jac:=sub(matrix makeJac(apply(F,i->sub(i,S)),xList),S);
    --print 0;
    topRow:=apply(#weight,i->sub(value("denomQ"),S)*uList_i-sub(value ("numerHB"),S)*xList_i*wList_i);
    M:=matrix{topRow}||jac;
    degRescale:={3}|(F/degree/first);
    maxDeg:=(max degRescale);
    --print degRescale;
    degRescale=apply(degRescale,i->maxDeg-i);
    --print degRescale;
    LV:=matrix{apply(lamList,degRescale,(lam,j)->if j==0 then sub(lam,S) else if j>0 then sub(lam,S)*(sub(value("numerHB"),S))^j else print "Error: Homogenized incorrectly.")};
    --print LV; 
    critEq:=flatten entries((LV*sub(M,S)));
--    restrictLam:=apply(#lamList-1-c,i->makeB'Section(drop(lamList,1)));
--    win:=restrictLam|F|critEq;
    win:=F|apply(#critEq-1,i->critEq_i+gamList_i*last critEq);
    theConstants:=(transpose{uList,data})|(transpose{wList,weight})|(transpose{gamList,apply(gamList,i->randCC())});
--Write input files. 
    makeB'InputFile(theDir,NameB'InputFile=>"inputRegenerationSuperSet",
	HomVariableGroup=>{xList,lamList},
	B'Configs=>{"UseRegeneration"=>1,"TrackType"=>0,"PrintPathProgress"=>1000},
	B'Polynomials=>win,
	B'Functions=>{value("numerHB")=>sum apply(uList,xList,(u,x)->u*x),value("denomQ")=>sum apply(wList,xList,(w,x)->w*x^2)},
	B'Constants=>theConstants
	);
    makeB'InputFile(theDir,NameB'InputFile=>"inputMemberY",
	AffVariableGroup=>flatten{xList,lamList},
	B'Configs=>{"UseRegeneration"=>1,"TrackType"=>1,"PrintPathProgress"=>1000},
	B'Polynomials=>topRow,
	B'Functions=>{value("numerHB")=>sum apply(uList,xList,(u,x)->u*x),value("denomQ")=>sum apply(wList,xList,(w,x)->w*x^2)},
	B'Constants=>theConstants
	);
    imt:=(s,k,bp)->makeB'InputFile(theDir,NameB'InputFile=>("inputMembershipTest"|s|toString k),
	AffVariableGroup=>flatten{xList,lamList},
	B'Configs=>{"UseRegeneration"=>1,"TrackType"=>k,"PrintPathProgress"=>1000},
	B'Polynomials=>bp,
	B'Functions=>{value("numerHB")=>sum apply(uList,xList,(u,x)->u*x),value("denomQ")=>sum apply(wList,xList,(w,x)->w*x^2)},
	B'Constants=>theConstants
	);
    imt("Residual",1,{last critEq});
    imt("Residual",3,{last critEq});
    imt("Degenerate",1,{"numerHB*denomQ*L0"});
    imt("Degenerate",3,{"numerHB*denomQ*L0"})    )  


writeLeftKernelProjectiveUnitEDDegree=method(Options=>{		})
writeLeftKernelProjectiveUnitEDDegree(String,List):= o->(theDir,F)->writeLeftKernelProjectiveWeightEDDegree(theDir,F,apply(#gens ring first F,i->1_(ring first F)))

writeLeftKernelProjectiveGenericEDDegree=method(Options=>{		})
writeLeftKernelProjectiveGenericEDDegree(String,List):= o->(theDir,F)->writeLeftKernelProjectiveWeightEDDegree(theDir,F,apply(#gens ring first F,i->randCC()))

runBertiniProjectiveEDDegree=method(Options=>{		})
runBertiniProjectiveEDDegree(String):= o->(storeBM2Files)->(
    runBertini(storeBM2Files,NameB'InputFile=>"inputRegenerationSuperSet");
    readFile(storeBM2Files); 
    moveB'File(storeBM2Files,"bertini_session.log","bertini_session0.log",CopyB'File => true);
    moveB'File(storeBM2Files,"nonsingular_solutions","member_points");
    runBertini(storeBM2Files,NameB'InputFile=>"inputMembershipTestResidual1");
    runBertini(storeBM2Files,NameB'InputFile=>"inputMembershipTestResidual3");
    imResidual:=importIncidenceMatrix(storeBM2Files);                                                                           
    runBertini(storeBM2Files,NameB'InputFile=>"inputMembershipTestDegenerate1");
    runBertini(storeBM2Files,NameB'InputFile=>"inputMembershipTestDegenerate3");
    imDegenerate:=importIncidenceMatrix(storeBM2Files);
    EDDeg:=0;
    scan(#imResidual,i->if imResidual_i=!={} and imDegenerate_i=={} then EDDeg=EDDeg+1);
    moveB'File(storeBM2Files,"bertini_session0.log","bertini_session.log");
    return EDDeg )                                                                           

end
R=CC[x,y,z]
f=x*y^2-z^3
writeLeftKernelProjectiveGenericEDDegree(theDir,{f})
runBertiniProjectiveEDDegree(theDir)



R=CC[x1,x2,x3,x4,x5,x6,x7,x8,x9]
f=det genericMatrix(R,3,3)
writeLeftKernelProjectiveGenericEDDegree(theDir,{f})--39
runBertiniProjectiveEDDegree(theDir)

writeLeftKernelProjectiveUnitEDDegree(theDir,{f})--3
runBertiniProjectiveEDDegree(theDir)
