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
