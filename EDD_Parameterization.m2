-- makeJac = (system,unknowns) -> ( for i in system list for j in unknowns list diff(j,i) )
rand := randomZZ

parameterizedWeightEDDegree = method()
parameterizedWeightEDDegree(List, List, List) := (F, U, W) -> (
    R := ring first F;
    coef := coefficientRing R;
    numX := #gens R;
    n := #F;

    S := R ** coef[apply(n - numX, i->"L"|i)];
    xList := flatten entries basis({1,0}, S);
    lamList := flatten entries basis({0,1}, S);
    M := sub(matrix{F}, S);

    -- Find a spanning set for ker(jacM)
    jacM := transpose sub(matrix makeJac(apply(F, i->sub(i,S)), xList), S);
    assert(rank jacM == numX);  -- ensure dim X = d
    columnVectors := gens kernel jacM;
    evalColumnVectors := sub(columnVectors, apply(xList, x -> x => random(1, 100)));

    A := matrix for i to n-1 list {};
    count := 0;
    selectColumns := {};
    while #selectColumns < n-numX and count < #columnVectors do (
        B := A | matrix(evalColumnVectors_count);
        if rank B > rank A then (
            A = B;
            selectColumns = selectColumns | {count}
        );
        count = count + 1
    );
    assert(#selectColumns == n - numX);
    unscaledMatrix := columnVectors_selectColumns;
    outMatrix := unscaledMatrix * transpose matrix{lamList};

    gradObjective := 2 * diagonalMatrix(W) * transpose(M - matrix{U});  -- 2w(f - u)
    I := ideal(outMatrix - gradObjective);
    degree I
)

parameterizedUnitEDDegree = method()
parameterizedUnitEDDegree(List) := (F) -> parameterizedWeightEDDegree(
    F,
    apply(#F, i->rand()),
    apply(#F, i->1_(ring first F))
)

parameterizedGenericEDDegree = method()
parameterizedGenericEDDegree(List) := (F) -> parameterizedWeightEDDegree(
    F,
    apply(#F, i->rand()),
    apply(#F, i->rand())
)

end

-*
n = 3
numx = 2
R = QQ[x_1,x_2,y_1..y_(n-numx), u_1..u_n]
M = matrix{{x_1^2 + 1, x_1 * x_2, x_2 - 1}}
assert(n==numcols M)
assert(#support M==numx)

jacM=diff(matrix for i from 1 to numx list {x_i}, M)
assert(numcols M==numcols jacM)--check

columnVectors = gens kernel jacM
jacM*columnVectors

--initialize
A = matrix for i to n-1 list {}
count =0
selectColumns={}
evalColumnVectors = sub(columnVectors,{x_1=>random(1,100),x_2=>random(1,100)})
currentRank=0;
while count+numx<n or currentRank>=n-numx do(    
    B:=A|matrix (evalColumnVectors_count);
    if rank B>currentRank then (
        A=B;
        selectColumns=selectColumns|{count}
    );
    count=count+1
)
unscaledMatrix=columnVectors_selectColumns

outMatrix = columnVectors_selectColumns*transpose matrix{for i to n-numx-1 list y_(i+1)}

gradObjective =  2*((transpose M) - matrix for i from 1 to n list {u_i})
fixU=for i from 1 to n list u_i=>random(1,100)
I = ideal sub(outMatrix-gradObjective,fixU)
{}==decompose minors(numcols unscaledMatrix,unscaledMatrix )

degree I--- I think this is a multiple of the ED degree
decompose I

imageModel = eliminate(support M, ideal(M - matrix{for i from 1 to n list u_i}))
determinantalUnitEDDegree((sub(imageModel,QQ[support imageModel]))_*)
*-