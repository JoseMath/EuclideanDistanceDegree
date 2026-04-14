--restart
--Projective formulation for intersections with linear spaces
rand := randomValue
(stageOne,stageTwo) = (1,2); 

--Assume ring is a complex inexact field
--G is a subset of F. 

NumericalComputationOptions = new Type of MutableHashTable

parameterKeys = {"StartWeight", "TargetWeight", "StartData", "TargetData", "GammaVector"}
jacKeys= {"JacobianWitnessModel", "JacobianStartSubmodel", "JacobianTargetSubmodel"}
modelKeys = {"Model","WitnessModel", "StartSubmodel", "TargetSubmodel"}
degreeKeys = {"DegreeWitnessModel", "DegreeSubmodel"}

bertiniKeys = {"BertiniStartFiberSolveConfiguration", "BertiniMembershipTestConfiguration", "BertiniSubstitute", "BertiniConstants"}
coordinateKeys = {"PrimalCoordinates", "HomogeneousVariableGroups", "AffineVariableGroups"}
planeKeys = {"Infinity", "PairGeneralHyperplaneList"}
variableKeys = {"JacobianVars", "GradientVars", "ScaleVars", "DataVars", "WeightVars"}

directoryKeys = {"Directory"} 
solutionKeys = {"TrackSolutions"}
outputKeys = {"OutputType", "FinerRestriction"}
fixValues = {"FixedData", "FixedWeight", "FixedSubmodel", "FixedJacobianSubmodel"}

nocKeys = parameterKeys|jacKeys|modelKeys|degreeKeys|bertiniKeys|coordinateKeys|planeKeys|variableKeys
nocKeys = nocKeys|directoryKeys|solutionKeys|outputKeys|fixValues

newNumericalComputationOptions = method(Options => { })
newNumericalComputationOptions(String, List, List, List) := o -> (dir, F, G, L) -> (
    NCO := new NumericalComputationOptions from apply(nocKeys, i -> i=>null);

    -- Temp directory for Bertini input file
    NCO#"Directory"=dir;

    -- Model keys
    NCO#"Model"=F;
    NCO#"WitnessModel"=G;
    NCO#"StartSubmodel"=L;
    NCO#"TargetSubmodel"=L;

    -- Degrees of models
    NCO#"DegreeSubmodel"=L/degree/first;
    NCO#"DegreeWitnessModel"=G/degree/first;

    -- Jacobians of models
    NCO#"JacobianWitnessModel"=diff(matrix{gens ring first G}, transpose matrix{G});
    NCO#"JacobianTargetSubmodel"=NCO#"JacobianStartSubmodel"=diff(matrix{gens ring first G}, transpose matrix {L});
    
    -- Data keys (used for other homotopy types)
    numX:=#gens ring first G;
    NCO#"TargetData"=NCO#"StartData"=apply(numX,i->random CC); 
    NCO#"TargetWeight"=apply(numX,i->1);
    NCO#"StartWeight"=apply(numX,i->random CC); 
    --NCO#"GammaVector"=apply(numX-1,i->random CC); 

    -- Variable group keys for Bertini
    scan(bertiniKeys,i->NCO#i={});
    NCO#"HomogeneousVariableGroups"={gens ring first F};
    NCO#"AffineVariableGroups"={};  
    NCO#"PrimalCoordinates"=gens ring first F;  -- This is different when working with a parameterization

    -- Fixed keys (which model to use for which run)
    NCO#"FixedData"="StartData";
    NCO#"FixedWeight"="StartWeight";
    NCO#"FixedSubmodel"="StartSubmodel";
    NCO#"FixedJacobianSubmodel"="JacobianStartSubmodel";

    -- Keys defining variable block names
    NCO#"JacobianVars"="jv";
    NCO#"GradientVars"="gv";
    NCO#"ScaleVars"="sv";
    NCO#"DataVars"="dv";
    NCO#"WeightVars"="wv";

    -- Hyperplanes
    NCO#"Infinity"=null;
    NCO#"PairGeneralHyperplaneList"=null;

    -- Flags
    NCO#"IsProjective"=false;  -- currently not used
    NCO#"OutputType"="Standard";
    NCO#"FinerRestriction"={};

    return NCO
)
newNumericalComputationOptions(String, List, List) := o -> (dir, F, G) -> newNumericalComputationOptions(dir, F, G, {})
newNumericalComputationOptions(List, List, List) := o -> (F, G, L) -> newNumericalComputationOptions(temporaryFileName(), F, G, L)
newNumericalComputationOptions(List, List) := o -> (F, G) -> newNumericalComputationOptions(F, G, {})

--##########################################################################--
-- Homotopy ED Degree method
--##########################################################################--

-- homotopyEDDegree
-- ------------------------------------------------------------------------
-- Purpose:
--   Orchestrates a two-stage parameter homotopy for ED degree computations
--   using Bertini, given a highly configurable NumericalComputationOptions (NCO).
--   The homotopy may vary weights ("Weight"), data ("Data"), or (intended)
--   submodel ("Submodel") parameters.
--
-- Inputs:
--   NCO : NumericalComputationOptions (a MutableHashTable) that contains:
--         - Model (F)                   : List of defining polynomials for the target variety
--         - WitnessModel (G)           : List of polynomials forming a CI containing F
--         - Start/TargetSubmodel (L)   : Additional equations (may be {})
--         - Fixed* keys                 : Which entries (Start vs Target) to use for a run
--         - Start/TargetData, Start/TargetWeight : parameter vectors (length = # variables)
--         - DegreeSubmodel, DegreeWitnessModel   : degrees of submodel and witness equations
--         - PrimalCoordinates (xList)            : variable tuple for the ambient space
--         - JacobianVars/GradientVars/ScaleVars/DataVars/WeightVars : symbol names for
--           the corresponding variable blocks in the symbolic ring used to build systems
--         - Infinity: optional Bertini section for homogenization (created if null)
--         - Directory: scratch/output directory for all Bertini IO
--         - (many other keys; see newNumericalComputationOptions)
--
--   ht         : String in {"Weight","Data","Submodel"} selecting the parameter homotopy.
--   isStageOne : Boolean; if true, run stage 1 solve.
--   isStageTwo : Boolean; if true, run stage 2 solve (and then filter/return counts).
--
-- Output:
--   The function returns the number of surviving critical points at the last stage
--   that was executed: stage 1 (index 1) or stage 2 (index 2), taken from
--   stageEDDegBound (a 3-entry MutableList).

possibleHT = {"Weight", "Data", "Submodel"}
stageOne=1
stageTwo=2
ht="Weight"
isStageOne=true
isStageTwo=true

homotopyEDDegree = method()
homotopyEDDegree(NumericalComputationOptions, String, Boolean, Boolean) := (NCO, ht, isStageOne, isStageTwo) -> (    
    ------------------------------------------------------------------------
    -- (CODE 1) Determine which parameter family the homotopy will vary.
    --           Exactly one of weightHomotopy/dataHomotopy/submodelHomotopy is true.
    ------------------------------------------------------------------------
    if not member(ht,possibleHT) then error("Argument 1 is not in "|toString possibleHT);
    weightHomotopy := ht === possibleHT_0;
    dataHomotopy := ht === possibleHT_1;
    submodelHomotopy := ht === possibleHT_2;

    -- Set up directory
    if not fileExists NCO#"Directory" then mkdir NCO#"Directory";

   ------------------------------------------------------------------------
    -- (CODE 2) Extract all required configuration from NCO.
    --           Fixed* keys point to which version ("Start*" or "Target*") to use.
    ------------------------------------------------------------------------
    jacL := NCO#(NCO#"FixedJacobianSubmodel");
    L := NCO#(NCO#"FixedSubmodel");
    startWeight := NCO#"StartWeight";
    targetWeight := NCO#"TargetWeight";
    startData := NCO#"StartData";
    targetData := NCO#"TargetData";
    
    -- Symbol names for special variables used in the homotopy machinery/inputs
    (lagMult,numerHB,denomQ,primal,tWeight) := ("lagMult","numerHB","denomQ","primal","tWeight"); 
    (jv,gv,sv) := (NCO#"JacobianVars",NCO#"GradientVars",NCO#"ScaleVars");   
    (dv,wv) := (NCO#"DataVars",NCO#"WeightVars");
    
    -- F is the "model" (target variety). V(G) ∩ V(L) is a CI inside V(F) ∩ V(L).
    (F,G,startL,targetL,jacG) := (NCO#"Model",NCO#"WitnessModel",NCO#"StartSubmodel",NCO#"TargetSubmodel",NCO#"JacobianWitnessModel");

    -- Optional gamma vector for randomization (not used here; keep for extensibility)
    -- randomGamma:=NCO#"GammaVector";

    -- List of primal coordinates (ambient variables)
    xList := NCO#"PrimalCoordinates";

    -- Ensure a hyperplane at infinity section is set; create it if missing.
    if NCO#"Infinity" === null then NCO#"Infinity" = makeB'Section(xList, NameB'Section=>"HX");

    ------------------------------------------------------------------------
    -- (FUNCTION 0) Helper: create a polynomial ring with numVars variables
    --              using symbol base 's' over coefficient ring 'kk'.
    ------------------------------------------------------------------------
    vRing := (numVars,s,kk) -> kk[apply(numVars,i->s|i)];   

    ------------------------------------------------------------------------
    -- (CODE 3) Build the extrinsic symbolic ring and block variables.
    --          This ring hosts:
    --            - Symbolic Jacobian entries (jv blocks indexed by rows/columns)
    --            - Scale variables (sv) for Lagrange multipliers
    --            - Gradient variables (gv)
    --            - Data (dv) and Weight (wv) variable placeholders
    --          Note: the 'basis({...}, ring)' calls below expect those blocks to be
    --          identifiable by degree; in practice, we are using them to select blocks.
    ------------------------------------------------------------------------
    nc := #xList;
    kk0 := QQ; 
    extrinsicRing := kk0[flatten transpose apply(#G+#L,i->apply(nc,j->jv|i|"v"|j))];
    scan({sv,gv,dv,wv}, {#G+#L+1,nc,nc,nc}, (s,numVars) -> extrinsicRing = extrinsicRing ** vRing(numVars,s,kk0));

    symbolicScaleMatrix := basis({0,1,0,0,0},extrinsicRing);  -- the sv's
    symbolicGradient := basis({0,0,1,0,0},extrinsicRing);  -- the gv's
    symbolicJac := genericMatrix(extrinsicRing,#G+#L,nc);

    -- symbolicSystem is the model system: (Scale) * (Gradient || Jacobian) that we want to solve
    -- to be specialized by pairing functions below.
    symbolicSystem := symbolicScaleMatrix*(symbolicGradient||symbolicJac);

    ------------------------------------------------------------------------
    -- (FUNCTION 1) Pair a row of a matrix with numeric/polynomial values.
    --   M1 : matrix of symbolic variables (e.g., a row of 'symbolicJac')
    --   M2 : matrix of actual polynomials/values to substitute
    --   hom: homogenizing variable name for Bertini (e.g., "HX")
    --   Behavior:
    --     - If val_i == 0, record arg_i => 0 in 'jacZero' and substitute into symbolicSystem.
    --     - Else, create a Bertini section binding NameB'Section=arg_i to val_i with one coefficient.
    ------------------------------------------------------------------------
    jacZero:={};
    pairJac:={};
    pairRowFunction := (M1, M2, hom) -> (
        arg := flatten entries M1;
        val := flatten entries M2;
        scan(#arg, i -> if val_i==0 then (
                jacZero = jacZero|{arg_i=>0};
                symbolicSystem = sub(symbolicSystem, {arg_i=>0})
            )
            else pairJac = pairJac | {makeB'Section(
                {val_i},
                B'NumberCoefficients=>{1},	
                B'Homogenization=>hom,	
                NameB'Section=>arg_i
            )}
        )
    );
    -- Example:
    --   M1 = matrix{{x,y},{z,w}} ; M2 = matrix{{1,2},{5,4}} ; pairRowFunction(M1,M2,"HX")

    ------------------------------------------------------------------------
    -- (FUNCTION 2) Pair a parameter vector for a parameter homotopy.
    --   p0,p1 : start and target parameter vectors (same length)
    --   r1,r2 : strings representing (1-T) and (T) coefficient symbols for Bertini
    --   sym   : the symbolic parameter vector (matches length of p0 and p1)
    --   bool  : whether to do a homotopy (true) or hold fixed at start (false)
    --   Returns:
    --     a list of B'Sections (one per component) binding sym_i to p0_i (and p1_i if homotopy).
    ------------------------------------------------------------------------
    pairParameterFunction := (p0,p1,r1,r2,sym,bool) -> (
	    if bool then apply(#p0,
            --- for each i we encode an equation like sym_i = p_0_i* r1+p1_i*r2,
            --- x = 2*r1+ 7*r2; Typically, r1 = 1-T and r2 = T.
            i->makeB'Section(
                {p0_i,p1_i},
                B'NumberCoefficients=>{r1,r2},		
                NameB'Section=>sym_i 
        ))
        else apply(#p0,
            --- for each i we encode an equation like sym_i = p0_i * 1,
            --- x = 2*1;  
            i->makeB'Section(
                {p0_i},
                B'NumberCoefficients=>{1},		
                NameB'Section=>sym_i 
        ))
	);	

    ------------------------------------------------------------------------
    -- (CODE 4) Pair symbolic Data/Weight with their start/target values.
    --          These define the parameter homotopies controlled by 'ht'.
    ------------------------------------------------------------------------
    weight := symbolicWeight := gens vRing(nc,wv,kk0);
    data := symbolicData := gens vRing(nc,dv,kk0);
    pairData := pairParameterFunction(startData,targetData, "(1-TData)","TData", symbolicData, dataHomotopy);	    
    pairWeight := pairParameterFunction(startWeight,targetWeight, "(1-TWeight)","TWeight", symbolicWeight, weightHomotopy);

    --print ("pairData", peek first pairData);
    --print ("pairWeight", peek first pairWeight);

    ------------------------------------------------------------------------
    -- (CODE 5) Pair the gradient vector:
    --   For each coordinate, define a section that encodes the ED gradient equation
    --   with ContainsPoint => {data_i}, Weight => {weight_i}, and homogenization "HX".
    --   topNumerHB / topDenomQ / topTWeight are extracted from a  parameter ring,
    --   but not used further (kept for possible extensions).
    ------------------------------------------------------------------------
    kk2 := ring first startWeight;
    topS := kk2[numerHB,denomQ,tWeight];
    (topNumerHB,topDenomQ,topTWeight) := toSequence flatten entries basis({1},topS);
    
    pairGradient := apply(#xList, i -> makeB'Section(
        {xList_i},
	    ContainsPoint=>{data_i},
	    B'NumberCoefficients=>{weight_i},
	    NameB'Section=>symbolicGradient_(0,i),
	    B'Homogenization=>"HX"
    ));

    ------------------------------------------------------------------------
    -- (CODE 6) Pair the Jacobian rows:
    --   - Build a homogenized copy of L|G in a ring with HX
    --   - Differentiate each homogenized equation w.r.t. each variable to form homogJac
    --   - Pair rows of symbolicJac with values from homogJac
    ------------------------------------------------------------------------
    jacLG := jacL||jacG;
    kk3 := coefficientRing ring first F;
    jacRing := kk3[gens ring first F|{"HX"}];
    HX := last gens jacRing;

    homogLG := homogenize(sub(matrix{L|G},jacRing), HX) // entries // flatten;
    homogJac := matrix apply(numrows jacLG,i->apply(numcols jacLG,j->diff((gens jacRing)_j,homogLG_i)));
    pairRowFunction(symbolicJac, homogJac, "HX");
    -- print("homogenized jacobian",homogJac);

    ------------------------------------------------------------------------
    -- (CODE 7) Pair scaling variables (Lagrange multipliers) for column homogenization:
    --   - Compute rescaling degrees per column to balance Jacobian columns
    --   - Create a list of new "generic hyperplanes" and pair them as Bertini sections
    --   - Build the diagonal scaling (pairScale) binding each scale slot with its product
    ------------------------------------------------------------------------
    (degSubmodel,degWitnessModel) := (NCO#"DegreeSubmodel",NCO#"DegreeWitnessModel");
    degAugJac := {1} | apply(degSubmodel|degWitnessModel, i -> i-1);
    maxDegree := degAugJac//max;
    degRescale := degAugJac/(i -> maxDegree-i);  -- how much to rescale each column
    bLagrangeVars := lagList := apply(#degRescale, i -> "L"|i);
    rescaling := new MutableList from apply(#degRescale, i -> lagList_i);
    
    -- Build generic linear forms H(i,v,j) used to homogenize per-column degrees
    --Homogenize cols by multiplying by a diagonal matrix of linear products on the left. 
    generalHyperplaneList:={};
    scan(#degRescale, i -> scan(
	    degRescale_i, 
	    j -> (
            hg:="H"|i|"v"|j;  --wants to be both
		    rescaling#i = (rescaling#i)|"*"|hg;
            generalHyperplaneList = generalHyperplaneList|{hg}
        )
    ));
    
    -- Reuse previously paired hyperplanes if provided; otherwise pair new ones now.
    if NCO#"PairGeneralHyperplaneList"=!=null then pairGeneralHyperplanes:=NCO#"PairGeneralHyperplaneList"
    else ( 
        pairGeneralHyperplanes=apply(#generalHyperplaneList, i->
            makeB'Section(xList|{"HX"},
            NameB'Section=>generalHyperplaneList_i)
        );
        NCO#"PairGeneralHyperplaneList"=pairGeneralHyperplanes
    );
    --print(peek first pairGeneralHyperplanes);

    -- Connect each scale slot with the corresponding homogenizing product.
    pairScale := apply(flatten entries symbolicScaleMatrix, rescaling, (i,j) -> i=>j);
    
    ------------------------------------------------------------------------
    -- (CODE 8) Compose Bertini inputs:
    --   - Variable groups, polynomials, constants, configuration (including regeneration)
    --   - BF is the full list of paired functions/constants used across writes
    ------------------------------------------------------------------------
    bModelVars := gens ring first F|{"HX"};
    bPoly := homogLG|flatten entries symbolicSystem;
    bConfiguration := {"TrackType"=>0, "PrintPathProgress"=>1000} | (NCO#"BertiniStartFiberSolveConfiguration");    
    BF := pairData|pairWeight|pairJac|pairGradient|pairGeneralHyperplanes|pairScale;

    ------------------------------------------------------------------------
    -- (FUNCTIONS 2) Solve helpers: write Bertini input and run it for a stage.
    --   Stage behavior:
    --     - Stage 1: ParameterHomotopy=1, UseRegeneration=1, parameters fixed at 0
    --     - Stage 2: ParameterHomotopy=2, UseRegeneration=0, with parameter path files
    --   Note: PG is a placeholder in stage 1 (no active parameter); BC sets TData/TWeight.
    ------------------------------------------------------------------------
    writeSolveInputFunction := (stage,nif) -> (
        if stage===stageOne then (
            PG:={"adfadfdf"};
            BC:={"TData"=>0,"TWeight"=>0}
        )
        else if stage===stageTwo then (
            BC={};
            if dataHomotopy then PG={"TData"}
            else if weightHomotopy then PG={"TWeight"}
        );

        if stage===stageOne then bConfiguration = bConfiguration | {"UseRegeneration"=>1};
        if stage===stageTwo then bConfiguration = bConfiguration | {"UseRegeneration"=>0};

        makeB'InputFile(
            NCO#"Directory",
            NameB'InputFile=>nif,
            HomVariableGroup=>{bLagrangeVars,bModelVars},
            BertiniInputConfiguration=>bConfiguration|{"ParameterHomotopy"=>stage},
            B'Polynomials=>bPoly,
            B'Constants=>jacZero|BC,
            ParameterGroup=>PG,
            B'Functions=>BF
        )
    );

    --our solution file will always be member points. 
    criticalPointName := "criticalPointFile";
    runSolveInputFunction := (stage,nif) -> (
        writeSolveInputFunction(stage,nif); 

        if stage==stageTwo then(
            writeParameterFile(NCO#"Directory",{0},NameParameterFile=>"start_parameters");
            writeParameterFile(NCO#"Directory",{1},NameParameterFile=>"final_parameters")
        );

        runBertini(NCO#"Directory",NameB'InputFile=>nif);
        readFile(NCO#"Directory");

        if stage==stageOne then(	
            moveB'File(NCO#"Directory","bertini_session.log","stageOne_log",CopyB'File=>true);
            moveB'File(NCO#"Directory","nonsingular_solutions","stageOne_solutions",CopyB'File=>true);
            moveB'File(NCO#"Directory","nonsingular_solutions","start",CopyB'File=>true);
            moveB'File(NCO#"Directory","nonsingular_solutions","member_points",CopyB'File=>true);
            moveB'File(NCO#"Directory","nonsingular_solutions",criticalPointName,CopyB'File=>true)
        );

        if stage==stageTwo then(
            moveB'File(NCO#"Directory","bertini_session.log","stageTwo_log",CopyB'File=>true);
            moveB'File(NCO#"Directory","nonsingular_solutions","stageTwo_solutions",CopyB'File=>true);
            moveB'File(NCO#"Directory","main_data","stageTwo_main_data",CopyB'File=>true);
            --moveB'File(NCO#"Directory","nonsingular_solutions","start",CopyB'File=>true);
            moveB'File(NCO#"Directory","nonsingular_solutions","member_points",CopyB'File=>true);
            moveB'File(NCO#"Directory","nonsingular_solutions",criticalPointName,CopyB'File=>true)
        );	    
    );

    ------------------------------------------------------------------------
    -- (FUNCTIONS 3) Membership tests and incidence matrices.
    --   We create and run two inputs per hypersurface:
    --     - TrackType=1 (positive-dimensional solve)
    --     - TrackType=3 (membership test to decide if points lie on the hypersurface)
    --   Output: an incidence matrix (list of lists) from importIncidenceMatrix.
    ------------------------------------------------------------------------    
    ttOne:=1;
    ttThree:=3;    
    nameFileFunction:=(stage,case,indexCase,hypersurface,theTrackType)->("input_first_MT_"|case|"_"|indexCase|"_"|theTrackType);

    writeIsMembershipFunction := (stage,case,indexCase,hypersurface,theTrackType) -> (
	    nif := nameFileFunction(stage,case,indexCase,hypersurface,theTrackType);
    	if stage===stageOne then BC:={"TData"=>0,"TWeight"=>0};
    	if stage===stageTwo then BC={"TData"=>1,"TWeight"=>1};
    	if not member(stage,{1,2}) then error"stage is in {1,2}";

    	makeB'InputFile(
            NCO#"Directory",
    	    NameB'InputFile=>nif,
	        AffVariableGroup=>flatten flatten {bLagrangeVars,bModelVars},
    	    BertiniInputConfiguration=>bConfiguration|{"TrackType"=>theTrackType},
	        B'Polynomials=>{hypersurface},
	        B'Constants=>jacZero|BC,
	        --ParameterGroup=>PG,
    	    B'Functions=>BF
	    )
    );

    -- Example test: isMembershipFunction(stageOne, "TT", 0, "x1*x2-x3*x4")
    isMembershipFunction:=(stage,case,indexCase,hypersurface)->(
        --Pos dim solve TrackType=>1
        writeIsMembershipFunction(stage,case,indexCase,hypersurface,ttOne);
        nif:=nameFileFunction(stage,case,indexCase,hypersurface,ttOne);
        runBertini(NCO#"Directory",NameB'InputFile=>nif);
        moveB'File(NCO#"Directory","bertini_session.log","bertini_session_"|nif|".log",CopyB'File => false);

        --MT TrackType=>3
        writeIsMembershipFunction(stage,case,indexCase,hypersurface,ttThree);
        nif=nameFileFunction(stage,case,indexCase,hypersurface,ttThree);
        runBertini(NCO#"Directory",NameB'InputFile=>nif);
        moveB'File(NCO#"Directory","bertini_session.log","bertini_session_"|nif|".log",CopyB'File => false);

        outIM := importIncidenceMatrix(NCO#"Directory");
        --print("Membership tests",nif);	
        --print outIM;
        return outIM
    );
    
    ------------------------------------------------------------------------
    -- (FUNCTIONS 4) Filtering using incidence matrix output.
    --   - filterSolutionFunction: rewrite "member_points" to only include selected
    --     solutions (by indices kp), preserving Bertini's line grouping.
    --   - positionFilterFunction: run membership test on a single polynomial, decide
    --     which solutions to KEEP (typeB) or DROP (typeA), rewrite "member_points",
    --     and return the number of retained solutions.
    ------------------------------------------------------------------------
    filterSolutionFunction:=(nsf,kp,ns,numCoords)->(
        -- member_points structure assumption: first line header, then groups
        -- of (1 + numCoords) lines per solution. Only selected solutions are kept.
	    
        -- print("RUN FILTER",kp=>numCoords);    
    	firstLine := true;
    	countSol := 0;
    	countLine := 0;
    	groupSize := 1+numCoords;
    	isSelected := null;

    	sf := openOut(NCO#"Directory"|"/"|nsf);
    	scanLineSolutionFunction := (ell) -> (
      	    if firstLine then (
                firstLine=false;
                sf << toString(#kp) << endl
            )
      	    else if countSol < ns then (
                if countLine==0 then isSelected = member(countSol,kp);
                countLine = countLine+1;
                if isSelected then sf << ell << endl;
                if countLine == groupSize then (
                    --print(countSol => isSelected);    	    	
                    --print (countLine,groupSize,"grp");
                    countLine = 0; 
                    countSol = countSol+1;
                )
            )
        );
        scanLines(scanLineSolutionFunction,(NCO#"Directory")|"/"|"member_points");      
        close sf;
        return (nsf)
    );

    -- positionFilterFunction:
    --   stage   : stage1 or stage2 (affects BC and file naming)
    --   case    : label (e.g. "SaturateH", "IntersectF")
    --   indexCase : index within the polyList
    --   hypersurface : polynomial to test
    --   bin     : "typeA" (drop if ON hypersurface) or "typeB" (keep only if ON)

    --filterSolutionFunction("T1",{1,2,3,4,5,6,7},8)
    --	saturateFunction=positionFunction=positionFilterFunction;
    positionFilterFunction := (stage,case,indexCase,hypersurface,bin) -> (
        --(stage,case,indexCase,hypersurface)
        --isMembershipFunction(stage,case,indexCase,hypersurface);
        --(kp,ns):=positionMembershipFunction(stage,case,indexCase,hypersurface);
    	if bin==="typeA" then isOffHypersurface := (m->(m==={}))
    	else if bin==="typeB" then isOffHypersurface = (m->(m=!={}))
	    else error"last argument is typeA or typeB";

	    imMT := isMembershipFunction(stage,case,indexCase,hypersurface);
    	kp := {};
    	scan(#imMT, i -> if isOffHypersurface(imMT_i) then kp=kp|{i});

    	ns:= #imMT;
	    (nsf,nc) := ("filterFile", #flatten {bLagrangeVars,bModelVars});
    	--print("Filter",kp,"num kp",#kp,"num sols",ns,"num coordinates",nc,bin);

	    filterSolutionFunction("filterFile",kp,ns,nc);
    	moveB'File(NCO#"Directory","filterFile","member_points",CopyB'File=>true);
	    return #kp
	);

    -- stageEDDegBound holds a label and the counts after stage1 and stage2 filtering.
    stageEDDegBound:=new MutableList from {"GEDvUED",null,null};
    
    ------------------------------------------------------------------------
    -- (FUNCTIONS 5)  filtering helpers.
    --   runSaturateUnionFunction(polyList, stage):
    --     For each polynomial in 'offPolyList', DROP solutions that lie on it (typeA).
    -- (FUNCTIONS 6) Functions to restrict to the variety 
    --   runRestrictIntersectionFunction(polyList, stage):
    --     For each polynomial in 'onPolyList', KEEP only solutions that lie on it (typeB).
    ------------------------------------------------------------------------
    runSaturateUnionFunction:=(polyList,stage)->(
    	--print("Remove critical points from member_points where any of these polynomials vanish",polyList);
    	(case,bin) := ("SaturateH","typeA");
    	scan(#polyList,i->(
            stageEDDegBound#stage = positionFilterFunction(stageOne,case,i,polyList_i,bin);
            --print(peek stageEDDegBound,"Saturate by polynomial"=>polyList_i)
        ));	 
    	--print(peek stageEDDegBound)
    );

    runRestrictIntersectionFunction:=(polyList,stage)->(
    	--print("Only keeping critical points from member_points where every one of these polynomials vanish",polyList);
    	(case, bin) := ("IntersectF", "typeB");
    	scan(#polyList, i -> (
		    stageEDDegBound#stage = positionFilterFunction(stageOne, case, i, polyList_i, bin);
		    --print(peek stageEDDegBound,"Vanish polynomial"=>polyList_i)
        ));	 
    	--print(peek stageEDDegBound)
    );
    
    ------------------------------------------------------------------------
    -- (FUNCTION 7) Stage runner:
    --   - Run solve for the stage (writes input, runs Bertini, organizes outputs)
    --   - Apply "off" saturation (drop residuals) and "on" restrictions (enforce F)
    --   - In stage 2, optionally apply finer restrictions (if provided)
    ------------------------------------------------------------------------
    runComputationStage:=(stage,offPolyList,onPolyList)->(
        if stage==stageOne then runSolveInputFunction(stageOne,"input_first_solve")
        else runSolveInputFunction(stageTwo,"input_second_solve");
        --print("offPolyList",offPolyList);

        if stage===stageOne or NCO#"OutputType"=!="TestHomotopyConjectureGEDvUED" then runSaturateUnionFunction(offPolyList,stage);
        --print("WIN","SATURATE");
        --moveB'File(NCO#"Directory","member_points","filterFile",CopyB'File=>true);
        --print("onPolyList",onPolyList);

        runRestrictIntersectionFunction(onPolyList,stage);
        --print("WIN","RESTRICT");

        if stage===stageTwo and NCO#"FinerRestriction"=!={} then(
            --print("In stage 2, keep the critical points where each of these polynomials vanish ",NCO#"FinerRestriction");
            runRestrictIntersectionFunction(NCO#"FinerRestriction",stage)
        );
        print("We have completed stage ",stage)
    ); 

    -- Polynomials that should vanish (onPolyList) or should not vanish (offPolyList).
    -- offPolyList includes:
    --   - HX (hyperplane at infinity)
    --   - "L0" (first Lagrange multiplier) to remove trivial solutions at infinity
    --   - all generic hyperplanes used in rescaling
    offPolyList := {HX,"L0"}|((pairGeneralHyperplanes/(i->i#NameB'Section)));

    -- onPolyList are the homogenized F equations (in the Jacobian ring with HX)
    onPolyList := F/(i->homogenize(sub(i,jacRing),HX));

    -- Execute the requested stages, then return the count for the last executed stage.
    if isStageOne then runComputationStage(stageOne,offPolyList,onPolyList);
    if isStageTwo then runComputationStage(stageTwo,offPolyList,onPolyList);
    if isStageTwo then return stageEDDegBound#2 else if isStageOne then return stageEDDegBound#1
)

--##########################################################################--
-- Numerical ED Degree Methods
--##########################################################################--
numericalOptions = { TempDirectory => null }
numericWeightEDDegree = method(Options => numericalOptions)
numericWeightEDDegree(List, List, List, List) := o -> (F, G, data, weight) -> (
    dir := "";
    if o.TempDirectory === null then dir = temporaryFileName() else dir = o.TempDirectory;
    NCO := newNumericalComputationOptions(dir, F, G);
    NCO#"StartWeight" = weight;
    NCO#"TargetData" = NCO#"StartData" = data;
    homotopyEDDegree(NCO, "Weight", true, false)
)

numericUnitEDDegree = method(Options => numericalOptions)
numericUnitEDDegree(List, List) := o -> (F, G) -> numericWeightEDDegree(
    F, G,
    apply(#gens ring first F, i->randCC()),
    apply(#gens ring first F, i->1_(ring first F)),
    o
)

numericGenericEDDegree = method(Options => numericalOptions)
numericGenericEDDegree(List, List) := o -> (F, G) -> numericWeightEDDegree(
    F, G,
    apply(#gens ring first F, i->randCC()),
    apply(#gens ring first F, i->randCC()),
    o
)

-- aED of cardiod (4.6) or ellipse (4.5) from Draisma et al.
averageNumericEDDegree = method(Options => numericalOptions)
averageNumericEDDegree(List, List, List, ZZ) := o -> (F, G, L, n) -> (
    -- Initial run
    R := ring first F;
    NCO := newNumericalComputationOptions(F, G, L);
    homotopyEDDegree(NCO, "Data", true, false);

    -- homotopy to sampled data and average
    avgEDDeg := 0;
    for i to n do (
        NCO#"StartData" = apply(#gens R, i -> randomRR());
        homotopyEDDegree(NCO, "Data", false, true);

        -- count real critical points
        limitPoints := importMainDataFile(NCO#"Directory", NameMainDataFile => "stageTwo_main_data");
        realEDDeg := number(limitPoints, P -> (
            X := drop(drop(P#Coordinates, #G+1), -1);
            all(X, x -> imaginaryPart x < 1e-8)
        ));
        avgEDDeg = avgEDDeg + realEDDeg;
    );
    avgEDDeg/n
)
averageNumericEDDegree(List, List, ZZ) := o -> (F, G, n) -> averageNumericEDDegree(
    F, G, {},
    n,
    o
)

vanishTally = method() 
vanishTally(NumericalComputationOptions, Ideal, RR) := (NCO, Z, setTolerance) -> (
    limitPoints := importMainDataFile(NCO#"Directory", NameMainDataFile => "stageTwo_main_data");
    (F, G) := (NCO#"Model", NCO#"WitnessModel");    
    S := CC[gens ring first F];
    return tally apply(#limitPoints, s -> (	    
        p := limitPoints#s;
        X := drop(drop(p#Coordinates, #G+1), -1);
        xSub := apply(gens S, X, (i,j) -> i=>j);
        if setTolerance > norm sub(sub(gens Z,S), xSub) then {p#PathNumber}|p#PathsWithSameEndpoint else null 
    ))
)
vanishTally(NumericalComputationOptions, Ideal) := (NCO, Z) -> vanishTally(NCO, Z, 1e-8)
vanishTally(NumericalComputationOptions, List, RR) := (NCO, F, setTolerance) -> vanishTally(NCO, ideal F, setTolerance)
vanishTally(NumericalComputationOptions, List) := (NCO, F) -> vanishTally(NCO, ideal F)

end