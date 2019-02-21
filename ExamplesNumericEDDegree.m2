--restart

--Open M2 in the directory containing NumericEDDegree.m2 or add the directory containing NumericEDDegree.m2 to the path. 
--path=prepend("/Users/jo/Dropbox/Euclidean distance degree/Projective varieties/ComputingEDDegree",path)
load"NumericEDDegree.m2"

--INPUTS: 

--storeBM2Files="/Users/jo/Desktop" 
----The input files are stored in the directory storeBM2Files which as a default is set to a temporary directory. Spaces are not allowed in the string. 
----If Bertini is installed to run with M2, then the runBertini commands and readFile will work. 
----The number of nonsingular-solutions equals the ED degree and found in bertini_session.log

--R 
----The coordinate ring of a codimension c variety X. 

--cod 
----The codimension of X 

--F 
----A list of polynomials that generate the ideal of X (assume X is generically reduced)


--Example 4.1.
--storeBM2Files="/Users/jo/Desktop/BertiniOutputFiles/Ex4.1"
R=QQ[x0,x1,x2];
F={x0^2*x2 - x1^2*(x1 + x2)};
cod=1--The codimension of the variety. 
numericEDDegree(storeBM2Files,cod,F)--ED degree is 7
--runBertini(storeBM2Files)
--readFile(storeBM2Files)
--runBertiniEDDegree(storeBM2Files)

--Example 4.2.
--storeBM2Files="/Users/jo/Desktop/BertiniOutputFiles/Ex4.2"
R=QQ[I][x0,x1,x2]; --I plays the role of the imaginary unit in Bertini. 
F={x0^2*x1 -(x1 - I*x2)^2*x2};
cod=1--The codimension of the variety. 
numericEDDegree(storeBM2Files,cod,F)--ED degree is 7
--runBertini(storeBM2Files)
--readFile(storeBM2Files)
--runBertiniEDDegree(storeBM2Files)

--Example 4.3.
--storeBM2Files="/Users/jo/Desktop/BertiniOutputFiles/Ex4.3"
R=QQ[I][x0,x1,x2]; --I plays the role of the imaginary unit in Bertini. 
F={x0^3 - (I*x0^2 + x1^2)*x2};
cod=1--The codimension of the variety. 
numericEDDegree(storeBM2Files,cod,F)--ED degree is 6
--runBertini(storeBM2Files)
--readFile(storeBM2Files)
--runBertiniEDDegree(storeBM2Files)

--Example 4.4.
--storeBM2Files="/Users/jo/Desktop/BertiniOutputFiles/Ex4.4"
R=QQ[x0,x1,x2,x3]; --I plays the role of the imaginary unit in Bertini. 
F={x0^2*x1-x2*x3^2};
cod=1--The codimension of the variety. 
numericEDDegree(storeBM2Files,cod,F)--ED degree is 10
--runBertini(storeBM2Files)
--readFile(storeBM2Files)
--runBertiniEDDegree(storeBM2Files)

--PROJECTIVE ED DEGREES.
----We can also use our new intrinsic formulation for projective ED Degrees.

--storeBM2Files="/Users/jo/Desktop/BertiniOutputFiles/ProjEx4.1"
R=QQ[x0,x1,x2,x3];
F={x0^2*x2 - x1^2*(x1 + x2)};
numericProjectiveEDDegree(storeBM2Files,F)--ED degree is 7
runBertiniProjectiveEDDegree(storeBM2Files)


--Projective Example 4.2.
--storeBM2Files="/Users/jo/Desktop/BertiniOutputFiles/ProjEx4.2"
R=QQ[I][x0,x1,x2]; --I plays the role of the imaginary unit in Bertini. 
F={x0^2*x1 -(x1 - I*x2)^2*x2};
numericProjectiveEDDegree(storeBM2Files,F)--ED degree is 7
runBertiniProjectiveEDDegree(storeBM2Files)

--Projective Example 4.3.
--storeBM2Files="/Users/jo/Desktop/BertiniOutputFiles/ProjEx4.3"
R=QQ[I][x0,x1,x2]; --I plays the role of the imaginary unit in Bertini. 
F={x0^3 - (I*x0^2 + x1^2)*x2};
numericProjectiveEDDegree(storeBM2Files,F)--ED degree is 6
runBertiniProjectiveEDDegree(storeBM2Files)

--runBertini(storeBM2Files)
--readFile(storeBM2Files)

--Projective Example 4.4.
--storeBM2Files="/Users/jo/Desktop/BertiniOutputFiles/ProjEx4.4"
R=QQ[x0,x1,x2,x3]; 
F={x0^2*x1-x2*x3^2};
numericProjectiveEDDegree(storeBM2Files,F)--ED degree is 10
runBertiniProjectiveEDDegree(storeBM2Files)

