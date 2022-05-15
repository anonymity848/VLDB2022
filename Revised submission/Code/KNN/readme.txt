Readme (Interactive Mining with Ordered and Unordered Attributes)
=========================
This package contains all source codes for 
a. Algorithm DI 
	1. It only works for the special case of IOU
	2. The code is in folder DI
b. Algorithm EDI 
	1. It works for the general case of IOU. 
	2. The code is in folder EDI
c. Algorithm BS 
	1. It works for the general case of IOU.
	2. The code is in folder BS
d. Algorithm Baseline
	1. It is an adaptation of BS which utilizes the Vonoroi Diagram. 
	2. The code is in folder Baseline. 
e. Algorithm RH 
	1. It is an adapted existing algorithm.
	2. The code is in folder RH.
f. Algorithm ActiveRanking 
	1. It is an adapted existing algorithm.
	2. The code is in folder ActiveRanking.
g. Algorithm UtilityApprox 
	1. It is an adapted existing algorithm.
	2. The code is in folder UtilityApprox.
h. Algorithm UH-Random 
	1. It is an adapted existing algorithm.
	2. The code is in folder UH.
i. Algorithm EDITopk
	1. It is an extension of algorithm EDI. It return k nearest points. 
	2. The code is in folder EDI.
j. Algorithm BSTopk
	1. It is an extension of algorithm BS. It return k nearest points. 
	2. The code is in folder BS.
k. Algorithm ActiveRankingTopk
	1. It is an extension of algorithm ActiveRanking. It return k nearest points. 
	2. The code is in folder ActiveRanking.

Make sure there is a folder called "input/", a folder called "output/" and a file called 
"config.txt" under the working directory.
They will be used for storing the input/output files, some intermediate results and the 
input parameters.

Usage Step
==========
a. Compilation
	mkdir build
	cd build
	cmake ..
	make

	You will need to install the GLPK package (for solving LPs) at first.
	See GLPK webpage <http://www.gnu.org/software/glpk/glpk.html>.
	Then update the path in CMakeLists.txt
		set(INC_DIR /usr/local/Cellar/glpk/5.0/include)
		set(LINK_DIR /usr/local/Cellar/glpk/5.0/lib)
	Update path "/usr/local/Cellar/glpk/5.0" to the path you install the GLPK package
	
b. Execution
	./run

c. Config
	The config file contains the input parameters (whose format will be described in Appendix A).

c. Input
	The input file contains the dataset (whose format will be described in Appendix B).
	
d. Output
	The output will be shown on the console (whose format will be described in Appendix C).

Example
=======
Sample input (input/Anti11d.txt) are provided. The dataset is described by one ordered 
attribute and one unordered attribute
Try: ./run



Appendix A. Format of Config File
------------------------------------
The format is: DatasetName type Beta gamma k e[1] e[2] ... e[d]
DatasetName - the name of the dataset
type - the output type (2 means it only shows the final result. 1 means it shows the middle result and the final result)
Beta - a parameter used in our algorithms
gamma - a parameter used in our algorithms
k - the number of returned points
e[1] - the first attribute value of the user's expected point 
e[2] - the second attribute value of the user's expected point
...
e[d] - the d-th attribute value of the user's expected point
For example, you might see
------------------------------------------
Anti11d.txt	 2  4  4  1  1000  120.018
------------------------------------------


Appendix B. Format of Input File
------------------------------------
The format of the first line is: n d_o d_u
n - the dataset size, integer
d_o - the number of ordered attributes in the dataset, integer
d_u - the number of unordered attributes in the dataset, integer
The format of the following n lines is
-----------------------------------------------------------------
<ordered attribute 1> <ordered attribute 2> ... <<ordered attribute d_o> 
<unordered attribute 1> <unordered attribute 2> ... <unordered attribute d_u> 
-----------------------------------------------------------------
Each line corresponds to a point.


Appendix C. Format of Console Output
-----------------------------------------------------------------
The format of the output is
--------------------------------------------------------------------------------------
|           Algorithm | # of Questions |  Preprocessing |    Interaction | Point #ID |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|        Ground Truth |              - |              - |              - |    PtID-0 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|            Baseline |      Q-count-1 |      PreTime-1 |    InterTime-1 |    PtID-1 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|                  DI |      Q-count-2 |      PreTime-2 |    InterTime-2 |    PtID-2 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|                 EDI |      Q-count-3 |      PreTime-3 |    InterTime-3 |    PtID-3 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|                  BS |      Q-count-4 |      PreTime-4 |    InterTime-4 |    PtID-4 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|                  RH |      Q-count-5 |      PreTime-5 |    InterTime-5 |    PtID-5 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|       ActiveRanking |      Q-count-6 |      PreTime-6 |    InterTime-6 |    PtID-6 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|       UtilityApprox |      Q-count-7 |      PreTime-7 |    InterTime-7 |    PtID-7 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|           UH-Random |      Q-count-8 |      PreTime-8 |    InterTime-8 |    PtID-8 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|              BSTopk |      Q-count-9 |      PreTime-9 |    InterTime-9 |    PtID-9 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|             EDITopk |      Q-count-10|      PreTime-10|    InterTime-10|    PtID-10|
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|   ActiveRankingTopk |      Q-count-11|      PreTime-11|    InterTime-11|    PtID-11|
--------------------------------------------------------------------------------------
where PtID-0 is the point ID of the user's favorite point (ground truth),
PtID-1 is the point ID of the user's favorite point returned by algorithm Baseline,
PtID-2 is the point ID of the user's favorite point returned by algorithm DI,
PtID-3 is the point ID of the user's favorite point returned by algorithm EDI,
PtID-4 is the point ID of the user's favorite point returned by algorithm BS,
PtID-5 is the point ID of the user's favorite point returned by algorithm RH,
PtID-6 is the point ID of the user's favorite point returned by algorithm ActiveRanking,
PtID-7 is the point ID of the user's favorite point returned by algorithm UtilityApprox,
PtID-8 is the point ID of the user's favorite point returned by algorithm UH-Random,
PtID-9 is the point ID of the user's favorite point returned by algorithm BSTopk,
PtID-10 is the point ID of the user's favorite point returned by algorithm EDITopk,
PtID-11 is the point ID of the user's favorite point returned by algorithm ActiveRankingTopk,
Q-count-1 is the number of questions asked by algorithm Baseline,
Q-count-2 is the number of questions asked by algorithm DI,
Q-count-3 is the number of questions asked by algorithm EDI,
Q-count-4 is the number of questions asked by algorithm BS,
Q-count-5 is the number of questions asked by algorithm RH,
Q-count-6 is the number of questions asked by algorithm ActiveRanking,
Q-count-7 is the number of questions asked by algorithm UtilityApprox,
Q-count-7 is the number of questions asked by algorithm UH-Random,
Q-count-9 is the number of questions asked by algorithm BSTopk,
Q-count-10 is the number of questions asked by algorithm EDITopk,
Q-count-11 is the number of questions asked by algorithm ActiveRankingTopk,
PreTime-1 is the preprocessing time of algorithm Baseline,
PreTime-2 is the preprocessing time of algorithm DI,
PreTime-3 is the preprocessing time of algorithm EDI,
PreTime-4 is the preprocessing time of algorithm BS,
PreTime-5 is the preprocessing time of algorithm RH,
PreTime-6 is the preprocessing time of algorithm ActiveRanking,
PreTime-7 is the preprocessing time of algorithm UtilityApprox,
PreTime-8 is the preprocessing time of algorithm UH-Random,
PreTime-9 is the preprocessing time of algorithm BSTopk,
PreTime-10 is the preprocessing time of algorithm EDITopk,
PreTime-11 is the preprocessing time of algorithm ActiveRankingTopk,
InterTime-1 is the interaction time of algorithm Baseline,
InterTime-2 is the interaction time of algorithm DI,
InterTime-3 is the interaction time of algorithm EDI,
InterTime-4 is the interaction time of algorithm BS,
InterTime-5 is the interaction time of algorithm RH,
InterTime-6 is the interaction time of algorithm ActiveRanking,
InterTime-7 is the interaction time of algorithm UtilityApprox,
InterTime-8 is the interaction time of algorithm UH-Random.
InterTime-9 is the interaction time of algorithm BSTopk,
InterTime-10 is the interaction time of algorithm EDITopk,
InterTime-11 is the interaction time of algorithm ActiveRankingTopk.

For example, you might see:
--------------------------------------------------------------------------------------
|           Algorithm | # of Questions |  Preprocessing |    Interaction | Point #ID |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|        Ground Truth |              - |              - |              - |     60625 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|            Baseline |              8 |      12.606395 |       0.000072 |     60625 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|                  DI |              6 |       7.764545 |       0.000022 |     60625 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|                 EDI |              6 |       6.953771 |       0.000051 |     60625 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|                  BS |              6 |       6.997851 |       0.000094 |     60625 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|                  RH |              9 |       8.422067 |       0.000735 |     60625 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|       ActiveRanking |             20 |       6.721187 |       0.012358 |     60625 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|       UtilityApprox |             10 |       0.010086 |       0.097760 |     60625 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|           UH-Random |              8 |       6.688584 |       0.028314 |     60625 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|              BSTopk |              6 |       6.785224 |       0.000232 |     60625 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|             EDITopk |              6 |       7.047608 |       0.000242 |     60625 |
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
|   ActiveRankingTopk |             20 |       6.868985 |       0.013542 |     60625 |
--------------------------------------------------------------------------------------