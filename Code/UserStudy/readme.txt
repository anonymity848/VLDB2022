Readme (The User Study of Interactive Mining with Ordered and Unordered Attributes)
=========================
This package contains all source codes for a user study on
a. Algorithm DI 
	1. It only works for the special case of IOU
	2. The code is in folder DI
b. Algorithm EDI 
	1. It works for the general case of IOU. 
	2. The code is in folder EDI
c. Algorithm BS 
	1. It works for the general case of IOU.
	2. The code is in folder BS
d. Algorithm RH 
	1. It is an adapted existing algorithm.
	2. The code is in folder RH.
e. Algorithm ActiveRanking 
	1. It is an adapted existing algorithm.
	2. The code is in folder ActiveRanking.
f. Algorithm UH-Random 
	1. It is an adapted existing algorithm.
	2. The code is in folder UH.

Make sure there are folders called "input/", "output/" and "Result/" under the working directory.
They will be used for storing the input/output files and some results.

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
	Update "/usr/local/Cellar/glpk/5.0" to the path where you install the GLPK package
	
b. Execution
	./run

c. Input
	The used car dataset is shown in input/carReal.txt (whose format will be described in Appendix A).
	
d. Output
	1. The console output will be shown on the console (whose format will be described in Appendix B).
	2. The user log could be found in Result/XXX.txt (whose format will be described in Appendix C), 
	where XXX is your name that is typed in at the beginning of the user study. 

Example
=======
Try: ./run
Sample user log is provided: Result/UserA.txt



Appendix A. Format of Input File
------------------------------------
For input/carReal.txt, the dataset is described by 3 ordered attributes (price, year of production, used mileage) 
and 3 unordered attributes (length, width and number of seats).
The format of the first line is: n d_o d_u
	n - the dataset size, integer
	d_o - the number of ordered attributes in the dataset, integer
	d_u - the number of unordered attributes in the dataset, integer
The format of the following n lines is
	<price> <year of production> <used mileage> <length> <width> <number of seats>.
Each line corresponds to a car.
	

Appendix B. Format of Console Output
------------------------------------
There are two cases of console output. The first case of the output comes from the algorithms, where you need to 
answer each question during the interaction. The second case of the output comes from the algorithms, where you do 
not have to answer each question during the interaction.

a. The first case (each question needs to be answered)
	The recommending car system processes 5 algorithms. In each algorithm, the console output consists of several components: (1) interaction, (2) satisfaction, (3) boredness and (4) rank. 
	

	(1). Interaction
		In each algorithm, you will be presented consecutive questions. Each question consists of two cars and asks 
		you to choose the one you favor more. For example, you might see:
		************************************************************************************************************
		The 1th question
		Please choose the car you favor more:
		------------------------------------------------------------------------------------------------------------
                   Price(¥)         Year         Used KM      Length(mm)      Width(mm)          Seats
		------------------------------------------------------------------------------------------------------------
  		Option 1           335k           2013          57000           4540           1640              4
		------------------------------------------------------------------------------------------------------------
  		Option 2           698k           2018          25000           3495           1665              4
		------------------------------------------------------------------------------------------------------------
		Your choice: 
		************************************************************************************************************

		At the end of each algorithm, you will see the number of questions asked by this algorithm and the returned 
		car. For example, you might see:
		************************************************************************************************************
		Round Finish.        The number of questions asked: 11
		------------------------------------------------------------------------------------------------------------
    	    ID         Price(¥)         Year         Used KM      Length(mm)      Width(mm)          Seats
		------------------------------------------------------------------------------------------------------------
    	  1009           564k           2019          19199           3640           1690              5
		------------------------------------------------------------------------------------------------------------
		************************************************************************************************************

	(2). Satisfaction
		At the end of each algorithm, it will ask you to give a score to indicate how satisfied you are with the 
		returned car. For example, you might see:
		*****************************************************************************
		Please give a number from 1 to 10 (i.e., 1, 2, .., 10) to indicate how 
		satisfy you are when seeing the recommended cars. (Note: 10 denotes that 
		you are very satisfied with the recommended cars and 1 denotes that you are
		unsatisfied with the recommended cars):
		*****************************************************************************
	
	(3). Boredness
		At the end of each algorithm, it will ask you to give a score to indicate how bored you feel when you are 
		asked with XX questions in order to obtain your favorite car. For example, you might see:
		*****************************************************************************
		Please give a number from 1 to 10 (i.e., 1, 2, .., 10) to indicate how 
		bored you feel when you are asked with 41 questions in this round in 
		order to obtain your favorite car (Note: 10 denotes that you feel the 
		most bored and 1 denotes that you feel the least bored.): 
		*****************************************************************************

	(4). Rank
		After all the algorithms finished, it will ask you to compare the recommended cars returned by different 
		algorithms. You will be presented all the recommended cars and you need to give an order of the recommended 
		cars based on your preference. For example, you might see:
		************************************************************************************************************
		The recommended tuples: 
		------------------------------------------------------------------------------------------------------------
	                   	Price(¥)         Year         Used KM      Length(mm)      Width(mm)          Seats
		------------------------------------------------------------------------------------------------------------
		         1         407.9k           2020              0           4250           1700              5
		------------------------------------------------------------------------------------------------------------
		         2           420k           2019          19500           4337           1700              5
		------------------------------------------------------------------------------------------------------------
		Please give an order for the shown used car(s) (e.g., 1 2 3 4 5), where the first one 
		is the most preferred tuple and the last one is the least preferred tuple: 
		************************************************************************************************************


b. The second case (you do not need to answer each question)
	The recommending car system processes 3 algorithms. In each algorithm, the console output consists of several components: (1) interaction, (2) satisfaction, (3) boredness and (4) rank. 

	(1). Interaction
		In each algorithm, you will be presented consecutive questions. Each question displays two cars and asks 
		you to choose the one you favor more. If you cannot decide which one is better, you can skip this question.
		For example, you might see:
		************************************************************************************************************
		The 1th question
		Please choose the car you favor more:
		------------------------------------------------------------------------------------------------------------
    	               Price(¥)         Year         Used KM      Length(mm)      Width(mm)          Seats
		------------------------------------------------------------------------------------------------------------
  		Option 1           564k           2019          19199           3640           1690              5
		------------------------------------------------------------------------------------------------------------
  		Option 2           495k           2018          24000           3740           1713              5
		------------------------------------------------------------------------------------------------------------
  		Option 3           Hard to decide which one is more preferred.
		------------------------------------------------------------------------------------------------------------
		Your choice:
		************************************************************************************************************

		At the end of each algorithm, you will see the number of questions asked by this algorithm and the returned cars. Note that since you may not answer all the questions, some algorithms may recommend more than one cars. For example, you might see:
		************************************************************************************************************
		Round Finish.        The number of questions asked: 11
		------------------------------------------------------------------------------------------------------------
      		 ID         Price(¥)         Year         Used KM      Length(mm)      Width(mm)          Seats
		------------------------------------------------------------------------------------------------------------
      		1134           510k           2013         100000           4213           2133              5
		------------------------------------------------------------------------------------------------------------
		************************************************************************************************************

	(2). Satisfaction
		At the end of each algorithm, it will ask you to give a score to indicate how satisfied you are with the returned cars. For example, you might see:
		*****************************************************************************
		Please give a number from 1 to 10 (i.e., 1, 2, .., 10) to indicate how 
		satisfy you are when seeing the recommended cars. (Note: 10 denotes that 
		you are very satisfied with the recommended cars and 1 denotes that you are
		unsatisfied with the recommended cars):
		*****************************************************************************
	
	(3). Boredness
		At the end of each algorithm, it will ask you to give a score to indicate how bored you feel when you are asked with XX questions in order to obtain your favorite car. For example, you might see:
		*****************************************************************************
		Please give a number from 1 to 10 (i.e., 1, 2, .., 10) to indicate how 
		bored you feel when you are asked with 41 questions in this round in 
		order to obtain your favorite car (Note: 10 denotes that you feel the 
		most bored and 1 denotes that you feel the least bored.): 
		*****************************************************************************

	(4). Rank
		After all the algorithms finished, it will ask you to compare the recommended cars returned by different 
		algorithms. Each algorithm will recommend a set of cars. You need to give an order of these sets based on 
		your preference. For example, you might see:
		************************************************************************************************************
		The recommended tuple sets: 
		------------------------------------------------------------------------------------------------------------
     		 Set 1         Price(¥)         Year         Used KM      Length(mm)      Width(mm)          Seats
		------------------------------------------------------------------------------------------------------------
		         1          1749k           2017          27971           4895           1928              9
		------------------------------------------------------------------------------------------------------------


		------------------------------------------------------------------------------------------------------------
		     Set 2         Price(¥)         Year         Used KM      Length(mm)      Width(mm)          Seats
		------------------------------------------------------------------------------------------------------------
		         1           564k           2019          19199           3640           1690              5
		------------------------------------------------------------------------------------------------------------


		------------------------------------------------------------------------------------------------------------
		     Set 3         Price(¥)         Year         Used KM      Length(mm)      Width(mm)          Seats
		------------------------------------------------------------------------------------------------------------
		         1           510k           2013         100000           4213           2133              5
		------------------------------------------------------------------------------------------------------------
		Please give an order for the shown sets (e.g., 1 2 3), where the first one 
		is the most preferred set and the last one is the least preferred set: 1
		************************************************************************************************************



Appendix C. Format of User Log
------------------------------------
	It contains two parts: (a) The first part contains the results of all the algorithms, where the questions asked 
	by each algorithm needs to be answered. (b) The second part contains the results of all the algorithms, where the questions can be unanswered. A sample of the user log is shown in Result/UserA.txt.


a. The first part contains (1) the result of algorithms and (2) the result of rank
	
	(1). The result of algorithms
		Lines 1-2 show the algorithm name and the number of questions asked.
		Lines 3-5 show the recommended car.
		Lines 6-8 show the evaluation result (the degree of satisfaction and the degree of boredness).
		For example, you might see:
		************************************************************************************************************
		Algorithm: ActiveRanking   Question: 41
		------------------------------------------------------------------------------------------------------------
		        ID         Price(¥)         Year         Used KM      Length(mm)      Width(mm)          Seats
		------------------------------------------------------------------------------------------------------------
		      1497         407.9k           2020              0           4250           1700              5
		------------------------------------------------------------------------------------------------------------
		Satisfaction: 3
		Boredness: 10
		************************************************************************************************************
	(2). The result of rank
		Line 1 is the title "order".
		The rest lines show the order of the recommended cars (which are returned by different algorithms) based on 
		the user preference, where the first one is the most preferred car and the last one is the least preferred 
		car. For example, you might see:
		************************************************************************************************************
		order: 
		------------------------------------------------------------------------------------------------------------
		         1         407.9k           2020              0           4250           1700              5
		------------------------------------------------------------------------------------------------------------
		         2           420k           2019          19500           4337           1700              5
		------------------------------------------------------------------------------------------------------------
		************************************************************************************************************


b. The second part contains (1) the result of algorithms and (2) the result of rank
	
	(1). The result of algorithms
		Lines 1-2 show the algorithm name and the number of questions asked.
		Lines 3-5 show the recommended car.
		Lines 6-8 show the evaluation result (the degree of satisfaction and the degree of boredness).
		For example, you might see:
		************************************************************************************************************
		Algorithm: BS   Question: 11
		------------------------------------------------------------------------------------------------------------
		        ID         Price(¥)         Year         Used KM      Length(mm)      Width(mm)          Seats
		------------------------------------------------------------------------------------------------------------
		      1497         407.9k           2020              0           4250           1700              5
		------------------------------------------------------------------------------------------------------------
		Satisfaction: 3
		Boredness: 2
		************************************************************************************************************
	
	(2). The result of rank
		Line 1 is the title "order".
		Line 2 shows the order of the sets given by the user. Each set contains the cars recommended by an 
		algorithm. The first one is the most preferred set and the last one is the least preferred set. For 
		example, you might see:
		************************************************************************************************************
		order: 1  2  3  
		************************************************************************************************************



