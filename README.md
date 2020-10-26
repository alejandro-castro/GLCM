<h1>Gray-Level Co-Occurence Matrix (GLCM)</h1>
This is an efficient implementation of Gray-Level Co-Occurence Matrix (GLCM) in C/C++ (OpenCV).  <br>

This software is a total reworked version of the work found here https://github.com/luizbuschetto/GLCM

This work is an hybrid combination of using only static memory allocation and the use of linked lists[4] to avoid processing points with zero probabilities.

It calculates 20 features¹ are:

Autocorrelation: [2]                      				(autoc)<br>
Contrast: [1,2]                    			            (contr)<br>
Correlation: [1,2]                        				(corr)<br>
Cluster Prominence: [2]                   			    (cprom)<br>
Cluster Shade: [2]                       	 			(cshad)<br>
Dissimilarity: [2]                        				(dissi)<br>
Energy: MATLAB / [1,2]                    				(energ)<br>
Entropy: [2]                              				(entro)<br>
Maximum probability: [2]                  			    (maxpr)<br>
Sum average [1]                           				(savgh)<br>
Sum variance [1]                         			 	(svarh)<br>
Sum entropy [1]                          				(senth)<br>
Difference variance [1]                   				(dvarh)<br>
Difference entropy [1]                    				(denth)<br>
Information measure of correlation1 [1]   	            (inf1h)<br>
Informaiton measure of correlation2 [1]   	            (inf2h)<br>
Inverse difference normalized (INN) [3]   		        (indnc)<br>
Inverse difference moment normalized [3]  	            (idmnc)<br>
Variance of X                                           (varX)<br>
Variance of Y                                           (varY)<br>

The last 2 features are reworked versions of the sum of squares proposed in [1].

¹ The number inside the brackets shows where the formula can be found (in the papers below). 

<h3>References:</h3>
[1] R. M. Haralick, K. Shanmugam, and I. Dinstein, Textural Features of Image Classification, IEEE Transactions on Systems, Man and Cybernetics, vol. SMC-3, no. 6, Nov. 1973;<br>
[2] L. Soh and C. Tsatsoulis, Texture Analysis of SAR Sea Ice Imagery Using Gray Level Co-Occurrence Matrices, IEEE Transactions on Geoscience and Remote Sensing, vol. 37, no. 2, March 1999;<br>
[3] D A. Clausi, An analysis of co-occurrence texture statistics as a function of grey level quantization, Can. J. Remote Sensing, vol. 28, no.1, pp. 45-62, 2002;<br>
[4] D. A. Clausi and M. E. Jernigan, "A fast method to determine co-occurrence texture features," in IEEE Transactions on Geoscience and Remote Sensing, vol. 36, no. 1, pp. 298-300, Jan. 1998;


<h2>How to use:</h2>
<i> ./glcm image_path number_of_levels </i>

Example: <i> ./glcm BadQuality1.png 8 </i>

The algorithm will move through the entire image extracting a window of 10 x 10 and it will calculate the features from each window with numLevels = 8 and with a vertical distance of 1 pixel.
The second parameter (number of levels) is optional. If you don't want to specify this number, the image will be processed with the default value: 256 levels.

<h2>Installation:</h2>
This software needs a minimum cmake of 3.17 and a working version of C++ 14. It has been tested using opencv 4.5.0.
It also needs MATLAB to run the test.m to check the Haralick features, but if it is not necessary, the MATLAB part can be erased in the CMakeLists.txt file.

To correctly use this software the following lines should be added to the .bash_profile or .bash_rc in the HOME directory. Take in account that the paths that are shown below should be changed with the correspondent ones in your system.
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:/Users/alejandro/Downloads/TesisMaestria/Software/SoftwareTerceros/opencv/install"

The following two lines are only need to save the texture images in mat file format to compare them with the Test.m script.
export PATH="/Applications/MATLAB_R2019b.app/bin:$PATH"
export DYLD_LIBRARY_PATH="/Applications/MATLAB_R2019b.app/bin/maci64:/Applications/MATLAB_R2019b.app/sys/os/maci64:$DYLD_LIBRARY_PATH"

NOTE: You should also change the absolute paths in the CMakeLists.txt file to the correct ones in your system.

Finally, if you has choose to not compare the Haralick features to the ones provided by MATLAB you should comment the line 3 y 55 of the file main.cpp.

Once you has finished the installation, you can process the cmake file using the command cmake . when you are in the folder where this project is. You don't need to do this more than one time.
Later you can compile the source files with the command make, again this works only if you are in the folder where this project is. With this configuration and compilation now you can use the software to process your images.