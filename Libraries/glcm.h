
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <map>
#include <string>
#include <iostream>
#define numberOfFeaturesToBeUsed 20

using namespace std;//I don't think it's a good practice, but it reduces lines of coding
using namespace cv;

const std::string haralickFeatureNames[numberOfFeaturesToBeUsed] = {"autoc", "contr", "corr", "cprom", "cshad", "dissi", "energ", "entro", "maxpr", "savgh", 
"svarh", "senth", "dvarh", "denth", "inf1h", "inf2h", "indnc", "idmnc", "varX", "varY"};


void glcm(cv::Mat img, int numLevels, int pos_y, int pos_x, double *imgResult, int dims[], bool print);

