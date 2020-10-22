
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <map>
#include <string>
#include <iostream>
#define numberOfFeaturesToBeUsed 22

using namespace std;//I don't think it's a good practice, but it reduces lines of coding
using namespace cv;

const std::string haralickFeatureNames[numberOfFeaturesToBeUsed] = {"autoc", "contr", "corrm", "corrp" , "cprom", "cshad", "dissi", "energ", "entro", 
"homom", "homop", "maxpr", "sosvh", "savgh", "svarh", "senth", "dvarh", "denth", "inf1h", "inf2h", "indnc", "idmnc"};


std::map<std::string, float>  glcm(cv::Mat &img, int numLevels);
void FillPixelWithHaralickFeatures(Mat img[], map <string, float> haralickFeatures, int y, int x);
void initializeArrayOfImages(Mat img[], int arraySize, int dimY, int dimX);