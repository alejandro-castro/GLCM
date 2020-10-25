#include "mat.h"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <cstring>
using namespace cv;

void saveFloatImageAsMatFile(double *image, int arraySize, int dimY, int dimX);