#include <cmath>
#include "glcm.h"
#include "CVMatlabInterface.h"

#define WINDOW_SIZE 10
#define defaultNumLevels 256

int main(int argc, char* argv[]){
  //Reading image
  Mat img = imread(argv[1]);
  if(img.empty()){
    cout << "No image!";
    return 0;
  }
  //Change to gray scale image
  cvtColor(img, img, COLOR_BGR2GRAY);

  //Getting number of levels
  int numLevels = defaultNumLevels;
  if (argc > 2) numLevels = atoi(argv[2]);

  // Normalizing
  float slope = float(numLevels) / 256; // we calculate the inverse of the number of points that should correspond to each new quantization level
  for(int i = 0; i < img.rows; i++){
    for(int j = 0; j < img.cols; j++){
      img.at<uchar>(i,j) = floor(slope * img.at<uchar>(i,j));
    }
  }

  //Getting dimensions of texture images
  int numberOfWindowsX = img.cols - WINDOW_SIZE + 1;//Without padding and with stride 1
  int numberOfWindowsY = img.rows - WINDOW_SIZE + 1;
  
  //Initialization of the texture images
  Mat imgResult[numberOfFeaturesToBeUsed];
  initializeArrayOfImages(imgResult, numberOfFeaturesToBeUsed, numberOfWindowsY, numberOfWindowsX);

  map <string, float>  ROIresult;
  
  for (int j=0; j<numberOfWindowsY; j++){
    for (int i=0; i<numberOfWindowsX; i++){
      Mat ROI = img(Range(j, j+WINDOW_SIZE), Range(i, i+WINDOW_SIZE)); 
      ROIresult = glcm(ROI, numLevels);
      FillPixelWithHaralickFeatures(imgResult, ROIresult, j, i);
    }
  }

  saveFloatImageAsMatFile(imgResult, numberOfFeaturesToBeUsed, numberOfWindowsY, numberOfWindowsX);
  return 0;
}

