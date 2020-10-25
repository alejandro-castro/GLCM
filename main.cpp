#include <cmath>
#include "glcm.h"
#include "CVMatlabInterface.h"
#include <chrono> 

#define WINDOW_SIZE 10
#define defaultNumLevels 256

int main(int argc, char* argv[]){
  std::chrono::time_point<std::chrono::high_resolution_clock>  start = std::chrono::high_resolution_clock::now(); 
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
      float new_level = slope * img.at<uint8_t>(i,j);
      img.at<uint8_t>(i,j) = floor(new_level);
    }
  }

  //Getting dimensions of texture images
  int numberOfWindowsX = img.cols - WINDOW_SIZE + 1;//Without padding and with stride 1
  int numberOfWindowsY = img.rows - WINDOW_SIZE + 1;
  
  //Initialization of the texture images
  double *imgResult;  int sizeOfImgResult = numberOfFeaturesToBeUsed * numberOfWindowsY * numberOfWindowsX;
  int dimImgResult[3] = {numberOfFeaturesToBeUsed, numberOfWindowsY, numberOfWindowsX};
  imgResult = (double*)malloc(sizeof(double) * sizeOfImgResult);

  double *temp = imgResult;
  for (int i=0; i < sizeOfImgResult;i ++){
    *temp = 0.0; 
    temp++;
  }


  for (int j=0; j<numberOfWindowsY; j++){
    for (int i=0; i<numberOfWindowsX; i++){
      Mat ROI = img(Range(j, j+WINDOW_SIZE), Range(i, i+WINDOW_SIZE)); 
      glcm(ROI, numLevels, j, i, imgResult, dimImgResult);
      }
  }

  std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now(); 
  std::chrono::duration<double, std::milli>  diffTime = end - start; 
  std::cout << "The running time of this script is: "<<diffTime.count() << " milliseconds."<<endl;

  saveFloatImageAsMatFile(imgResult, numberOfFeaturesToBeUsed, numberOfWindowsY, numberOfWindowsX);

  return 0;
}

