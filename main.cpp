#include <iostream>
#include <math.h>
#include <string>
#include <map>
#include <vector>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "mat.h"

#define WINDOW_SIZE 10
#define defaultNumLevels 256
#define numberOfFeaturesToBeUsed 22
using namespace std;
using namespace cv;

const string haralickFeatureNames[22] = {"autoc", "contr", "corrm", "corrp" , "cprom", "cshad", "dissi", "energ", "entro", 
"homom", "homop", "maxpr", "sosvh", "savgh", "svarh", "senth", "dvarh", "denth", "inf1h", "inf2h", "indnc", "idmnc"};

map <String, float>  glcm(Mat &img, int numLevels); //map <String, float> 
void FillPixelWithHaralickFeatures(Mat img, map <String, float> haralickFeatures, int y, int x);
int address(int i, int j, int k, MatStep step){ return i*step[0]+j*step[1]+k*step[2];}
void saveFloatImageAsMatFile(double* image, int size[]);


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
   
  //Getting dimensions of texture images
  int numberOfWindowsX = img.cols - WINDOW_SIZE + 1;//Without padding and with stride 1
  int numberOfWindowsY = img.rows - WINDOW_SIZE + 1;
  
  int size[] = {22, numberOfWindowsY, numberOfWindowsX};
  Mat imgResult(3, size, CV_64F, Scalar::all(0));

  map <String, float>  ROIresult;
  

  for (int j=0; j<numberOfWindowsY; j++){
    for (int i=0; i<numberOfWindowsX; i++){
      Mat ROI = img(Range(j, j+WINDOW_SIZE), Range(i, i+WINDOW_SIZE));
      ROIresult = glcm(ROI, numLevels);
      cout << ROIresult["corrm"]<<endl;
      FillPixelWithHaralickFeatures(imgResult, ROIresult, j, i);
    }
  }
//Range ranges[3] = {Range::all(), Range::all(), Range(0,1)};
//Mat autocimg = imgResult(ranges);

  double* dataPointer = (double*)imgResult.data;
  saveFloatImageAsMatFile(dataPointer, size);
return 0;
}

void saveFloatImageAsMatFile(double* image, int size[]){
  const char *file = "Image1.mat";
  
  //Creacion de archivo mat
  MATFile *pmat = matOpen(file, "w");
  mxArray *pa1;

  for(int i=0;i < size[0];i++){
    //Getting the Haralick feature one by one
    //Range ranges[3] = {Range::all(), Range::all(), Range(i,i+1)};
    //Mat img = image(ranges);
    
    //Creating the matrix to be save in mat file and copying the data from image
    pa1 = mxCreateDoubleMatrix(size[1],size[2],mxREAL);
    memcpy((void *)(mxGetPr(pa1)), (void *)(image+i*size[1]*size[2]), sizeof(double)*size[1]*size[2]); 
    const char *name = haralickFeatureNames[i].c_str();
    matPutVariable(pmat, name, pa1);
  }
  mxDestroyArray(pa1);
  matClose(pmat);
}
void FillPixelWithHaralickFeatures(Mat img, map <String, float> haralickFeatures, int y, int x){
  for (int i=0;i<22;i++){
    string name = haralickFeatureNames[i];
    img.at<double>(i, y, x)  = haralickFeatures[name];
  }
}

map <String, float>  glcm(Mat &img, int numLevels){
vector <String> haralickFeatureNames = {"autoc", "contr", "corrm", "corrp" , "cprom", "cshad", "dissi", "energ", "entro", 
"homom", "homop", "maxpr", "sosvh", "savgh", "svarh", "senth", "dvarh", "denth", "inf1h", "inf2h", "indnc", "idmnc"}; //Va a ser usado para decidir si se calcula ese feature o no

map <String, float> haralickFeaturesValues;


  Mat glcm = Mat::zeros(numLevels, numLevels, CV_32F);


  // Normalizing
  float slope = float(numLevels) / 255; // "255" is the maximum pixel value for the image type -> 0 ~ 255
  float intercept = 1 - (slope * 0);    // "0"   is the minimum pixel value for the image type -> 0 ~ 255

  for(int i = 0; i < img.rows; i++){
    for(int j = 0; j < img.cols; j++){
      img.at<uchar>(i,j) = floor((slope * img.at<uchar>(i,j)) + 1);
    }
}

  // Creating GLCM matrix with "numLevels", radius = 1 and in the horizontal direction
  for(int i = 0; i < img.rows-1; i++){
    for(int j = 0; j < img.cols; j++)
    {
      glcm.at<float>(img.at<uchar>(i,j) - 1, img.at<uchar>(i+1,j) - 1) += 1; //glcm.at<float>(img.at<uchar>(i,j) - 1, img.at<uchar>(i,j+1) - 1) += 1;
    }      
  }

  // Normalizing GLCM matrix for parameter determination
  glcm = glcm / sum(glcm)[0];

  // Means
  float mu_i = 0, mu_j = 0;
  float mu = 0;
  for(int i = 0; i < numLevels; i++){
    for(int j = 0; j < numLevels; j++)
    {
      mu = mu + glcm.at<float>(i,j);
      mu_i = mu_i + ((i+1) * glcm.at<float>(i,j));
      mu_j = mu_j + ((j+1) * glcm.at<float>(i,j));
    }
  }

  mu = mu / pow(numLevels, 2);

  // Standard Deviation
  float sigma_i = 0, sigma_j = 0;
  vector<float> sum_pixels((2 * numLevels)-1); // sum_pixels[i + j] = p_{x + y}(k), where k = i + j
  vector<float> dif_pixels(2 * numLevels);     // dif_pixels[i + j] = p_{x - y}(k), where k = |i + j|
  vector<float> px(numLevels), py(numLevels);  // Marginal-probability matrix obtained by summing the rows of p(i, j) (Matrix gl) [2]

  for(int i = 0; i < numLevels; i++)
  {
    for(int j = 0; j < numLevels; j++)
    {
      sigma_i = sigma_i + (((i+1) - mu_i) * ((i+1) - mu_i) * glcm.at<float>(i,j));
      sigma_j = sigma_j + (((i+1) - mu_j) * ((i+1) - mu_j) * glcm.at<float>(i,j));

      // Implementing savgh - Sum Average [1]
      if (i + j >= 2)
        sum_pixels[i + j] = sum_pixels[i + j] + glcm.at<float>(i,j);

      dif_pixels[abs(i - j)] = dif_pixels[abs(i - j)] + glcm.at<float>(i,j);
      px[j] = px[j] + glcm.at<float>(i,j);
      py[i] = py[i] + glcm.at<float>(i,j);
    }
  }

  sigma_i = sqrt(sigma_i);
  sigma_j = sqrt(sigma_j);

  float savgh = 0, senth = 0;
  for(int i = 0; i < (2 * numLevels) - 1; i++)
  {
    savgh = savgh + ((i+2) * sum_pixels[i]); 
    senth = senth - (sum_pixels[i] * log(sum_pixels[i] + 0.0000000000001)); //¹
  }

  float svarh = 0;
  for(int i = 0; i < (2 * numLevels) - 1; i++)
  {
    svarh = svarh + pow((i+2) - senth, 2) * sum_pixels[i]; 
  }

  float dvarh = 0, denth = 0;
  for (int i = 0; i < numLevels; i++)
  {
    dvarh = dvarh + (pow(i, 2) * dif_pixels[i]);
    denth = denth - dif_pixels[i] * log(dif_pixels[i] + 0.0000000000001); //¹
  }

float energ = 0, contr = 0, homom = 0, homop = 0,  entro = 0, corrm = 0, corrp = 0, autoc = 0,  cprom = 0, cshad = 0, dissi = 0, maxpr = 0,
  sosvh = 0, inf1h = 0, inf2h = 0, indnc = 0,   idmnc = 0;

  float HX = 0, HY = 0, HXY = 0, HXY1 = 0, HXY2 = 0;

  for(int i = 0; i < numLevels; i++)
  {
    HX = HX - px[i] * log(px[i] + 0.0000000000001);
    HY = HY - py[i] * log(py[i] + 0.0000000000001);

    for(int j = 0; j < numLevels; j++)
    {
      autoc = autoc + (i+1) * (j+1) * glcm.at<float>(i,j); 
      contr = contr + (abs(i-j) * abs(i-j) * glcm.at<float>(i,j));
      corrm = corrm + ((i+1) - mu_i) * ((j+1) - mu_j) * glcm.at<float>(i,j); 
      cprom = cprom + pow(((i+1) + (j+1) - mu_i - mu_j), 4) * glcm.at<float>(i,j); 
      cshad = cshad + pow(((i+1) + (j+1) - mu_i - mu_j), 3) * glcm.at<float>(i,j);
      dissi = dissi + (abs(i - j) * glcm.at<float>(i,j));
      energ = energ + glcm.at<float>(i,j) * glcm.at<float>(i,j);

      if(glcm.at<float>(i,j) != 0)
        entro = entro - glcm.at<float>(i,j) * log(glcm.at<float>(i,j));

      homom = homom + glcm.at<float>(i,j) / (1 + abs(i-j));
      homop = homop + (glcm.at<float>(i,j) / (1 + ((i - j) * (i - j))));

      if (glcm.at<float>(i,j) > maxpr)
        maxpr = glcm.at<float>(i,j);

      sosvh = sosvh + (pow((i+1) - mu, 2) * glcm.at<float>(i,j));

      HXY = HXY - glcm.at<float>(i,j) * log(glcm.at<float>(i,j) + 0.0000000000001);
      HXY1 = HXY1 - glcm.at<float>(i,j) * log(px[i] * py[j] + 0.0000000000001);
      HXY2 = HXY2 - px[i] * py[j] * log(px[i] * py[j] + 0.0000000000001);

      indnc = indnc + (glcm.at<float>(i,j) / (1 + (float(abs(i - j)) / numLevels)));
      idmnc = idmnc + (glcm.at<float>(i,j) / (1 + (pow(i - j, 2)) / pow(numLevels, 2)));
    }
  }

  corrp = (autoc - mu_i * mu_j) / (sigma_i * sigma_j);
  corrm = corrm / (sigma_i * sigma_j);

  float valueH = 0;

  if (HX >= HY)
    valueH = HX;
  else
    valueH = HY;

  inf1h = (HXY - HXY1) / valueH;
  inf2h = pow((1 - exp(-2 * (HXY2 - HXY))), 0.5);


haralickFeaturesValues["autoc"]=autoc; haralickFeaturesValues["contr"]=contr; haralickFeaturesValues["corrm"]=corrm; haralickFeaturesValues["corrp"]=corrp;
haralickFeaturesValues["cprom"]=cprom; haralickFeaturesValues["cshad"]=cshad; haralickFeaturesValues["dissi"]=dissi; haralickFeaturesValues["energ"]=energ;
haralickFeaturesValues["entro"]=entro; haralickFeaturesValues["homom"]=homom; haralickFeaturesValues["homop"]=homop; haralickFeaturesValues["maxpr"]=maxpr;
haralickFeaturesValues["sosvh"]=sosvh; haralickFeaturesValues["savgh"]=savgh; haralickFeaturesValues["svarh"]=svarh; haralickFeaturesValues["senth"]=senth;
haralickFeaturesValues["dvarh"]=dvarh; haralickFeaturesValues["denth"]=denth; haralickFeaturesValues["inf1h"]=inf1h; haralickFeaturesValues["inf2h"]=inf2h;
haralickFeaturesValues["indnc"]=indnc; haralickFeaturesValues["idmnc"]=idmnc; 


  return haralickFeaturesValues;
}
