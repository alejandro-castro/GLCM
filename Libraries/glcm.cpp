#include "glcm.h"


map<String, float>  glcm(Mat &img, int numLevels){
  map <String, float> haralickFeaturesValues;
  Mat glcm = Mat::zeros(numLevels, numLevels, CV_32F);
  
  // Creating GLCM matrix with "numLevels", radius = 1 and in the horizontal direction
  for(int i = 0; i < img.rows-1; i++){
    for(int j = 0; j < img.cols; j++){
      glcm.at<float>(img.at<uchar>(i,j) , img.at<uchar>(i+1,j) ) += 1; //glcm.at<float>(img.at<uchar>(i,j) - 1, img.at<uchar>(i,j+1) - 1) += 1;
    }      
  }

  // Normalizing GLCM matrix for parameter determination
  glcm = glcm / sum(glcm)[0];

  // Means
  float mu_i = 0, mu_j = 0;
  float mu = 0;
  for(int i = 0; i < numLevels; i++){
    for(int j = 0; j < numLevels; j++){
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

  for(int i = 0; i < numLevels; i++){
    for(int j = 0; j < numLevels; j++){
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
  for(int i = 0; i < (2 * numLevels) - 1; i++){
    savgh = savgh + ((i+2) * sum_pixels[i]); 
    senth = senth - (sum_pixels[i] * log(sum_pixels[i] + 0.0000000000001)); //¹
  }

  float svarh = 0;
  for(int i = 0; i < (2 * numLevels) - 1; i++){
    svarh = svarh + pow((i+2) - senth, 2) * sum_pixels[i]; 
  }

  float dvarh = 0, denth = 0;
  for (int i = 0; i < numLevels; i++){
    dvarh = dvarh + (pow(i, 2) * dif_pixels[i]);
    denth = denth - dif_pixels[i] * log(dif_pixels[i] + 0.0000000000001); //¹
  }
float corrm = 0;
float energ = 0, contr = 0, homom = 0, homop = 0,  entro = 0, corrp = 0, autoc = 0,  cprom = 0, cshad = 0, dissi = 0, maxpr = 0,
  sosvh = 0, inf1h = 0, inf2h = 0, indnc = 0,   idmnc = 0;

  float HX = 0, HY = 0, HXY = 0, HXY1 = 0, HXY2 = 0;

  for(int i = 0; i < numLevels; i++){
    HX = HX - px[i] * log(px[i] + 0.0000000000001);
    HY = HY - py[i] * log(py[i] + 0.0000000000001);

    for(int j = 0; j < numLevels; j++){
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

haralickFeaturesValues["autoc"]=autoc; 
haralickFeaturesValues["contr"]=contr; 
haralickFeaturesValues["corrm"]=corrm; haralickFeaturesValues["corrp"]=corrp;
haralickFeaturesValues["cprom"]=cprom; haralickFeaturesValues["cshad"]=cshad; haralickFeaturesValues["dissi"]=dissi; haralickFeaturesValues["energ"]=energ;
haralickFeaturesValues["entro"]=entro; haralickFeaturesValues["homom"]=homom; haralickFeaturesValues["homop"]=homop; haralickFeaturesValues["maxpr"]=maxpr;
haralickFeaturesValues["sosvh"]=sosvh; haralickFeaturesValues["savgh"]=savgh; haralickFeaturesValues["svarh"]=svarh; haralickFeaturesValues["senth"]=senth;
haralickFeaturesValues["dvarh"]=dvarh; haralickFeaturesValues["denth"]=denth; haralickFeaturesValues["inf1h"]=inf1h; haralickFeaturesValues["inf2h"]=inf2h;
haralickFeaturesValues["indnc"]=indnc; haralickFeaturesValues["idmnc"]=idmnc; 
return haralickFeaturesValues;

}



void FillPixelWithHaralickFeatures(Mat img[], map <string, float> haralickFeatures, int y, int x){
  for (int i=0;i<numberOfFeaturesToBeUsed;i++){
    Mat textureImg = img[i];
    string name = haralickFeatureNames[i];
    textureImg.at<double>(y, x)  = haralickFeatures[name];
  }
}

void initializeArrayOfImages(Mat img[], int arraySize, int dimY, int dimX){
  for(int i=0;i<arraySize; i++){
    img[i] = Mat::zeros(dimY, dimX, CV_64F); 
  }
}