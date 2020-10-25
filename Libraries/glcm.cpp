#include "glcm.h"


void glcm(Mat img, int numLevels, int pos_y, int pos_x, Mat imgResult[], bool print){
  Mat glcm = Mat::zeros(numLevels, numLevels, CV_32F);
  
  // Creating GLCM matrix with "numLevels", radius = 1 and in the horizontal direction
  for(int i = 0; i < img.rows-1; i++){
    for(int j = 0; j < img.cols; j++){
      int pos_y = img.at<uint8_t>(i,j);
      int pos_x = img.at<uint8_t>(i+1,j);
      glcm.at<float>(pos_y, pos_x) = glcm.at<float>(pos_y, pos_x) + 1; 
    }      
  }

  // Normalizing GLCM matrix for parameter determination
  glcm = glcm / sum(glcm)[0];

  // Marginal pmf of X, pmf of Y, pmf of x + y, pmf of the |x - y| and expected value of X*Y
  vector <float> p_x(numLevels, 0.0), p_y(numLevels, 0.0), p_x_plus_y(2*numLevels-1, 0), p_x_minus_y(numLevels, 0.0);
  float mu_x_times_y = 0.0;
  for(int y = 0; y < numLevels; y++){
    for(int x = 0; x < numLevels; x++){
      float prob = glcm.at<float>(y,x);
      p_y[y] = p_y[y] + prob;
      p_x[x] = p_x[x] + prob;
      p_x_plus_y[x+y] = p_x_plus_y[x+y] + prob;
      p_x_minus_y[abs(x-y)] = p_x_minus_y[abs(x-y)] + prob;
      mu_x_times_y = mu_x_times_y + x * y * prob;
    }
  }

  // Means and expected values of X and Y squared
  float mu_x = 0.0, mu_y = 0.0, mu_x_squared = 0.0, mu_y_squared=0.0;
  for(int x=0; x < numLevels; x++){
    float prob = p_x[x];
    mu_x = mu_x + x * prob;
    mu_x_squared = mu_x_squared + x * x * prob;
  }
  for(int y=0; y < numLevels; y++){
    float prob = p_y[y];
    mu_y = mu_y + y * prob;
    mu_y_squared = mu_y_squared + y * y * prob;
  }

  // Expected value of X+Y and (X+Y)^2
  float mu_x_plus_y = 0.0, mu_x_plus_y_squared = 0.0;
  for(int z=0; z < 2*numLevels-1; z++){
    float prob = p_x_plus_y[z];
    mu_x_plus_y = mu_x_plus_y + z * prob;
    mu_x_plus_y_squared = mu_x_plus_y_squared + z * z * prob;
  }

  // Expected value of |X-Y| and (X-Y)^2
  float mu_x_minus_y = 0.0, mu_x_minus_y_squared = 0.0;
  for(int z=0; z < numLevels; z++){
    float prob = p_x_minus_y[z];
    mu_x_minus_y = mu_x_minus_y + z * prob;
    mu_x_minus_y_squared = mu_x_minus_y_squared + z * z * prob;
  }

  // Variances and correlation
  float var_x = mu_x_squared - mu_x * mu_x; float var_y = mu_y_squared - mu_y * mu_y;
  float corr = (mu_x_times_y - mu_x * mu_y)/sqrt(var_x * var_y);

  // Variance of X+Y, variance of |X-Y|
  float var_x_plus_y = mu_x_plus_y_squared - mu_x_plus_y * mu_x_plus_y;
  float var_x_minus_y = mu_x_minus_y_squared - mu_x_minus_y * mu_x_minus_y;

  // SE PUEDE HACER TODO EN UN SOLO BUCLE PARA EVITAR ACCESOS ADICIONALES A MEMORIA
  // Energy, inverse difference moment normalized, inverse difference normalized, entropy
  // The unnormalized inverse different moment is also called homogeneity, so I going to use only the normalized version
  float energy = 0.0, inverse_diff_moment_norm = 0.0, inverse_diff_norm = 0.0, entropy = 0.0;
  for(int y = 0; y < numLevels; y++){
    for(int x = 0; x < numLevels; x++){
      float prob =  glcm.at<float>(y,x);
      float inverse_moment_value_norm = 1 + (x - y)*(x - y)/(float)(numLevels * numLevels);
      float inverse_diff_value_norm = 1 + ((float)abs(x - y))/((float)numLevels) ;

      energy = energy + prob * prob;
      inverse_diff_moment_norm = inverse_diff_moment_norm + prob/inverse_moment_value_norm;
      inverse_diff_norm = inverse_diff_norm + prob/inverse_diff_value_norm;
      entropy = entropy - prob * log(prob + 0000000000001);
    }
  }

  //Entropy of X+Y, entropy of |X-Y|, se puede juntar con el bucle de x+y para reducir accesos a memoria
  float entropy_x_plus_y = 0.0, entropy_x_minus_y=0.0;
  for (int z=0; z <2*numLevels-1;z++){
    float prob = p_x_plus_y[z];
    entropy_x_plus_y = entropy_x_plus_y - prob *  log(prob + 0.0000000000001);
  }
  for (int z=0; z <numLevels;z++){
    float prob = p_x_minus_y[z];
    entropy_x_minus_y = entropy_x_minus_y - prob *  log(prob + 0.0000000000001);
  }

  // Entropy of X, entropy of Y
  float entropy_x = 0.0, entropy_y = 0.0;
  for (int i=0; i <numLevels;i++){
    float prob_x = p_x[i], prob_y = p_y[i];
    entropy_x = entropy_x - prob_x *  log(prob_x + 0.0000000000001);
    entropy_y = entropy_y - prob_y *  log(prob_y + 0.0000000000001);
  }

  // Getting the coefficients used in the information measures of correlation
  float coeff_XY_1 = 0.0, coeff_XY_2 = 0.0;
  for (int y=0; y< numLevels; y++){
    for (int x=0; x <numLevels;x++){
      float prob_x = p_x[x], prob_y = p_y[y];
      coeff_XY_1 = coeff_XY_1 - glcm.at<float>(y,x) * log(prob_x*prob_y + 0.0000000000001);
      coeff_XY_2 = coeff_XY_2 - prob_x * prob_y * log(prob_x*prob_y + 0.0000000000001);
    }
  }

  // Information of correlation measures
  float information_measure_corr_1 = (entropy - coeff_XY_1)/max(entropy_x, entropy_y);
  float information_measure_corr_2 = sqrt(1 - exp(-2*(coeff_XY_2 - entropy) ));

  // Third and Fourth central moment of X+Y, aka cluster shade and cluster prominence respectively
  float cluster_shade = 0.0, cluster_prominence = 0.0;
  for(int z=0; z < 2*numLevels-1; z++){
    float prob = p_x_plus_y[z];
    cluster_shade = cluster_shade + pow(z - mu_x_plus_y, 3) * prob;
    cluster_prominence = cluster_prominence + pow(z - mu_x_plus_y, 4) * prob;
  }


  // Maximum probability, AGAIN THIS COULD BE MERGED WITH ANOTHER BUCLE TO DECREASE NUMBER OF MEMORY ACCESSES
  float maximum_prob = 0.0;
  for(int y = 0; y < numLevels; y++){
    for(int x = 0; x < numLevels; x++){
      float prob = glcm.at<float>(y,x);
      if (prob > maximum_prob)  maximum_prob = prob;
    }
  }
  //if (print) cout << var_x<<" "<< var_y<< " "<<mu_x << " "<<mu_y<<" "<<mu_x_times_y<<" "<<corr<<endl;
  // We don't use the two versions of homogeneity because we use the normalized version that is inverse difference moment normalized and inverse difference normalized

  // Now we save the haralick features instead of returning them
  imgResult[0].at<double>(pos_y, pos_x) = mu_x_times_y; imgResult[1].at<double>(pos_y, pos_x) = mu_x_minus_y_squared; imgResult[2].at<double>(pos_y, pos_x) = corr;
  imgResult[3].at<double>(pos_y, pos_x) = cluster_prominence; imgResult[4].at<double>(pos_y, pos_x) = cluster_shade; imgResult[5].at<double>(pos_y, pos_x) = mu_x_minus_y;
  imgResult[6].at<double>(pos_y, pos_x) = energy; imgResult[7].at<double>(pos_y, pos_x) = entropy; imgResult[8].at<double>(pos_y, pos_x) = maximum_prob;
  imgResult[9].at<double>(pos_y, pos_x) = mu_x_plus_y; imgResult[10].at<double>(pos_y, pos_x) = var_x_plus_y; imgResult[11].at<double>(pos_y, pos_x) = entropy_x_plus_y;

  imgResult[12].at<double>(pos_y, pos_x) = var_x_minus_y; imgResult[13].at<double>(pos_y, pos_x) = entropy_x_minus_y; imgResult[14].at<double>(pos_y, pos_x) = information_measure_corr_1;
  imgResult[15].at<double>(pos_y, pos_x) = information_measure_corr_2; imgResult[16].at<double>(pos_y, pos_x) = inverse_diff_norm; imgResult[17].at<double>(pos_y, pos_x) = inverse_diff_moment_norm;
  imgResult[18].at<double>(pos_y, pos_x) = var_x; imgResult[19].at<double>(pos_y, pos_x) = var_y; 
}



void FillPixelWithHaralickFeatures(Mat img[], map <string, float> haralickFeatures, int y, int x, bool print){
  for (int i=0;i<numberOfFeaturesToBeUsed;i++){
    Mat textureImg = img[i];
    string name = haralickFeatureNames[i];
    textureImg.at<double>(y, x)  = haralickFeatures[name];
    if ((print)&&(i==2)){
      cout <<"Correlation: "<< haralickFeatures[name]<<endl;
      
    }
    //cout <<"Feature "<<i<<": "<<&(textureImg.data)<<endl;
  }
}

void initializeArrayOfImages(Mat img[], int arraySize, int dimY, int dimX){
  for(int i=0;i<arraySize; i++){
    img[i] = Mat::zeros(dimY, dimX, CV_64F); 

    //cout << img[i];
  }
}