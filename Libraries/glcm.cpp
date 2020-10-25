#include "glcm.h"


void glcm(Mat img, int numLevels, int pos_y, int pos_x, double *imgResult, int dims[]){
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
  glcm = glcm/(img.cols * (img.rows-1));//Should be changed the rescaling factor if we used another direction vector

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
  // Energy, inverse difference moment normalized, inverse difference normalized, entropy, Maximum probability
  // The unnormalized inverse different moment is also called homogeneity, so I going to use only the normalized version
  float energy = 0.0, inverse_diff_moment_norm = 0.0, inverse_diff_norm = 0.0, entropy = 0.0, maximum_prob = 0.0;
  for(int y = 0; y < numLevels; y++){
    for(int x = 0; x < numLevels; x++){
      float prob =  glcm.at<float>(y,x);
      float inverse_moment_value_norm = 1 + (x - y)*(x - y)/(float)(numLevels * numLevels);
      float inverse_diff_value_norm = 1 + ((float)abs(x - y))/((float)numLevels) ;

      energy = energy + prob * prob;
      inverse_diff_moment_norm = inverse_diff_moment_norm + prob/inverse_moment_value_norm;
      inverse_diff_norm = inverse_diff_norm + prob/inverse_diff_value_norm;
      entropy = entropy - prob * log(prob + 0000000000001);
      if (prob > maximum_prob)  maximum_prob = prob;
    }
  }

  //Entropy of X+Y, entropy of |X-Y|, third and fourth central moment of X+Y aka cluster shade and cluster prominence respectively
  float entropy_x_plus_y = 0.0, entropy_x_minus_y=0.0, cluster_shade = 0.0, cluster_prominence = 0.0;
  for (int z=0; z <2*numLevels-1;z++){
    float prob = p_x_plus_y[z];
    entropy_x_plus_y = entropy_x_plus_y - prob *  log(prob + 0.0000000000001);
    cluster_shade = cluster_shade + pow(z - mu_x_plus_y, 3) * prob;
    cluster_prominence = cluster_prominence + pow(z - mu_x_plus_y, 4) * prob;
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

  // We don't use the two versions of homogeneity because we use the normalized version that is inverse difference moment normalized and inverse difference normalized

  // Now we save the haralick features instead of returning them
  int dim_y = dims[1], dim_x = dims[2];
  int offset = pos_y + pos_x * dim_y; int sizeTextureImg = dim_x * dim_y; // Saving in col-major order
  //These Haralick features are saved in the same order as the constant haralickFeatureNames in glcm.h
  imgResult[offset + 0 * sizeTextureImg] = mu_x_times_y; imgResult[offset + 1 * sizeTextureImg] = mu_x_minus_y_squared; imgResult[offset + 2 * sizeTextureImg] = corr;
  imgResult[offset + 3 * sizeTextureImg] = cluster_prominence; imgResult[offset + 4 * sizeTextureImg] = cluster_shade; imgResult[offset + 5 * sizeTextureImg] = mu_x_minus_y;
  imgResult[offset + 6 * sizeTextureImg] = energy; imgResult[offset + 7 * sizeTextureImg] = entropy; imgResult[offset + 8 * sizeTextureImg] = maximum_prob;
  imgResult[offset + 9 * sizeTextureImg] = mu_x_plus_y; imgResult[offset + 10 * sizeTextureImg] = var_x_plus_y; imgResult[offset + 11 * sizeTextureImg] = entropy_x_plus_y;

  imgResult[offset + 12 * sizeTextureImg] = var_x_minus_y; imgResult[offset + 13 * sizeTextureImg]= entropy_x_minus_y; imgResult[offset + 14 * sizeTextureImg] = information_measure_corr_1;
  imgResult[offset + 15 * sizeTextureImg] = information_measure_corr_2; imgResult[offset + 16 * sizeTextureImg] = inverse_diff_norm; imgResult[offset + 17 * sizeTextureImg] = inverse_diff_moment_norm;
  imgResult[offset + 18 * sizeTextureImg] = var_x; imgResult[offset + 19 * sizeTextureImg] = var_y; 
}

