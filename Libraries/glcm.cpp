#include "glcm.h"
#include "LinkedList.h"

void glcm(cv::Mat img, int numLevels, int pos_y, int pos_x, double imgResult[], int dims[]){
  Node2D *glcm=NULL;  //Head of linked list
  Node2D nodesToBeSaved[(img.rows-1)*img.cols];

  for(int i = 0; i < img.rows-1; i++){
    for(int j = 0; j < img.cols; j++){
      int pos_y = img.at<uint8_t>(i,j);
      int pos_x = img.at<uint8_t>(i+1,j);
      nodesToBeSaved[j+i*img.cols] = Node2D(pos_x, pos_y);
    }
  }

  // Creating GLCM matrix with "numLevels", radius = 1 and in the horizontal direction
  int sizeGLCM = 0; //This will indicates how many pairs with non zero probability the glcm has
  for(int i = 0; i < img.rows-1; i++){
    for(int j = 0; j < img.cols; j++){
      sizeGLCM = sizeGLCM + Node2D::addNewNode(&glcm, &nodesToBeSaved[j+i*img.cols]);
    }      
  }

  // Normalizing GLCM matrix for parameter determination
  glcm->Normalize(img.cols * (img.rows-1));//Should be changed the rescaling factor if we used another direction vector
  
  // Marginal pmf of X, pmf of Y, pmf of x + y, pmf of the |x - y| and expected value of X*Y
  Node1D *p_x=NULL, *p_y=NULL, *p_x_plus_y=NULL, *p_x_minus_y=NULL; //Maybe using a linked list for the marginal probabilities is an overkill
  Node1D p_xArray[sizeGLCM], p_yArray[sizeGLCM], p_x_plus_yArray[sizeGLCM], p_x_minus_yArray[sizeGLCM]; // Static memory allocation
  float mu_x_times_y = 0.0;
  
  //Preloading the Node1D values to the array, so we can create the linked list only pointing to some of them
  Node2D *temp = glcm;
  for(int i=0;i<sizeGLCM;i++){
    float prob = temp->getProb();
    int x = temp->x;
    int y = temp->y;

    p_xArray[i] = Node1D(x, prob);
    p_yArray[i] = Node1D(y, prob);
    p_x_plus_yArray[i] = Node1D(x+y, prob);
    p_x_minus_yArray[i] = Node1D(abs(x-y), prob);
    mu_x_times_y = mu_x_times_y + x * y * prob;

    temp = temp->next;
  }

  for(int i=0;i<sizeGLCM;i++){
    Node1D::addNewNode(&p_x, &p_xArray[i]);
    Node1D::addNewNode(&p_y, &p_yArray[i]);
    Node1D::addNewNode(&p_x_plus_y, &p_x_plus_yArray[i]);
    Node1D::addNewNode(&p_x_minus_y, &p_x_minus_yArray[i]);
  }

  
  // Means and expected values of X and Y squared
  Node1D *aux = p_x;
  float mu_x = 0.0, mu_y = 0.0, mu_x_squared = 0.0, mu_y_squared=0.0;

  while(aux!=NULL){
    float prob = aux->getProb();
    int x = aux->z;
    mu_x = mu_x + x * prob;
    mu_x_squared = mu_x_squared + x * x * prob;

    aux = aux->next;
  }

  aux = p_y;
  while(aux!=NULL){
    float prob = aux->getProb();
    int y = aux->z;
    mu_y = mu_y + y * prob;
    mu_y_squared = mu_y_squared + y * y * prob;

    aux = aux->next;
  }


  // Expected value of X+Y and (X+Y)^2
  float mu_x_plus_y = 0.0, mu_x_plus_y_squared = 0.0;
  aux = p_x_plus_y;
  while(aux!=NULL){
    float prob = aux->getProb();
    int z = aux->z;
    mu_x_plus_y = mu_x_plus_y + z * prob;
    mu_x_plus_y_squared = mu_x_plus_y_squared + z * z * prob;
    
    aux = aux->next;
  }


  // Expected value of |X-Y| and (X-Y)^2
  float mu_x_minus_y = 0.0, mu_x_minus_y_squared = 0.0;
  aux = p_x_minus_y;
  while(aux!=NULL){
    float prob = aux->getProb();
    int z = aux->z;
    mu_x_minus_y = mu_x_minus_y + z * prob;
    mu_x_minus_y_squared = mu_x_minus_y_squared + z * z * prob;
    
    aux = aux->next;
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

  temp = glcm;
  while(temp!=NULL){
    float prob = temp->getProb();
    int x = temp->x;
    int y = temp->y;
    float inverse_moment_value_norm = 1 + (x - y)*(x - y)/(float)(numLevels * numLevels);
    float inverse_diff_value_norm = 1 + ((float)abs(x - y))/((float)numLevels);

    energy = energy + prob * prob;
    inverse_diff_moment_norm = inverse_diff_moment_norm + prob/inverse_moment_value_norm;
    inverse_diff_norm = inverse_diff_norm + prob/inverse_diff_value_norm;
    entropy = entropy - prob * log(prob + 0000000000001);
    if (prob > maximum_prob)  maximum_prob = prob;

    temp = temp->next;
  }


  //Entropy of X+Y, entropy of |X-Y|, third and fourth central moment of X+Y aka cluster shade and cluster prominence respectively
  float entropy_x_plus_y = 0.0, entropy_x_minus_y=0.0, cluster_shade = 0.0, cluster_prominence = 0.0;

  aux = p_x_plus_y;
  while(aux!=NULL){
    float prob = aux->getProb();
    int z = aux->z;
    entropy_x_plus_y = entropy_x_plus_y - prob *  log(prob + 0.0000000000001);
    cluster_shade = cluster_shade + pow(z - mu_x_plus_y, 3) * prob;
    cluster_prominence = cluster_prominence + pow(z - mu_x_plus_y, 4) * prob;
    
    aux = aux->next;
  }

  aux = p_x_minus_y;
  while(aux!=NULL){
    float prob = aux->getProb();
    entropy_x_minus_y = entropy_x_minus_y - prob *  log(prob + 0.0000000000001);
    
    aux = aux->next;
  }


  // Entropy of X, entropy of Y
  float entropy_x = 0.0, entropy_y = 0.0;
  aux = p_x;
  while(aux!=NULL){
    float prob_x = aux->getProb();
    entropy_x = entropy_x - prob_x *  log(prob_x + 0.0000000000001);
    
    aux = aux->next;
  }

  aux = p_y;
  while(aux!=NULL){
    float prob_y = aux->getProb();
    entropy_y = entropy_y - prob_y *  log(prob_y + 0.0000000000001);
    
    aux = aux->next;
  }


  // Getting the coefficients used in the information measures of correlation
  float coeff_XY_1 = 0.0, coeff_XY_2 = 0.0;
  Node1D *aux2 = p_y;
  while(aux2!=NULL){
    aux = p_x;
    while(aux!=NULL){
      float prob_x = aux->getProb();
      float prob_y = aux2->getProb();
      float prob_x_y = Node2D::getProb(glcm, aux->z, aux2->z);

      coeff_XY_1 = coeff_XY_1 - prob_x_y * log(prob_x*prob_y + 0.0000000000001);
      coeff_XY_2 = coeff_XY_2 - prob_x * prob_y * log(prob_x*prob_y + 0.0000000000001);

      aux = aux->next;
    }
    aux2 = aux2->next;
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

