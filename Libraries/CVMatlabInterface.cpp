#include "CVMatlabInterface.h"
#include "glcm.h"

void saveFloatImageAsMatFile(double *image, int arraySize, int dimY, int dimX){
  const char *file = "Image1.mat";
  //Creacion de archivo mat
  MATFile *pmat = matOpen(file, "w");

  int sizeOfImage = dimX*dimY;
  double *textureImg;
  for(int i=0;i < arraySize;i++){
    //Getting the Haralick feature one by one
    textureImg = image + i*sizeOfImage;

    //Creating the matrix to be save in mat file and copying the data from  
    mxArray *pa1 = mxCreateDoubleMatrix(dimY, dimX, mxREAL);
    memcpy(mxGetDoubles(pa1), textureImg, 8*sizeOfImage); 
    // Saving the texture image in the mat file
    const char *name = haralickFeatureNames[i].c_str();
    matPutVariable(pmat, name, pa1);
    mxDestroyArray(pa1);
  }
}
