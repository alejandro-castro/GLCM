#include "CVMatlabInterface.h"
#include "glcm.h"

void saveFloatImageAsMatFile(Mat image[], int arraySize, int dimY, int dimX){
  const char *file = "Image1.mat";
  
  //Creacion de archivo mat
  MATFile *pmat = matOpen(file, "w");
  mxArray *pa1;

  for(int i=0;i < arraySize;i++){
    //Getting the Haralick feature one by one
    Mat textureImg = image[i];
    
    //Creating the matrix to be save in mat file and copying the data from  
    pa1 = mxCreateDoubleMatrix(dimY, dimX, mxREAL);
    if (textureImg.isContinuous()){
        memcpy((void *)(mxGetPr(pa1)), (void *)(textureImg.data), textureImg.elemSize()*dimX*dimY); 
    }    
    else{
        cout <<"Imagen en memoria no continua, no se puede copiar directamente la data";
    }
    // Saving the texture image in the mat file
    const char *name = haralickFeatureNames[i].c_str();
    matPutVariable(pmat, name, pa1);
  }
  mxDestroyArray(pa1);
  matClose(pmat);
}
