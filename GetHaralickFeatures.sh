#!/bin/bash

images=(`ls Images`)
rm Outputs/*

for ((i=0; i<${#images[@]}; i++)); do
#do something to each element of array
  ./glcm "Images/${images[$i]}" 8 >> "$(echo "Outputs/${images[$i]}" | cut -f 1 -d '.').csv"
  #sleep 1
done
