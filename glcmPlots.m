clear; close all; clc;
sizes = [13, 15];

qualities = ["LowQuality", "HighQuality"];
features = containers.Map({'LowQuality', 'HighQuality'}, {containers.Map(), containers.Map()});
for quality = qualities
    qualityMap = features(quality);
    for i = 1:13
        data = readtable('Outputs/'+quality+string(i)+'.csv');

        for j = 1:22
           key = string(data(j,1).Variables);
           value = data(j,2).Variables;
           if isKey(qualityMap, key)
                qualityMap(key) = [qualityMap(key), value];
           else
                qualityMap(key) = [value];
           end       
        end
    end
end

for key = features("LowQuality").keys()
    figure;
    hold on;
    data = [];
    qualityIndexes = [];
    for quality = qualities
        qualityMap = features(quality); 
        data = [data, qualityMap(string(key))];
        nrep = length(qualityMap(string(key)));
        qualityIndexes = [qualityIndexes, repmat(quality,1,nrep)];
    end
    boxplot(data, qualityIndexes)
    pause(2);       

    %xlabel(quality)
    ylabel(key)
    title('Miles per Gallon for All Vehicles')    
end
    
