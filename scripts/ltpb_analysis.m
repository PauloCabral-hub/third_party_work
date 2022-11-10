% LTPB collaboration analysis.

addpath(genpath('/home/roberto/Documents/Pos-Doc/BagOfCodesAndData/Git')) 
T = readtable('/home/roberto/Documents/Pos-Doc/BagOfCodesAndData/ExperimentData/ColetaLTPB/DataframeFull_20102022.csv');


OriData = zeros(size(T,1),size(T,2));
for k = 1:size(T,2)
    for a = 1:size(T,1)
        if k < 12
           OriData(a,k) = T{a,k};
        else
           if strcmp(T{a,12}{1,1},'LTPB')   % LTPB will be labeled 1
              OriData(a,k) = 1; 
           else                             % Control will be labed 2
              OriData(a,k) = 2;   
           end
           if strcmp(T{a,13}{1,1},'Online') % Online will be labeled 1
              OriData(a,k) = 1; 
           else                             % Presential will be labed 2
              OriData(a,k) = 2;   
           end           
        end
    end
end


