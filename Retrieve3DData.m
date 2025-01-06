clear;
close all;
clc;
%% User Input
Path = 'D:\Analysis_done\AuNPs@mSi@PEI\3D\48u\MCF7';
Folder = dir(Path);

SpheroidInt = [];
IntMatrix = [-200:300].';
for j = 3:size(Folder,1)
    if Folder(j).isdir == 1
        try
            file.path = append(Folder(j).folder, filesep, Folder(j).name);
    
            PathToFile = append(file.path, filesep, 'IntensityInFunctionOfDepth.mat');
            IntDepth = load(PathToFile);
            IntDepth = IntDepth.IntensityInFunctionOfDepth;
    
            for l = 1:size(IntDepth, 1)
                pos = IntDepth(l, 1);  
                value = IntDepth(l, 2); 
                rowIndex = find(IntMatrix == pos, 1);  
                if ~isempty(rowIndex)                   
                    IntMatrix(rowIndex, j-1) = value;   
                end
            end
    
            %% Get full spheroid intensity uptake
            PathToFile2 = append(file.path, filesep, 'PartSph.mat');
            PartSph = load(PathToFile2);
            PartSph = PartSph.PartSph;
            Elevation = min(abs(PartSph(:,2)));
            AreaFactor = 2./(1-Elevation);
    
            PathToFile3 = append(file.path, filesep, 'IntListIn.mat');
            IntListIn = load(PathToFile3);
            IntListIn = IntListIn.IntListIn;     
            PartInt = sum(IntListIn(:,1));
            
            TotInt = PartInt*AreaFactor;
        
            SpheroidInt = [SpheroidInt, TotInt];
        catch
        end
    end
end
filename = append(Path, filesep,'IntMatrix.mat');
save(filename, 'IntMatrix');
filename = append(Path, filesep,'SpheroidIntTotal.mat');
save(filename, 'SpheroidInt');
