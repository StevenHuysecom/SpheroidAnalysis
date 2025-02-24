clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
MainFolder = {'D:\Steven\Au@mSi'};
DimensionFolders = {'3D'};
HourFolders = {'48hour'};
ParticleFolders = {'HeLa'};

%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.ch01 = 'Membrane';
chan.ch02 = 'Particles';
chan.ch03 = 'ignore';
chan.ch04 = 'ignore';

for m = 1:numel(DimensionFolders)
    DimensionFolder = DimensionFolders{m};
    for a = 1:numel(HourFolders)
        HourFolder = HourFolders{a};
        for r = 1:numel(ParticleFolders)
            try
                ParticleFolder = ParticleFolders{r};
                Path = append(MainFolder, filesep, DimensionFolder, filesep, HourFolder,...
                    filesep, ParticleFolder);
                % Path = append(MainFolder, filesep, ParticleFolder);
                file.path = Path{1,1};
    
                Load.Movie.lif.LoadImages(file, chan);
                CurrentFolder = dir(file.path);
                CurrentFolder(1:2) = [];
                isDirColumn = [CurrentFolder.isdir]';
    
                for i = 1:size(CurrentFolder,1)
                    try
                        if isDirColumn(i,1) == 1
                            SubFolder = dir(append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name));
                            SubFolder(1:2) = [];
                            isSubDirColumn = [SubFolder.isdir]';
                            SpheroidInt = [];
                            IntMatrix = [-100:300].';
                            for j = 1:size(SubFolder,1)
                                try
                                    if isSubDirColumn(j,1) == 1
                                        file.path = append(SubFolder(j).folder, filesep, SubFolder(j).name);
            
                                        PxSize = load(append(file.path, filesep, 'PxSizes.mat'));
                                        info.pxSizeXY = PxSize.PxSizes(1);
                                        info.pxSizeZ  = PxSize.PxSizes(3);
                                        stack = Core.Spheroid3D(file,info);
                                        stack.loadDataBioform(chan);
                                        stack.showChannel;
                                        
                                        %% Find center of spheroid
                                        stack.findCenter;
                                        
                                        %% Change to ellipsoid coordinates
                                        stack.CartToEllipsoid;
                        
                                        %% Integrate over r
                                        [IntDepth] = stack.IntegrateR;
                                        for l = 1:size(IntDepth, 1)
                                            pos = IntDepth(l, 1);  
                                            value = IntDepth(l, 2); 
                                            rowIndex = find(IntMatrix == pos, 1);  
                                            if ~isempty(rowIndex)                   
                                                IntMatrix(rowIndex, j+1) = value;   
                                            end
                                        end
            
                                        %% Get full spheroid intensity uptake
                                        [TotInt] = stack.GetFullInt;
                                        SpheroidInt = [SpheroidInt, TotInt];
                                    end
                                catch 
                                end
                            end
                            filename = append(SubFolder(1).folder, filesep,'IntMatrix.mat');
                            save(filename, 'IntMatrix');
                            filename = append(SubFolder(1).folder, filesep,'SpheroidIntTotal.mat');
                            save(filename, 'SpheroidInt');
                        end
                    catch
                    end
                end
            catch
            end
        end
    end
end