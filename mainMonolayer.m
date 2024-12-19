clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
file.runSegmentation = 'load'; %load or run

% info.pxSizeXY = 454.5; 
% info.pxSizeZ  = 1000;
info.Membrane = 'excluded'; %included or excluded

MainFolder = {'D:\Data Uptake\AuNP@mSi@PEI'};
DimensionFolders = {'2D'};
HourFolders = {'48hour'};
CellineFolders = {'KM12C'};

%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.ch01 = 'Membrane';
chan.ch02 = 'Particles';
chan.ch03 = 'ignore';
chan.ch04 = 'ignore';

%% Loading data

for m = 1:numel(DimensionFolders)
    DimensionFolder = DimensionFolders{m};
    for a = 1:numel(HourFolders)
        HourFolder = HourFolders{a};
        for r = 1:numel(CellineFolders)
            try
                CellineFolder = CellineFolders{r};
                Path = append(MainFolder, filesep, DimensionFolder, filesep, HourFolder,...
                    filesep, CellineFolder);
                file.path = Path{1,1};
    
                Load.Movie.lif.LoadImages(file, chan);
    
                CurrentFolder = dir(file.path);
                CurrentFolder(1:2) = [];
                isDirColumn = [CurrentFolder.isdir]';
                for i = 1:size(CurrentFolder,1)
                    if isDirColumn(i,1) == 1
                        try
                            SubFolder = dir(append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name));
                            SubFolder(1:2) = [];
                            isSubDirColumn = [SubFolder.isdir]';
                            CellInt = [];
                            MembrInt = [];
                            for j = 1:size(SubFolder,1)
                                try
                                    if isSubDirColumn(j,1) == 1
                                        file.path = append(SubFolder(j).folder, filesep, SubFolder(j).name);
                        
                                        stack = Core.MonolayerSegmentation(file,info);
                                        stack.loadDataBioform(chan);
                                        stack.showChannel;
                        
                                        %% DL segmentation per plane
                                        Ws = stack.segmentMembraneDL;
                        
                                        %% Particle Intensity
                                        [Int, MembraneInt]= stack.getIntensityInside();
                                        CellInt = [CellInt; Int];
                                        MembrInt = [MembrInt; MembraneInt];
                                    end
                                catch
                                    continue
                                end
                            end
                            filenameCell = append(SubFolder(1).folder, filesep,'CellInt.mat');
                            filenameMembr = append(SubFolder(1).folder, filesep,'MembrInt.mat');
                            save(filenameCell, 'CellInt');
                            save(filenameMembr, 'MembrInt');
                            close all
                        catch
                            continue
                        end
                    end
                end
            catch
                continue
            end
        end
    end
end