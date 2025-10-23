clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
file.runSegmentation = 'run'; %load or run

% info.pxSizeXY = 454.5; 
% info.pxSizeZ  = 1000;
info.Membrane = 'excluded'; %included or excluded

MainFolder = {'D:'};
HourFolders = {'mini'};%'3hour', '24hour', '48hour'
CellineFolders = {'lamp1 confocal test for steven'};

%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.ch01 = 'Lysotracker';
chan.ch02 = 'IetAnders';
chan.ch03 = 'Particles';
chan.ch04 = 'ignore';

%% Loading data

for a = 1:numel(HourFolders)
    HourFolder = HourFolders{a};
    for r = 1:numel(CellineFolders)
            CellineFolder = CellineFolders{r};
            Path = append(MainFolder, filesep, HourFolder, filesep, CellineFolder);
            file.path = Path{1,1};

            Load.Movie.lif.LoadImages(file, chan);

            CurrentFolder = dir(file.path);
            CurrentFolder(1:2) = [];
            isDirColumn = [CurrentFolder.isdir]';
            for i = 1:size(CurrentFolder,1)
                if isDirColumn(i,1) == 1
                        SubFolder = dir(append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name));
                        SubFolder(1:2) = [];
                        isSubDirColumn = [SubFolder.isdir]';
                        CellInt = [];
                        MembrInt = [];
                        CellDens = [];
                        MembrDens = [];
                        VolumeList = [];
                        for j = 1:size(SubFolder,1)
                                if isSubDirColumn(j,1) == 1
                                    file.path = append(SubFolder(j).folder, filesep, SubFolder(j).name);

                                    stack = Core.LysosomeSegmentation(file,info);
                                    stack.loadDataBioform(chan);
                                    stack.showChannel;
                                end
                            close all
                        end
                end
            end
    end
end