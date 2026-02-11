clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
file.runSegmentation = 'run'; %load or run
file.drawROI = 'off'; %'off', or channel name

MainFolder = {'D:\mini\A549 Lamp1 Gal3'};
HourFolders = {'LysoTracker'};
CellineFolders = {'mSiPEI'};

%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.ch01 = 'Lamp1';
chan.ch02 = 'Lysotracker';
chan.ch03 = 'Galactin';
chan.ch04 = 'ignore';

%some parameters
slice = 1; %which slice of the 3D stack to select the ROI on
Threshold = [0.10, 0.10, 0.10]; %[0-1], keep it under 0.15, intensity threshold for lysosomes (high = throw away dim/out-of-focus lysosomes)

%% Loading data

for a = 1:numel(HourFolders)
    HourFolder = HourFolders{a};
    for r = 1:numel(CellineFolders)
        CellineFolder = CellineFolders{r};
        Path = append(MainFolder, filesep, HourFolder, filesep, CellineFolder);
        file.path = Path{1,1};

        Load.Movie.lif.LoadImages(file, chan);
        Load.Movie.DrawApplyROI(file, chan, slice);
    end
end


for a = 1:numel(HourFolders)
    HourFolder = HourFolders{a};
    for r = 1:numel(CellineFolders)
        CellineFolder = CellineFolders{r};
        Path = append(MainFolder, filesep, HourFolder, filesep, CellineFolder);
        file.path = Path{1,1};
        CurrentFolder = dir(file.path);
        CurrentFolder(1:2) = [];
        isDirColumn = [CurrentFolder.isdir]';
        for i = 1:size(CurrentFolder,1)
            if isDirColumn(i,1) == 1
                SubFolder = dir(append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name));
                SubFolder(1:2) = [];
                isSubDirColumn = [SubFolder.isdir]';

                BigResults = [];
                FileNames = [];

                for j = 1:size(SubFolder,1)
                    if isSubDirColumn(j,1) == 1
                        file.path = append(SubFolder(j).folder, filesep, SubFolder(j).name);

                        stack = Core.LysosomeSegmentation(file);
                        stack.loadDataBioform(chan);
                        stack.showChannel;
                        stack.SegmentChannels(Threshold);

                        [Results] = stack.GalectinAnalysis;
                    end

                    BigResults = [BigResults; Results*100];

                    FileNames{j,1} = file.path;
                    close all
                end

                varNames = {'FileNames', CellineFolder};

                FileNamesPadded = strings(size(BigResults, 1),1);
                FileNamesPadded(1:size(BigResults, 1)) = string(FileNames);
                T = table(FileNamesPadded, BigResults ,'VariableNames', varNames);

            end
        end
        writetable(T, append(MainFolder{1}, filesep, HourFolder, filesep, 'GalactinAnalysis.xlsx'), 'Sheet', CellineFolder);
    end
end