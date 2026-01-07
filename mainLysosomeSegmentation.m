clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
file.runSegmentation = 'run'; %load or run
file.drawROI = 'off'; %'off', or channel name

MainFolder = {'E:\mini\A549 Lamp1 Gal3'};
HourFolders = {'LysoTracker'};
CellineFolders = {'mSiPEI', 'CQ', 'Control'};

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
                CellInt = [];
                MembrInt = [];
                CellDens = [];
                MembrDens = [];
                VolumeList = [];
                FileNames = [];
                for j = 1:size(SubFolder,1)
                    if isSubDirColumn(j,1) == 1
                        file.path = append(SubFolder(j).folder, filesep, SubFolder(j).name);

                        stack = Core.LysosomeSegmentation(file);
                        stack.loadDataBioform(chan);
                        stack.showChannel;
                        stack.SegmentChannels(Threshold);

                        [Results] = stack.CrossCorrelation;
                    end

                    PairNames = fieldnames(Results);
                    for k = 1:size(PairNames, 1)
                        if j == 1
                            BigResults.(PairNames{k}) = Results.(PairNames{k});
                        else
                            BigResults.(PairNames{k}) = [BigResults.(PairNames{k}); Results.(PairNames{k})];
                        end
                        FileNames{j,1} = file.path;
                    end
                    close all
                end

                fields = fieldnames(BigResults);
                ResultsTables = struct();
                for s = 1:numel(fields)
                    field = fields{s};
                    data = BigResults.(field);   % n × 5 matrix

                    parts = split(field, '_');
                    nameA = string(parts{1});
                    nameB = string(parts{2});

                    varNames = {'FileNames', ...
                        sprintf('PxOverlap fraction_of_%s_containing_%s', nameA, nameB), ...
                        sprintf('PxOverlap fraction_of_%s_containing_%s', nameB, nameA), ...
                        sprintf('Manders fraction_of_%s_containing_%s', nameA, nameB), ...
                        sprintf('Manders fraction_of_%s_containing_%s', nameB, nameA), ...
                        'Pearson' ...
                    };

                    FileNamesPadded = strings(size(data, 1),1);
                    FileNamesPadded(1:size(data, 1)) = string(FileNames);
                    T = table(FileNamesPadded, data(:,1), data(:,2),  data(:,3),  data(:,4), data(:,5),'VariableNames', varNames);
                    % 
                    % T = array2table(data, 'VariableNames', varNames);
                    ResultsTables.(field) = T;

                    writetable(T, append(Path{1}, filesep, 'ResultsCorrelation.xlsx'), 'Sheet', field);
                end
            end
        end
    end
end