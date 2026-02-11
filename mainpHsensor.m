clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
file.runSegmentation = 'load'; %load or run
file.drawROI = 'off'; %'off', or channel name

MainFolder = {'D:\mini\FIREphly sensor'};
HourFolders = {'20260204', '20260115'};
CellineFolders = {'mSi', 'mSiPEI', 'Control', 'Chloroquine'}; %'mSi', 'mSiPEI', 'Control', , 'Chloroquine'

%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.ch01 = 'mTFP1';
chan.ch02 = 'mCherry';
chan.ch03 = 'Particles';
chan.ch04 = 'Lysotracker';

%some parameters
slice = 1; %which slice of the 3D stack to select the ROI on
Threshold = [0.01, 0.1, 0.10, 0.20]; %[0-1], keep it under 0.15, intensity threshold for lysosomes (high = throw away dim/out-of-focus lysosomes)

%% Loading data

for a = 1:numel(HourFolders)
    HourFolder = HourFolders{a};
    for r = 1:numel(CellineFolders)
        try
            CellineFolder = CellineFolders{r};
            Path = append(MainFolder, filesep, HourFolder, filesep, CellineFolder);
            file.path = Path{1,1};
    
            Load.Movie.lif.LoadImages(file, chan);
            Load.Movie.DrawApplyROI(file, chan, slice);
        catch
        end
    end
end


for a = 1:numel(HourFolders)
    HourFolder = HourFolders{a};
    for r = 1:numel(CellineFolders)
        try
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
    
    
                    Results  = [];
                    Filename = {};
                    k = 0;   % index for valid entries
                    
                    for j = 1:size(SubFolder,1)
                    
                        if isSubDirColumn(j,1) == 1
                            k = k + 1;
                    
                            file.path = append(SubFolder(j).folder, filesep, SubFolder(j).name);
                    
                            stack = Core.LysosomeSegmentation(file);
                            stack.loadDataBioform(chan);
                            stack.showChannel;
                            stack.SegmentChannels(Threshold);
    
                            [Res, Trend] = stack.CalculatepH;
    
                            Results(k,:)  = Res;
                            Filename{k,1} = file.path;

                            Trendline(:,1) = Trend(:,1);
                            Trendline(:,end+1) = Trend(:,2);
                    
                            close all
                        end
                    end
    
    
                    ResultTable = table(string(Filename), Results(:,1), Results(:,2), Results(:,3), Results(:,4), Results(:,5), Results(:,6),  Results(:,7),  Results(:,8), Results(:,9), Results(:,10), Results(:,11), Results(:,12),...
                        Results(:,13), Results(:,14),  Results(:,15),  Results(:,16),  Results(:,17), Results(:,18), Results(:,19), Results(:,20), 'VariableNames', {'Filename', 'Signal total Ratio', 'SignalInLysotracker', 'SignalOutsideLysotracker', 'NPdensityInside', 'NPdensityOutside', 'NPintensityIn', 'NPintensityOut',...
                        'Lysotracker Fraction', 'DensityAcidic', 'DensityBasic', 'OverlapAcidic', 'OverlapBasic', 'AreaAcidic', 'AreaBasic', 'TotalParticleIntensityAcidic', 'TotalParticleIntensityBasic',...
                        'FractionOfParticlesAcidic', 'FractionOfParticlesBasic', 'FractionOfLysotrackerInAcidic', 'FractionOfLysotrackerInBasic'});
    
                    writetable(ResultTable, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'ResultspHSensor.xlsx'));
    
                    FinalResults.(CellineFolder) = [Results(:,1), Results(:,8), Results(:,9), Results(:,10)];

                    x = Trendline(:,1);
                    Y = Trendline(:,2:end);   % 200x20
                    y_mean = nanmean(Y, 2);        % 200x1
                    y_std  = nanstd(Y, 0, 2);      % 200x1
                    Fig2 = figure; hold on;
                    errorbar(x, y_mean, y_std, ...
                        'o-', ...
                        'LineWidth', 1.5, ...
                        'MarkerSize', 5, ...
                        'CapSize', 3);
                    xlabel('Relative intensity Lysotracker');
                    ylabel('Ratio mTFP/mCherry (Mean ± SD)');
                    title('Mean of 20 samples');
                    box on;
                    ylim([0 2.5])

                    saveas(Fig2, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'Trendline.png'));
                    saveas(Fig2, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'Trendline.svg'));
                    writematrix(Trendline, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'TrendLinepHSensor.xlsx'));

                    disp(append(CellineFolder, ' particle density acidic: ',num2str(mean(ResultTable.DensityAcidic))));
                    disp(append(CellineFolder, ' particle density basic: ',num2str(mean(ResultTable.DensityBasic))));
                    ResUlsts{r,1} = CellineFolder;
                    ResUlsts{r,2} = mean(ResultTable.OverlapAcidic);
                    ResUlsts{r,3} = mean(ResultTable.OverlapBasic);
    
                    writetable(ResultTable, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'ResultspHSensor.xlsx'));

                    Trendline = [];
                end
            end
        catch
        end
    end
end