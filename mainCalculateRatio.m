clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
file.runSegmentation = 'load'; %load or run
file.drawROI = 'off'; %'off', or channel name

MainFolder = {'D:\STEVEN - DATA RITA\ROS\13May_ROS Generation', 'D:\STEVEN - DATA RITA\ROS\11May_ROS Generation'};
SubFolders = {'NPs highC, 2min irradiation','NO NPs, 2min irradiation', 'No NPs, NO irradiation', ...
    'NPs highC, NO irradiation', 'NPs lowC, 2 min irradiation', 'NPs lowC, NO irradiation',...
    '11May_NewNPs_2min_15mW', '11May_NewNPs_5min_15mW', '11May_NewNPs_15min_15mW', '11May_NewNPs_NoLight', ...
    '11May_NoNPs_2min_15mW', '11May_NoNPs_5min_15mW', '11May_NoNPs_NoLight', '11May_OldNPs_2min_15mW', ...
    '11May_OldNPs_5min_15mW', '11May_OldNPs_15min_15mW', '11May_OldNPs_NoLight'};
SubsubFolders = {'0_min', '2_min', '5_min', '10_min'}; %, 'mSi', 'Control', 'CQ', 

ImageSize = [512, 512]; 
RatiobarLim = [0 1];
NPThreshold = 50;
Threshold = [0.01, 0.0001]; %[0-1], keep it under 0.15, intensity threshold for lysosomes (high = throw away dim/out-of-focus lysosomes)

%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.ch01 = 'SigChannel';
chan.ch02 = 'ignore';
chan.ch03 = 'RefChannel';
chan.ch04 = 'ignore';

%some parameters


%% Loading data

for aa = 1:numel(MainFolder)
    for a = 1:numel(SubFolders)
        SubFolder = SubFolders{a};
        for r = 1:numel(SubsubFolders)
            try
                SubsubFolder = SubsubFolders{r};
                Path = append(MainFolder{aa}, filesep, SubFolder, filesep, SubsubFolder);
                file.path = Path;
        
                Load.Movie.lif.LoadImages(file, chan);
                Load.Movie.DrawApplyROI(file, chan, 1);
            catch
    
            end
        end
    end
end


for aa = 1:numel(MainFolder)
    for a = 1:numel(SubFolders)
        SubFolder = SubFolders{a};
        if exist(append(MainFolder{aa}, filesep, SubFolder))
            for r = 1:numel(SubsubFolders)
                try
                    SubsubFolder = SubsubFolders{r};
                    Path = append(MainFolder{aa}, filesep, SubFolder, filesep, SubsubFolder);
                    file.path = Path;
                    CurrentFolder = dir(file.path);
                    CurrentFolder(1:2) = [];
                    isDirColumn = [CurrentFolder.isdir]';
                    for i = 1:size(CurrentFolder,1)
                        if isDirColumn(i,1) == 1
                            SubcurrentFolder = dir(append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name));
                            SubcurrentFolder(1:2) = [];
                            isSubDirColumn = [SubcurrentFolder.isdir]';
            
                            Results  = nan([ImageSize(1)^2, size(SubcurrentFolder,1)]);
                            Filename = {};
                            k = 1;   % index for valid entries
                            
                            for j = 1:size(SubcurrentFolder,1)
                            
                                if isSubDirColumn(j,1) == 1
                                    try
                                        file.path = append(SubcurrentFolder(j).folder, filesep, SubcurrentFolder(j).name);
                                
                                        stack = Core.LysosomeSegmentation(file);
                                        stack.loadDataBioform(chan);
                                        stack.showChannel;
                                        stack.Segmentation2;
                
                                        [Ratios] = stack.CalculateRatio(RatiobarLim, NPThreshold);
                                        Results(1:size(Ratios, 1), j) = Ratios;
                                        Ratios = [];
         
                                        close all
                                    catch
                                        disp('fail')
                                    end
                                end
                            end
        
                            Results(all(isnan(Results),2),:) = [];
                            AvResults = median(Results, 1, 'omitnan');
            
                            Range = append('A3:ZZ', num2str(size(Results, 1)+3));
                            % writematrix(Results, append(MainFolder{1}, filesep, SubFolder, filesep, 'ResultsRatio.xlsx'), 'Sheet', SubsubFolder, 'Range', Range);
                            writematrix(AvResults, append(MainFolder{aa}, filesep, SubFolder, filesep, 'ResultsRatio.xlsx'), 'Sheet', SubsubFolder, 'Range', 'A1:ZZ1');
                        end
                    end
        
                    [Hist,~] = histcounts(Results(:), [0:0.005:1]);
                    HistNorm = Hist./max(Hist);
                    HistBins = [0.005:0.005:1];
                    HistCell{r} = HistNorm;
                    HistLabel{r} = strrep(SubsubFolder, '_', ' ');
            
                    Histmatrix(:,r) = HistNorm';
                    Hist = [];
                    HistNorm = [];
        
                catch
                end
            end
            Fig = figure();
            fails = zeros([1 r]);
            for r = 1:numel(SubsubFolders)
                try
                    plot(HistBins, HistCell{r});        
                catch
                end
                hold on
            end
            try
                HistLabel  = HistLabel(~cellfun('isempty',HistLabel));
                legend(HistLabel)
            catch
            end
            ylabel('Occurance')
            xlabel('Ratio')
            title('Histogram of ratios')
            saveas(Fig, append(MainFolder{aa}, filesep, SubFolder, filesep, 'HistogramRatios.png'))
            saveas(Fig, append(MainFolder{aa}, filesep, SubFolder, filesep, 'HistogramRatios.svg'))
        
            HistCell = [];
            HistLabel = [];
            writematrix(Histmatrix, append(MainFolder{aa}, filesep, SubFolder, filesep, 'ResultsRatioHistogram.xlsx'));
        end
    end
end