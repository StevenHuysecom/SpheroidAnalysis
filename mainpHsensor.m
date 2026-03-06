clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
file.runSegmentation = 'load'; %load or run
file.drawROI = 'off'; %'off', or channel name

MainFolder = {'D:\mini\FIREphly sensor'};
HourFolders = {'20260218'};
CellineFolders = {'mSiPEI'}; %, 'mSi', 'Control', 'CQ', 

%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.ch01 = 'mTFP1';
chan.ch02 = 'mCherry';
chan.ch03 = 'Particles';
chan.ch04 = 'Lysotracker';

%some parameters
slice = 1; %which slice of the 3D stack to select the ROI on
Threshold = [0.01, 0.05, 0.10, 0.20]; %[0-1], keep it under 0.15, intensity threshold for lysosomes (high = throw away dim/out-of-focus lysosomes)

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
                    k = 1;   % index for valid entries
                    
                    for j = 1:size(SubFolder,1)
                    
                        if isSubDirColumn(j,1) == 1
                            try
                                
                        
                                file.path = append(SubFolder(j).folder, filesep, SubFolder(j).name);
                        
                                stack = Core.LysosomeSegmentation(file);
                                stack.loadDataBioform(chan);
                                stack.showChannel;
                                stack.SegmentChannels(Threshold);
        
                                [Res, Trend, TrendParticle, PxCounts] = stack.CalculatepH;
        
                                Results(k,:)  = Res;
                                Filename{k,1} = file.path;
    
                                Trendline(:,1) = Trend(:,1);
                                Trendline(:,end+1) = medfilt1(Trend(:,2), 35);
    
                                TrendParticleline(:,1) = TrendParticle(:,1);
                                TrendParticleline(:,end+1) = medfilt1(TrendParticle(:,2), 35);

                                PxCountsAll(:,1) = TrendParticle(:,1);
                                PxCountsAll(:,end+1) = PxCounts;
                    
                                close all
                                k = k + 1;
                            catch
                                disp('fail')
                            end
                        end
                    end
    
    
                    ResultTable = table(string(Filename), Results(:,1), Results(:,2), Results(:,3), Results(:,4), Results(:,5), Results(:,6),  Results(:,7),  Results(:,8), Results(:,9), Results(:,10), Results(:,11), Results(:,12),...
                        Results(:,13), Results(:,14),  Results(:,15),  Results(:,16),  Results(:,17), Results(:,18), Results(:,19), Results(:,20), 'VariableNames', {'Filename', 'Signal total Ratio', 'SignalInLysotracker', 'SignalOutsideLysotracker', 'NPdensityInside', 'NPdensityOutside', 'NPintensityIn', 'NPintensityOut',...
                        'Lysotracker Fraction', 'DensityAcidic', 'DensityBasic', 'OverlapAcidic', 'OverlapBasic', 'AreaAcidic', 'AreaBasic', 'TotalParticleIntensityAcidic', 'TotalParticleIntensityBasic',...
                        'FractionOfParticlesAcidic', 'FractionOfParticlesBasic', 'FractionOfLysotrackerInAcidic', 'FractionOfLysotrackerInBasic'});
    
                    writetable(ResultTable, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'ResultspHSensor.xlsx'));
    
                    FinalResults.(CellineFolder) = [Results(:,1), Results(:,8), Results(:,9), Results(:,10)];

                    x = Trendline(:,1);
                    
                    Trendline(:, nanmax(Trendline, [], 1) > 0.5) = [];
                    for hh = 1:size(Trendline, 2)
                        [~, PeakLoc(hh)] = max(Trendline(:,hh));
                    end
                    Trendline(:, PeakLoc' > 60) = [];
                  
                    PxCountsAv = mean(PxCountsAll(:, 2:end), 2, 'omitnan');
                    PxCountsNorm = PxCountsAv./max(PxCountsAv); 
                    [~, MaxPixLoc] = max(PxCountsAv);
                    Trendline(1:MaxPixLoc, :) = NaN;

                    Y = Trendline(:,2:end);
                    y_mean = smooth(nanmean(Y, 2));
                    y_std  = smooth(nanstd(Y, 0, 2));

                    alphaVals = log(PxCountsAv)./max(log(PxCountsAv));
                    alphaVals(alphaVals < 0) = 0;
                    alphaVals(isinf(abs(alphaVals))) = 0;
                    alphaVals(isnan(alphaVals)) = 0;
                    
                    Fig3 = figure;
                    subplot(5,1,[1 4])
                    set(gcf,'Renderer','opengl');
                    hold on;

                    dx = nanmean(diff(x));
                    halfWidth = dx/2;
                    
                    %% ---- Plot SD as thick vertical patches (light blue) ----
                    for ii = 1:length(x)
                        x_patch = [x(ii)-halfWidth, x(ii)+halfWidth, ...
                                   x(ii)+halfWidth, x(ii)-halfWidth];
                        
                        y_patch = [y_mean(ii)-y_std(ii), ...
                                   y_mean(ii)-y_std(ii), ...
                                   y_mean(ii)+y_std(ii), ...
                                   y_mean(ii)+y_std(ii)];
                        
                        p = patch(x_patch, y_patch, [0.6 0.85 1], ...
                            'EdgeColor','none');
                        
                        p.FaceAlpha = alphaVals(ii);
                    end
                    
                    %% ---- Plot mean trend line (dark blue, variable alpha) ----
                    for ii = 1:length(x)-1
                        c = [0 0 0.5 alphaVals(ii)];  % dark blue with alpha
                        
                        line(x(ii:ii+1), y_mean(ii:ii+1), ...
                            'Color', c, ...
                            'LineWidth', 2);
                    end

                    ylabel('Relative intensity lysotracker');
                    title(['Mean of ', num2str(size(Trendline,2)-1), ' samples']);
                    ylim([0 0.35])
                    xlim([0 2.5])
                    box on;

                    subplot(5,1,5)
                    plot(x, PxCountsAv)
                    xlabel('Ratio mTFP/mCherry (Mean ± SD)');
                    xlim([0 3])
                    ylabel('Pixel Count');

                    saveas(Fig3, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'Trendline.png'));
                    saveas(Fig3, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'Trendline.svg'));
                    writematrix(Trendline, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'TrendLinepHSensor.xlsx'));








                    x = TrendParticleline(:,1);
                    
                    TrendParticleline(:, nanmax(TrendParticleline, [], 1) > 0.5) = [];

                    PxCountsAv = mean(PxCountsAll(:, 2:end), 2, 'omitnan');
                    PxCountsNorm = PxCountsAv./max(PxCountsAv); 
                    [~, MaxPixLoc] = max(PxCountsAv);
                    Trendline(1:MaxPixLoc, :) = NaN;

                    Y = TrendParticleline(:,2:end);
                    y_mean = smooth(nanmean(Y, 2));
                    y_std  = smooth(nanstd(Y, 0, 2));

                    alphaVals = log(PxCountsAv)./max(log(PxCountsAv));
                    alphaVals(alphaVals < 0) = 0;
                    alphaVals(isinf(abs(alphaVals))) = 0;
                    alphaVals(isnan(alphaVals)) = 0;
                    
                    Fig3 = figure;
                    subplot(5,1,[1 4])
                    set(gcf,'Renderer','opengl');
                    hold on;

                    dx = nanmean(diff(x));
                    halfWidth = dx/2;
                    
                    %% ---- Plot SD as thick vertical patches (light blue) ----
                    for ii = 1:length(x)
                        x_patch = [x(ii)-halfWidth, x(ii)+halfWidth, ...
                                   x(ii)+halfWidth, x(ii)-halfWidth];
                        
                        y_patch = [y_mean(ii)-y_std(ii), ...
                                   y_mean(ii)-y_std(ii), ...
                                   y_mean(ii)+y_std(ii), ...
                                   y_mean(ii)+y_std(ii)];
                        
                        p = patch(x_patch, y_patch, [0.6 0.85 1], ...
                            'EdgeColor','none');
                        
                        p.FaceAlpha = alphaVals(ii);
                    end
                    
                    %% ---- Plot mean trend line (dark blue, variable alpha) ----
                    for ii = 1:length(x)-1
                        c = [0 0 0.5 alphaVals(ii)];  % dark blue with alpha
                        
                        line(x(ii:ii+1), y_mean(ii:ii+1), ...
                            'Color', c, ...
                            'LineWidth', 2);
                    end

                    ylabel('Relative intensity particles');
                    title(['Mean of ', num2str(size(Trendline,2)-1), ' samples']);
                    ylim([0 0.45])
                    xlim([0 2.5])
                    box on;

                    subplot(5,1,5)
                    plot(x, PxCountsAv)
                    xlabel('Ratio mTFP/mCherry (Mean ± SD)');
                    xlim([0 2.5])
                    ylabel('Pixel Count');
                    
                    saveas(Fig3, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'TrendlineParticleNew.png'));
                    saveas(Fig3, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'TrendlineParticleNew.svg'));
                    writematrix(TrendParticleline, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'TrendLinepHSensorParticleNew.xlsx'));

                    TrendParticleline(:, nanmax(TrendParticleline, [], 1) > 10) = [];

                    disp(append(CellineFolder, ' particle density acidic: ',num2str(mean(ResultTable.DensityAcidic))));
                    disp(append(CellineFolder, ' particle density basic: ',num2str(mean(ResultTable.DensityBasic))));
                    ResUlsts{r,1} = CellineFolder;
                    ResUlsts{r,2} = mean(ResultTable.OverlapAcidic);
                    ResUlsts{r,3} = mean(ResultTable.OverlapBasic);
    
                    writetable(ResultTable, append(CurrentFolder(i).folder, filesep, CurrentFolder(i).name, filesep, 'ResultspHSensor.xlsx'));

                    Trendline = [];
                    TrendParticleLine = [];
                end
            end
        catch
        end
    end
end