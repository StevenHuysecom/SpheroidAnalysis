clear;
close all;
clc;
%% User Input
file.ext  = '.lif';

MainFolder = {'D:\STEVEN - DATA RITA'};
DimensionFolders = {'life_dead'};
HourFolders = {'13May Dead Staining Spheroids', '24April Dead Staining Spheroids'};
ParticleFolders = {'PI_45mW_NewNPs', 'PI_45mW_No_NPs', 'PI_45mW_Old_NPs', 'PI_NoIrr_NewNPs',...
    'PI_NoIrr_No_NPs', 'PI_NoIrr_Old_NPs', 'PI_Irr_No_NPs', 'PI_Irr_NPs', 'PI_NoIrr_No_NPs',...
    'PI_NoIrr_NPs'};


%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.ch01 = 'ignore';
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
                            IntMatrix = [-100:200].';
                            for j = 1:size(SubFolder,1)
                                try
                                    if isSubDirColumn(j,1) == 1
                                        file.path = append(SubFolder(j).folder, filesep, SubFolder(j).name);
            
                                        PxSize = load(append(file.path, filesep, 'PxSizes.mat'));
                                        info.pxSizeXY = PxSize.PxSizes(1);
                                        info.pxSizeZ  = PxSize.PxSizes(3);
                                        stack = Core.Spheroid3D(file,info);
                                        stack.loadDataBioformSingleChan(chan);
                                        stack.showChannel;

                                        %% Get full spheroid intensity uptake
                                        [TotInt] = stack.GetFullIntSpheroids;
                                        ToAdd = [num2cell(TotInt), repmat({SubFolder(j).name}, size(TotInt, 1), 1)];
                                        SpheroidInt = [SpheroidInt; ToAdd];
                                        ToAdd = [];
                                        TotInt = [];
                                    end
                                catch 
                                end
                            end

                            ResultsTable = cell2table(SpheroidInt, 'VariableNames', {'Intensity','Volume (px)','Position'});
                            writetable(ResultsTable, append(MainFolder{1}, filesep, DimensionFolder,...
                                filesep, HourFolder, filesep,'IntensitySpheroids.xlsx'), 'Sheet', ParticleFolder);
                            close all
                        end
                    catch
                    end
                end
            catch
            end
        end
    end
end