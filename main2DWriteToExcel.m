clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
MainFolder = {'D:\Data Uptake\AuNP@mSi@PEI'};
DimensionFolders = {'2D'};
HourFolders = {'3hour', '6hour', '24hour', '48hour'};
ParticleFolders = {'A549', 'HeLa', 'KM12C', 'MCF7'};

BigMatrix = [];

for m = 1:numel(DimensionFolders)
    DimensionFolder = DimensionFolders{m};
    for a = 1:numel(HourFolders)
        HourFolder = HourFolders{a};
        for r = 1:numel(ParticleFolders)
            try
                ParticleFolder = ParticleFolders{r};
                Path = append(MainFolder, filesep, DimensionFolder, filesep, HourFolder,...
                    filesep, ParticleFolder);
                file.path = Path{1,1};

                Folder = dir(file.path);

                NumCell = [];
                CellVolumes = [];
                for i = 3:size(Folder,1)
                    if Folder(i).isdir == 1
                        FileName = append(Folder(i).folder, filesep, Folder(i).name, filesep, 'MembraneSegmentation.mat');
                        MembrSegm = load(FileName);
                        MembrSegm = MembrSegm.MembraneSegmentation;

                        stats = regionprops3(MembrSegm, 'Volume');
                        NumCell = [NumCell; size(stats, 1)];
                        CellVolumes = [CellVolumes; stats];
                    else
                    end 
                end
    
                if strcmp(ParticleFolder, 'A549')
                    Range1 = 'A1:A1';
                    Range2 = 'B2:B50';
                    Range3 = 'C1:C2500';
                    Range4 = 'D2:D2500';
                    Range5 = 'E2:E2500';
                    Range6 = 'B1:B1';
                    Range7 = 'D1:D1';
                    Range8 = 'E1:E1';
                elseif strcmp(ParticleFolder, 'HeLa')
                    Range1 = 'AA1:AA1';
                    Range2 = 'AB2:AB50';
                    Range3 = 'AC1:AC2500';
                    Range4 = 'AD2:AD2500';
                    Range5 = 'AE2:AE2500';
                    Range6 = 'AB1:AB1';
                    Range7 = 'AD1:AD1';
                    Range8 = 'AE1:AE1';
                elseif strcmp(ParticleFolder, 'KM12C')
                    Range1 = 'BA1:BA1';
                    Range2 = 'BB2:BB50';
                    Range3 = 'BC1:BC2500';
                    Range4 = 'BD2:BD2500';
                    Range5 = 'BE2:BE2500';
                    Range6 = 'BB1:BB1';
                    Range7 = 'BD1:BD1';
                    Range8 = 'BE1:BE1';
                elseif strcmp(ParticleFolder, 'MCF7')
                    Range1 = 'CA1:CA1';
                    Range2 = 'CB2:CB50';
                    Range3 = 'CC1:CC2500';
                    Range4 = 'CD2:CD2500';
                    Range5 = 'CE2:CE2500';
                    Range6 = 'CB1:CB1';
                    Range7 = 'CD1:CD1';
                    Range8 = 'CE1:CE1';
                end
    

                writecell({ParticleFolder}, append(MainFolder{1}, filesep, 'Results2D.xlsx'), 'Sheet', HourFolder, 'Range', Range1);
                writematrix(NumCell, append(MainFolder{1}, filesep, 'Results2D.xlsx'), 'Sheet', HourFolder, 'Range', Range2);
                writetable(CellVolumes, append(MainFolder{1}, filesep, 'Results2D.xlsx'), 'Sheet', HourFolder, 'Range', Range3);

                Path2Load = append(file.path, filesep, 'CellInt.mat');
                CellInt = load(Path2Load);
                CellInt = CellInt.CellInt;
                CellInt = CellInt(1:size(CellVolumes,1), 1);
                writematrix(CellInt, append(MainFolder{1}, filesep, 'Results2D.xlsx'), 'Sheet', HourFolder, 'Range', Range4);

                Path2Load2 = append(file.path, filesep, 'MembrInt.mat');
                MembrInt = load(Path2Load2);
                MembrInt = MembrInt.MembrInt;
                writematrix(MembrInt, append(MainFolder{1}, filesep, 'Results2D.xlsx'), 'Sheet', HourFolder, 'Range', Range5);

                writecell({'Number of cells'}, append(MainFolder{1}, filesep, 'Results2D.xlsx'), 'Sheet', HourFolder, 'Range', Range6);
                writecell({'Cell Intensities'}, append(MainFolder{1}, filesep, 'Results2D.xlsx'), 'Sheet', HourFolder, 'Range', Range7);
                writecell({'Membrane Intensities'}, append(MainFolder{1}, filesep, 'Results2D.xlsx'), 'Sheet', HourFolder, 'Range', Range8);

    
                
            catch
            end
        end
    end
end