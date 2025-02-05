clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
MainFolder = {'D:\Steven'};
DimensionFolders = {'2D'};
HourFolders = {'3hour', '24hour', '48hour'}; %
ParticleFolders = {'A549', 'HeLa', 'KM12C', 'MCF7'};

BigMatrix = [];

for m = 1:numel(DimensionFolders)
    DimensionFolder = DimensionFolders{m};
    for a = 1:numel(HourFolders)
        HourFolder = HourFolders{a};
        SliceMatrix = readmatrix(append(MainFolder{1,1}, filesep, DimensionFolder, filesep, 'Slices.xlsx'), 'Sheet', HourFolder);
        for r = 1:numel(ParticleFolders)
            try
                ParticleFolder = ParticleFolders{r};
                Path = append(MainFolder, filesep, DimensionFolder, filesep, HourFolder,...
                    filesep, ParticleFolder);
                file.path = Path{1,1};

                Folder = dir(file.path);

                NumCell = [];
                CellVolumes = [];
                CellInt = [];
                CellDens = [];
                DeleteValues = [];
                for i = 3:size(Folder,1)
                    if Folder(i).isdir == 1
                        NewFolder = append(Folder(i).folder, filesep, Folder(i).name);
                        NewFolder = dir(NewFolder);
                        for z = 3:size(NewFolder, 1)
                            if NewFolder(z).isdir == 1
                                FileName = append(NewFolder(z).folder, filesep, NewFolder(z).name, filesep, 'MembraneSegmentation.mat');
                                MembrSegm = load(FileName);
                                MembrSegm = MembrSegm.ws;

                                Position = NewFolder(z).name;
                                x = str2num(Position(9:end));

                                Slice = SliceMatrix(x, r);

                                if isnan(Slice)
                                    DeleteValues(end+1,1) = 1;
                                else
                                    DeleteValues(end+1,1) = 0;
                                end
                                

                                stats = regionprops(MembrSegm, 'Area');
                                VolList = (struct2array(stats)).';
                                VolList(VolList < 2700) = [];
                                NumCell = [NumCell; size(VolList, 1)];
                                CellVolumes = [CellVolumes; VolList];
                            end
                        end 
                        
                        NumCell(DeleteValues == 1) = [];
                        Check = sum(NumCell);
                        CellVolumes = load(append(NewFolder(1).folder, filesep, "VolumeList.mat"));
                        CellVolumes = CellVolumes.VolumeList;

                        if strcmp(ParticleFolder, 'A549')
                            Range1 = 'A1:A1';
                            Range2 = 'B2:B50';
                            Range3 = 'C2:C2500';
                            Range4 = 'D2:D2500';
                            Range5 = 'E2:E2500';
                            Range6 = 'B1:B1';
                            Range7 = 'D1:D1';
                            Range8 = 'E1:E1';
                            Range9 = 'C1:C1';
                        elseif strcmp(ParticleFolder, 'HeLa')
                            Range1 = 'AA1:AA1';
                            Range2 = 'AB2:AB50';
                            Range3 = 'AC2:AC2500';
                            Range4 = 'AD2:AD2500';
                            Range5 = 'AE2:AE2500';
                            Range6 = 'AB1:AB1';
                            Range7 = 'AD1:AD1';
                            Range8 = 'AE1:AE1';
                            Range9 = 'AC1:AC1';
                        elseif strcmp(ParticleFolder, 'KM12C')
                            Range1 = 'BA1:BA1';
                            Range2 = 'BB2:BB50';
                            Range3 = 'BC2:BC2500';
                            Range4 = 'BD2:BD2500';
                            Range5 = 'BE2:BE2500';
                            Range6 = 'BB1:BB1';
                            Range7 = 'BD1:BD1';
                            Range8 = 'BE1:BE1';
                            Range9 = 'BC1:BC1';
                        elseif strcmp(ParticleFolder, 'MCF7')
                            Range1 = 'CA1:CA1';
                            Range2 = 'CB2:CB50';
                            Range3 = 'CC2:CC2500';
                            Range4 = 'CD2:CD2500';
                            Range5 = 'CE2:CE2500';
                            Range6 = 'CB1:CB1';
                            Range7 = 'CD1:CD1';
                            Range8 = 'CE1:CE1';
                            Range9 = 'CC1:CC1';
                        end
            
                        %% abs intensities
                        writecell({ParticleFolder}, append(MainFolder{1}, filesep, 'Results2D_absInt.xlsx'), 'Sheet', HourFolder, 'Range', Range1);
                        writematrix(NumCell, append(MainFolder{1}, filesep, 'Results2D_absInt.xlsx'), 'Sheet', HourFolder, 'Range', Range2);
                        writematrix(CellVolumes, append(MainFolder{1}, filesep, 'Results2D_absInt.xlsx'), 'Sheet', HourFolder, 'Range', Range3);
        
                        Path2Load = append(NewFolder(1).folder, filesep, 'CellInt.mat');
                        CellInt = load(Path2Load);
                        CellInt = CellInt.CellInt;
                        CellInt = CellInt(1:end, 1);
                        writematrix(CellInt, append(MainFolder{1}, filesep, 'Results2D_absInt.xlsx'), 'Sheet', HourFolder, 'Range', Range4);
        
                        Path2Load2 = append(NewFolder(1).folder, filesep, 'MembrInt.mat');
                        MembrInt = load(Path2Load2);
                        MembrInt = MembrInt.MembrInt;
                        writematrix(MembrInt, append(MainFolder{1}, filesep, 'Results2D_absInt.xlsx'), 'Sheet', HourFolder, 'Range', Range5);
        
                        writecell({'Number of cells'}, append(MainFolder{1}, filesep, 'Results2D_absInt.xlsx'), 'Sheet', HourFolder, 'Range', Range6);
                        writecell({'Cell Intensities'}, append(MainFolder{1}, filesep, 'Results2D_absInt.xlsx'), 'Sheet', HourFolder, 'Range', Range7);
                        writecell({'Membrane Intensities'}, append(MainFolder{1}, filesep, 'Results2D_absInt.xlsx'), 'Sheet', HourFolder, 'Range', Range8);
                        writecell({'Cell Volumes'}, append(MainFolder{1}, filesep, 'Results2D_absInt.xlsx'), 'Sheet', HourFolder, 'Range', Range9);

                        %% intensity density
                        writecell({ParticleFolder}, append(MainFolder{1}, filesep, 'Results2D_density.xlsx'), 'Sheet', HourFolder, 'Range', Range1);
                        writematrix(NumCell, append(MainFolder{1}, filesep, 'Results2D_density.xlsx'), 'Sheet', HourFolder, 'Range', Range2);
                        writematrix(CellVolumes, append(MainFolder{1}, filesep, 'Results2D_density.xlsx'), 'Sheet', HourFolder, 'Range', Range3);
        
                        Path2Load = append(NewFolder(1).folder, filesep, 'CellDens.mat');
                        CellDens = load(Path2Load);
                        CellDens = CellDens.CellDens;
                        CellDens = CellDens(1:size(CellVolumes,1), 1);
                        writematrix(CellDens, append(MainFolder{1}, filesep, 'Results2D_density.xlsx'), 'Sheet', HourFolder, 'Range', Range4);
        
                        Path2Load2 = append(NewFolder(1).folder, filesep, 'MembrDens.mat');
                        MembrDens = load(Path2Load2);
                        MembrDens = MembrDens.MembrDens;
                        writematrix(MembrDens, append(MainFolder{1}, filesep, 'Results2D_density.xlsx'), 'Sheet', HourFolder, 'Range', Range5);
        
                        writecell({'Number of cells'}, append(MainFolder{1}, filesep, 'Results2D_density.xlsx'), 'Sheet', HourFolder, 'Range', Range6);
                        writecell({'Cell Intensities'}, append(MainFolder{1}, filesep, 'Results2D_density.xlsx'), 'Sheet', HourFolder, 'Range', Range7);
                        writecell({'Membrane Intensities'}, append(MainFolder{1}, filesep, 'Results2D_density.xlsx'), 'Sheet', HourFolder, 'Range', Range8);
                        writecell({'Cell Volumes'}, append(MainFolder{1}, filesep, 'Results2D_density.xlsx'), 'Sheet', HourFolder, 'Range', Range9);

                    else
                    end 
                end     
            catch
            end
        end
    end
end