clc
clear all
close all


MainFolder = {'D:\Data Uptake\AuNP@mSi@PEI'};
DimensionFolders = {'2D'};
HourFolders = {'3hour', '6hour', '24hour', '48hour'};
CellineFolders = {'A549', 'HeLa', 'KM12C', 'MCF7'};

figure()
for m = 1:numel(DimensionFolders)
    for a = 1:numel(HourFolders)
        CellIntAll = [];
        subplot(2,2,a)

        for r = 1:numel(CellineFolders)
            FolderPath = append(MainFolder, filesep, DimensionFolders{m}, filesep, HourFolders{a},...
                filesep, CellineFolders{r});
            Folder = dir(FolderPath{1,1});
            for i = 3:size(Folder)
                if Folder(i).isdir == true
                    FilePath = append(Folder(i).folder, filesep, Folder(i).name, filesep, 'CellInt.mat');
                    break
                else
                    continue
                end
            end
            CellInt = load(FilePath);
            CellInt = CellInt.CellInt;
            
            if size(CellInt, 1) < size(CellIntAll, 1)
                CellInt = [CellInt; nan(size(CellIntAll, 1) - size(CellInt, 1), 1)];
            elseif size(CellInt, 1) > size(CellIntAll, 1)
                CellIntAll = [CellIntAll; nan(size(CellInt, 1) - size(CellIntAll, 1), size(CellIntAll, 2))];
            end
            CellIntAll = [CellIntAll, CellInt];
        end
        boxplot(CellIntAll,'Labels',CellineFolders);
        ylim([0 1*10^6])
        title(HourFolders{a})
        ylabel('Intensity per cell (a.u.)')
        hold on
    end
end

figure()
for m = 1:numel(DimensionFolders)
    for a = 1:numel(HourFolders)
        MembrIntAll = [];
        subplot(2,2,a)

        for r = 1:numel(CellineFolders)
            FolderPath = append(MainFolder, filesep, DimensionFolders{m}, filesep, HourFolders{a},...
                filesep, CellineFolders{r});
            Folder = dir(FolderPath{1,1});
            for i = 3:size(Folder)
                if Folder(i).isdir == true
                    FilePath = append(Folder(i).folder, filesep, Folder(i).name, filesep, 'MembrInt.mat');
                    break
                else
                    continue
                end
            end
            MembrInt = load(FilePath);
            MembrInt = MembrInt.MembrInt;
            
            if size(MembrInt, 1) < size(MembrIntAll, 1)
                MembrInt = [MembrInt; nan(size(MembrIntAll, 1) - size(MembrInt, 1), 1)];
            elseif size(MembrInt, 1) > size(MembrIntAll, 1)
                MembrIntAll = [MembrIntAll; nan(size(MembrInt, 1) - size(MembrIntAll, 1), size(MembrIntAll, 2))];
            end
            MembrIntAll = [MembrIntAll, MembrInt];
        end
        boxplot(MembrIntAll,'Labels',CellineFolders);
        title(HourFolders{a})
        ylabel('Intensity between cells (a.u.)')
        hold on
    end
end
