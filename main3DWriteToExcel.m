clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
MainFolder = {'D:\Analysis_done\AuNPs@mSi@PEI'};
DimensionFolders = {'3D'};
HourFolders = {'3u', '6u', '24u', '48u'};
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
    
                if strcmp(ParticleFolder, 'A549')
                    Range1 = 'A3:Z503';
                    Range2 = 'B1:Z1';
                    Range3 = 'A1:A1';
                elseif strcmp(ParticleFolder, 'HeLa')
                    Range1 = 'AA3:AZ503';
                    Range2 = 'AB1:AZ1';
                    Range3 = 'AA1:AA1';
                elseif strcmp(ParticleFolder, 'KM12C')
                    Range1 = 'BA3:BZ503';
                    Range2 = 'BB1:BZ1';
                    Range3 = 'BA1:BA1';
                elseif strcmp(ParticleFolder, 'MCF7')
                    Range1 = 'CA3:CZ503';
                    Range2 = 'CB1:CZ1';
                    Range3 = 'CA1:CA1';
                end
    
                IntProfile = load(append(file.path, filesep, 'IntMatrix.mat'));
                IntProfile = IntProfile.IntMatrix;   
                writematrix(IntProfile, append(MainFolder{1}, filesep, 'Results3D.xlsx'), 'Sheet', HourFolder, 'Range', Range1);
                IntTot = load(append(file.path, filesep, 'SpheroidIntTotal.mat'));
                IntTot = IntTot.SpheroidInt;
                writematrix(IntTot, append(MainFolder{1}, filesep, 'Results3D.xlsx'), 'Sheet', HourFolder, 'Range', Range2);
                writecell({ParticleFolder}, append(MainFolder{1}, filesep, 'Results3D.xlsx'), 'Sheet', HourFolder, 'Range', Range3);
    
                
            catch
            end
        end
    end
end
