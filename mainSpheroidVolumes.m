clear;
close all;
clc;
%% User Input
file.ext  = '.lif';
MainFolder = {'F:\Data Uptake\AuNP@mSi@PEI'};
DimensionFolders = {'3D'};
HourFolders = {'3hour', '6hour', '24hour', '48hour'};
ParticleFolders = {'A549', 'Hela', 'KM12C', 'MCF7'};

for m = 1:numel(DimensionFolders)
    for a = 1:numel(HourFolders)
        for r = 1:numel(ParticleFolders)
            Path = append(MainFolder, filesep, DimensionFolders{m}, filesep, HourFolders{a}, filesep, ParticleFolders{r});

            Folder = dir(Path{1,1});

            for z = 3:size(Folder, 1)
                if Folder(z).isdir == 1
                        NewFolder =  dir(append(Folder(z).folder, filesep, Folder(z).name));
                        Volume = [];
                        msg = append(HourFolders{a}, ' ', ParticleFolders{r});
                        h = waitbar(0, msg);
                        for q = 3:size(NewFolder, 1)
                            waitbar(q./size(NewFolder,1), h, msg);
                            try
                                if NewFolder(q).isdir == 1
                                    f = waitbar(0, 'loading membrane');
                                    FilePath = append(NewFolder(q).folder, filesep, NewFolder(q).name);    
                                    membrane = load(append(FilePath, filesep, 'Membrane.mat'));
                                    membrane = membrane.Membrane;
                                    waitbar(0.2, f, 'Thresholding');
                                    membrane(membrane < 10) = 0;
                                    waitbar(0.3, f, 'dilation');
                                    se = strel('cube', 2);
                                    membrane = imdilate(membrane, se);
                                    waitbar(0.5, f, 'medfilt');
                                    membrane = medfilt3(membrane, [5 5 5]);
                                    waitbar(0.6 ,f, 'closing gaps');
                                    membrane = bwareaopen(membrane, 500000);
                                    for i = 1:size(membrane, 3)
                                        membrane(:,:,i) = imfill(membrane(:,:,i), "holes");
                                        stats = regionprops(membrane(:,:,i), 'Area');
                                        if max(struct2array(stats)) > 5000
                                            membrane(:,:,i) = bwareaopen(membrane(:,:,i), 5000);
                                        end
                                    end
                                    waitbar(0.8 ,f, 'filling holes');
                                    membrane = imfill(membrane, "holes");    
                                    waitbar(0.9 ,f, 'calculating volume');
                                    membraneList = membrane(:);
                                    PixelNumber = sum(membraneList);
                                    
                                    PxSizes = load(append(FilePath, filesep, 'PxSizes.mat'));
                                    PxSizes = PxSizes.PxSizes;
                                    PxVolume = PxSizes(1)*PxSizes(2)*PxSizes(3); 
                                    waitbar(0.9 ,f, 'calculating volume');
                                    Volume(1, end+1) = PixelNumber*PxVolume;
                                    close(f)
                                end
                            catch
                            end
                        end
                        close(h)
    
                        if a == 1
                            Sheet = '3hour';
                        elseif a == 2
                            Sheet = '6hour';
                        elseif a == 3
                            Sheet = '24hour';
                        elseif a == 4
                            Sheet = '48hour';
                        end
    
                        if r == 1
                            Range = 'B1:Z1';
                        elseif r == 2
                            Range = 'AB1:AZ1';
                        elseif r == 3
                            Range = 'BB1:BZ1';
                        elseif r == 4
                            Range = 'CB1:CZ1';
                        end
    
                        Filename = append(MainFolder, filesep, 'ResultsVolumes.xlsx');
                        Filename = Filename{1,1};
                        writematrix(Volume, Filename,'Sheet', Sheet, 'Range', Range);
                end
            end
        end
    end
end







