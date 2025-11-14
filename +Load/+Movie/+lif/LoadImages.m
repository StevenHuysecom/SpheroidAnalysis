function LoadImages(file, chan)
    Folder = dir(file.path);
    %% Check if there was already a loading
    CheckForDir = [];
    for i = 3:size(Folder, 1)
        CheckForDir(end+1,1) = Folder(i).isdir;
    end

    if max(CheckForDir >= 1)
        disp('Data was already loaded and stored in folders')
    else
        %% make subfolder per sample       
        LifFiles = {};
        h = waitbar(0, 'initializing');
        for i = 3:size(Folder, 1)
            if and(strcmp(Folder(i).name(end-2:end), 'lif'), ~strcmp(Folder(i).name(1), '.'))
                LifFiles{end+1,1} = Folder(i).name;
                subfolder = Folder(i).name(1:end-4);
                mkdir([file.path,filesep,subfolder])
    
                FileName = append([file.path,filesep,subfolder,file.ext]);
                bfI = BioformatsImage(FileName);
    
                 for k = 1:bfI.seriesCount
                    bfI.series = k;
                    ImageName = append('Position', num2str(k));
                    mkdir([file.path,filesep,subfolder,filesep,ImageName])

                    for l = 1:bfI.sizeC
                        if strcmp(chan.(append('ch0', num2str(l))), 'ignore')
                            disp(append('Ignored channel ', num2str(l)))
                        elseif strcmp(chan.(append('ch0', num2str(l))), 'Membrane')
                            disp(append('Channel ', num2str(l), ' equals Membrane'))
                            for j = 1:bfI.sizeZ
                                waitbar(k./bfI.seriesCount, h, append('Extracting position ', num2str(k), ' of ', num2str(bfI.seriesCount),...
                                    ' / channel ', num2str(l), ' / slice ', num2str(j), ' of ', num2str(bfI.sizeZ)));
                                stack = double(getPlane(bfI, j, l, 1));
                                if all(size(stack) ~= [512, 512])
                                    if j == 1
                                        Scale = 512 ./ size(stack,1);
                                    end
                                    Membrane(:,:,j) = imresize(stack, Scale);
                                else
                                    Scale = 1;
                                    Membrane(:,:,j) = stack;
                                end
                            end
                            waitbar(k./bfI.seriesCount, h, append('Extracting position ', num2str(k), ' of ', num2str(bfI.seriesCount),...
                                    ' / channel ', num2str(l), ' / saving....'));
                            MatFileName = append(file.path,filesep,subfolder,filesep,ImageName,filesep,'Membrane.mat');
                            save(MatFileName, 'Membrane', '-v7.3');
                        elseif strcmp(chan.(append('ch0', num2str(l))), 'Particles')
                            disp(append('Channel ', num2str(l), ' equals Particles'))
                            for j = 1:bfI.sizeZ
                                waitbar(k./bfI.seriesCount, h, append('Extracting position ', num2str(k), ' of ', num2str(bfI.seriesCount),...
                                    ' / channel ', num2str(l), ' / slice ', num2str(j), ' of ', num2str(bfI.sizeZ)));
                                stack = double(getPlane(bfI, j, l, 1));
                                if all(size(stack) ~= [512, 512])
                                    if j == 1
                                        Scale = 512 ./ size(stack,1);
                                    end
                                    Particles(:,:,j) = imresize(stack, Scale);
                                else
                                    Scale = 1;
                                    Particles(:,:,j) = stack;
                                end
                            end
                            waitbar(k./bfI.seriesCount, h, append('Extracting position ', num2str(k), ' of ', num2str(bfI.seriesCount),...
                                    ' / channel ', num2str(l), ' / saving....'));
                            MatFileName = append(file.path,filesep,subfolder,filesep,ImageName,filesep,'Particles.mat');
                            save(MatFileName, 'Particles', '-v7.3');
                        else
                            disp(append('Channel ', num2str(l), ' equals ', chan.(append('ch0', num2str(l)))))
                            for j = 1:bfI.sizeZ
                                waitbar(k./bfI.seriesCount, h, append('Extracting position ', num2str(k), ' of ', num2str(bfI.seriesCount),...
                                    ' / channel ', num2str(l), ' / slice ', num2str(j), ' of ', num2str(bfI.sizeZ)));
                                stack = double(getPlane(bfI, j, l, 1));
                                if all(size(stack) ~= [512, 512])
                                    if j == 1
                                        Scale = 512 ./ size(stack,1);
                                    end
                                    Channel(:,:,j) = imresize(stack, Scale);
                                else
                                    Scale = 1;
                                    Channel(:,:,j) = stack;
                                end
                            end
                            waitbar(k./bfI.seriesCount, h, append('Extracting position ', num2str(k), ' of ', num2str(bfI.seriesCount),...
                                    ' / channel ', num2str(l), ' / saving....'));
                            MatFileName = append(file.path,filesep,subfolder,filesep,ImageName,filesep, chan.(append('ch0', num2str(l))), '.mat');
                            save(MatFileName, 'Channel', '-v7.3');
                        end
                    end

                    PxSizes = bfI.pxSize;
                    PxSizes(1:2) = PxSizes(1,2)./Scale;
                    MatFileName = append(file.path,filesep,subfolder,filesep,ImageName,filesep,'PxSizes.mat');
                    save(MatFileName, 'PxSizes')

                    clear Membrane Particles
                end
            end
        end
        close(h)
    end
end