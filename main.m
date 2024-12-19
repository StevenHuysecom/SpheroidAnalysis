clear;
close all;
clc;
%% User Input
file.path = 'C:\Users\Windows 11\OneDrive - KU Leuven\Documents\KU Leuven\PhD\data\Testdata_Maria\KM12C_1';
file.ext  = '.tif';

info.pxSizeXY = 454.5; 
info.pxSizeZ  = 1000;

%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.ch00 = 'Membrane';
chan.ch02 = 'Particles';
chan.ch03 = 'ignore';
chan.ch04 = 'ignore';

%% Loading data
%rendering3D.compile3DRendering();

stack = Core.OrganoidSegmentation(file,info);

stack.loadData(chan);

stack.showChannel;

%% Segmentation of nuclei

%segmented = stack.segmentNuclei;

%% Watershed attempt 1
gAdapt = stack.segmentNucleiWatershed;



%% Intensity

[Intensity, IntensityImage] = stack.getIntensityAround();



%% make video
contour = stack.results.cellContour;

figure(1)
colormap(jet) 
vidfile = VideoWriter('testmovie.mp4','MPEG-4');
open(vidfile);
for i = 1:size(IntensityImage,3)
    
    imagesc(IntensityImage(:,:,i)),
    hold on
    currContour = contour{i};
    for j = 1:length(currContour)
        plot(currContour{j}(:,2),currContour{j}(:,1),'w')
    end
    caxis([0 max(IntensityImage(:))])
    colorbar
    axis image
    drawnow
    F(i) = getframe(gcf); 
    writeVideo(vidfile,F(i));
end
close(vidfile)

%% Plotting
%segmentation
fr = 49;
figure 
subplot(1,3,1)
imagesc(stack.channels.cell(:,:,fr))
axis image
subplot(1,3,2)
imagesc(stack.channels.nucleus(:,:,fr))
axis image
subplot(1,3,3)
imagesc(stack.results.cellMask(:,:,fr))
axis image

%overlay
figure 
subplot(1,2,1)
imagesc(stack.channels.cell(:,:,fr))
contour = stack.results.cellContour{fr};
hold on
for i = 1: length(contour)
    plot(contour{i}(:,2),contour{i}(:,1),'w','Linewidth',2)
end
axis image
subplot(1,2,2)
imagesc(stack.channels.nucleus(:,:,fr))
hold on
for i = 1: length(contour)
    plot(contour{i}(:,2),contour{i}(:,1),'w','Linewidth',2)
end
axis image

% 3D visualization
rendering3D.compile3DRendering();

stack.renderCell3D(1);


%% plot segmentation
test = bwlabeln(stack.results.cellMask);
figure
imagesc(test(:,:,30))
colormap('colorcube')


%% Intensity calculation






