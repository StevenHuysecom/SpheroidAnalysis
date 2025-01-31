OutputFolder = 'E:\Steven\2D';
CellLine = 'A549';
fr = 20;

close all

MembraneImg = figure();
imagesc(Membrane(:,:,fr))
colormap('hot')
axis off
Filename = append(OutputFolder, filesep, CellLine, '_2D_24hour_', 'MembraneImage.png');
saveas(MembraneImg, Filename);

ParticleImage = figure();
imagesc(Particles(:,:,fr))
colormap('hot')
axis off
Filename = append(OutputFolder, filesep, CellLine, '_2D_24hour_', 'ParticleImage.png');
saveas(ParticleImage, Filename);

MembraneSegmentImage = figure();
imagesc(MembraneSegment(:,:,fr))
axis off
Filename = append(OutputFolder, filesep, CellLine, '_2D_24hour_', 'MembraneSegment.png');
saveas(MembraneSegmentImage, Filename);

CellSegmentImage = figure();
imagesc(ws(:,:,fr))
axis off
Filename = append(OutputFolder, filesep, CellLine, '_2D_24hour_', 'CellSegment.png');
saveas(CellSegmentImage, Filename);

gBW = ws;
contour = cell(1,size(gBW,3));
for i = 1:size(gBW,3)
    currBW = gBW(:,:,i);
    
    %Get the largest area
    cBWarea = regionprops(currBW,'Area');
    [~,idx2BiggestArea] = max(cell2mat({cBWarea.Area}));
    
    if isempty(idx2BiggestArea)
    else
        
        [pContour] = bwboundaries(currBW);
        for j = 1:length(pContour)
            contour{i}{j} = pContour{j};
        end
    end
end
cContour = contour{1,fr};
ParticleContour = figure();
imagesc(Particles(:,:,fr))
colormap('hot')
axis image
for i = 1:length(cContour)
  hold on
  plot(cContour{i}(:,2),cContour{i}(:,1),'w','Linewidth',1)     
end
axis off
Filename = append(OutputFolder, filesep, CellLine, '_2D_24hour_', 'ParticleContour.png');
exportgraphics(ParticleContour, Filename);

gBW = ws;
contour = cell(1,size(gBW,3));
for i = 1:size(gBW,3)
    currBW = gBW(:,:,i);
    
    %Get the largest area
    cBWarea = regionprops(currBW,'Area');
    [~,idx2BiggestArea] = max(cell2mat({cBWarea.Area}));
    
    if isempty(idx2BiggestArea)
    else
        
        [pContour] = bwboundaries(currBW);
        for j = 1:length(pContour)
            contour{i}{j} = pContour{j};
        end
    end
end
cContour = contour{1,fr};
CellContour = figure();
imagesc(Membrane(:,:,fr))
colormap('hot')
axis image
for i = 1:length(cContour)
  hold on
  plot(cContour{i}(:,2),cContour{i}(:,1),'w','Linewidth',1)     
end
axis off
Filename = append(OutputFolder, filesep, CellLine, '_2D_24hour_', 'CellContour.png');
exportgraphics(CellContour, Filename);

gBW = MembraneSegment;
contour = cell(1,size(gBW,3));
for i = 1:size(gBW,3)
    currBW = gBW(:,:,i);
    
    %Get the largest area
    cBWarea = regionprops(currBW,'Area');
    [~,idx2BiggestArea] = max(cell2mat({cBWarea.Area}));
    
    if isempty(idx2BiggestArea)
    else
        
        [pContour] = bwboundaries(currBW);
        for j = 1:length(pContour)
            contour{i}{j} = pContour{j};
        end
    end
end
cContour = contour{1,fr};
MembraneContour = figure();
imagesc(Particles(:,:,fr))
colormap('hot')
axis image
for i = 1:length(cContour)
  hold on
  plot(cContour{i}(:,2),cContour{i}(:,1),'w','Linewidth',1)     
end
axis off
Filename = append(OutputFolder, filesep, CellLine, '_2D_24hour_', 'MembraneContour.png');
exportgraphics(MembraneContour, Filename);