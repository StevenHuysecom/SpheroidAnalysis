t = {'G:\multicolor_polarization\polarisation\20241112_AuBPs\2D_150nm\snaps\0_degrees_1\0_degrees_1_MMStack_Pos0.ome.tif',...
'G:\multicolor_polarization\polarisation\20241112_AuBPs\2D_150nm\snaps\20_degrees_1\20_degrees_1_MMStack_Pos0.ome.tif',...
'G:\multicolor_polarization\polarisation\20241112_AuBPs\2D_150nm\snaps\40_degrees_1\40_degrees_1_MMStack_Pos0.ome.tif',...
'G:\multicolor_polarization\polarisation\20241112_AuBPs\2D_150nm\snaps\60_degrees_1\60_degrees_1_MMStack_Pos0.ome.tif',...
'G:\multicolor_polarization\polarisation\20241112_AuBPs\2D_150nm\snaps\80_degrees_1\80_degrees_1_MMStack_Pos0.ome.tif',...
'G:\multicolor_polarization\polarisation\20241112_AuBPs\2D_150nm\snaps\100_degrees_1\100_degrees_1_MMStack_Pos0.ome.tif',...
'G:\multicolor_polarization\polarisation\20241112_AuBPs\2D_150nm\snaps\120_degrees_1\120_degrees_1_MMStack_Pos0.ome.tif',...
'G:\multicolor_polarization\polarisation\20241112_AuBPs\2D_150nm\snaps\140_degrees_1\140_degrees_1_MMStack_Pos0.ome.tif',...
'G:\multicolor_polarization\polarisation\20241112_AuBPs\2D_150nm\snaps\160_degrees_1\160_degrees_1_MMStack_Pos0.ome.tif',...
'G:\multicolor_polarization\polarisation\20241112_AuBPs\2D_150nm\snaps\180_degrees_1\180_degrees_1_MMStack_Pos0.ome.tif'};

for i = 1:size(t, 2)
    bfI = BioformatsImage(t{1,i});
    Image(:,:,i) = getPlane(bfI, 1,1,1);
end

degrees = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180];

close all
%% Particle 1
for i = 1:10
    ParticleFrame = double(Image(282:293, 1816:1830, i));
    ParticleFrame(ParticleFrame < 1200) = NaN;
    Int(i) = nanmean(ParticleFrame, 'all');
end
figure()
subplot(2,4,1)
plot(degrees, Int)
xlabel('Rotation Lamba/2 polarizer (degrees)');
ylabel('average Px Int in particle spot');
title('Particle 1')
hold on

for i = 1:10
    ParticleFrame = double(Image(883:894, 1796:1810, i));
    ParticleFrame(ParticleFrame < 800) = NaN;
    Int(i) = nanmean(ParticleFrame, 'all');
end
subplot(2,4,5)
plot(degrees, Int)
xlabel('Rotation Lamba/2 polarizer (degrees)');
ylabel('average Px Int in particle spot');
hold on

%% Particle 2
for i = 1:10
    ParticleFrame = double(Image(306:317, 1759:1779, i));
    ParticleFrame(ParticleFrame < 1200) = NaN;
    Int(i) = nanmean(ParticleFrame, 'all');
end
subplot(2,4,2)
plot(degrees, Int)
xlabel('Rotation Lamba/2 polarizer (degrees)');
ylabel('average Px Int in particle spot');
title('Particle 2')
hold on

for i = 1:10
    ParticleFrame = double(Image(306+596:317+596, 1759-25:1779-25, i));
    ParticleFrame(ParticleFrame < 800) = NaN;
    Int(i) = nanmean(ParticleFrame, 'all');
    % figure()
    % imagesc(ParticleFrame)
end
subplot(2,4,6)
plot(degrees, Int)
xlabel('Rotation Lamba/2 polarizer (degrees)');
ylabel('average Px Int in particle spot');
hold on


%% Particle 3
for i = 1:10
    ParticleFrame = double(Image(323:334, 1720:1732, i));
    ParticleFrame(ParticleFrame < 1000) = NaN;
    Int(i) = nanmean(ParticleFrame, 'all');
    
end
subplot(2,4,3)
plot(degrees, Int)
xlabel('Rotation Lamba/2 polarizer (degrees)');
ylabel('average Px Int in particle spot');
title('Particle 3')
hold on

for i = 1:10
    ParticleFrame = double(Image(323+596:334+596, 1720-22:1732-22, i));
    ParticleFrame(ParticleFrame < 600) = NaN;
    Int(i) = nanmean(ParticleFrame, 'all');
end
subplot(2,4,7)
plot(degrees, Int)
xlabel('Rotation Lamba/2 polarizer (degrees)');
ylabel('average Px Int in particle spot');
hold on

%% Background
ImageBg = double(Image(120:599, 1670:1960, :));
ImageBg(ImageBg > 850) = NaN;
for i = 1:10
    Frame = ImageBg(:,:,i);
    Int(i) = nanmean(Frame, 'all');
end
subplot(2,4,4)
plot(degrees, Int)
xlabel('Rotation Lamba/2 polarizer (degrees)');
ylabel('average Px Int in particle spot');
title('background')
hold on

ImageBg = double(Image(730:1180, 1650:1935, :));
ImageBg(ImageBg > 500) = NaN;
for i = 1:10
    Frame = ImageBg(:,:,i);
    Int(i) = nanmean(Frame, 'all');
end
subplot(2,4,8)
plot(degrees, Int)
xlabel('Rotation Lamba/2 polarizer (degrees)');
ylabel('average Px Int in particle spot');
title('background')
hold on
