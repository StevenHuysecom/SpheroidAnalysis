close all
OutputFolder = "E:\Steven";

%% save membrane figure
f = figure()
set(f, 'Position', get(0, 'Screensize'));
imagesc(rot90(squeeze(Membrane(:, 512, :)), -1));
colormap('gray')
axis image
saveas(f, append(OutputFolder, filesep, 'Membrane.svg'));

%%% save particle figure
g = figure()
set(g, 'Position', get(0, 'Screensize'));
imagesc(rot90(squeeze(Particles(:, 512, :)), -1));
colormap('gray')
axis image
saveas(g, append(OutputFolder, filesep, 'Particles.svg'));

%% save filled segment figure
membrane = Membrane;
membrane(membrane < 10) = 0;
se = strel('cube', 2);
membrane = imdilate(membrane, se);
membrane = medfilt3(membrane, [5 5 5]);
membrane = bwareaopen(membrane, 500000);
for i = 1:size(membrane, 3)
    membrane(:,:,i) = imfill(membrane(:,:,i), "holes");
    stats = regionprops(membrane(:,:,i), 'Area');
    if max(struct2array(stats)) > 5000
        membrane(:,:,i) = bwareaopen(membrane(:,:,i), 5000);
    end
end
membrane = imfill(membrane, "holes");
h = figure()
set(h, 'Position', get(0, 'Screensize'));
bg = rot90(squeeze(Membrane(:, 512, :)), -1);
h1 = imagesc(bg);
colormap('gray');
hold on;
overlay = rot90(squeeze(membrane(:, 512, :)), -1);
yellowOverlay = cat(3, ones(size(overlay)), ones(size(overlay)), zeros(size(overlay)));
j2 = imagesc(yellowOverlay);
set(j2, 'AlphaData', overlay);
axis image;
hold off;
saveas(h, append(OutputFolder, filesep, 'Membrane_filled.svg'));

%% save edge image + with center
SpheroidEdge = edge3(membrane, "sobel", 0.5);
for i = 1:size(membrane, 3)
    Area(i) = sum(membrane(:,:,i), 'all');
end
[~, MaxPlaneIdx] = max(Area);
Coords = [];
if MaxPlaneIdx < size(membrane,3)-20
    RefPlane = MaxPlaneIdx;
else
    RefPlane = size(membrane,3)-20;
end
for i = (RefPlane - 20) : (RefPlane + 20)
    MaxPlane = membrane(:,:,i);
    props = regionprops(MaxPlane, 'Centroid', 'Area');
    [~, largestIdx] = max([props.Area]);
    Coords(end+1,:) = props(largestIdx).Centroid;
end
Center3D = round([mean(Coords, 1), MaxPlaneIdx]);
SegmentPlane = rot90(squeeze(membrane(:,512,:)), -1);
SegmentPlane = bwareaopen(SegmentPlane, 30000);
stats = regionprops(SegmentPlane, "Centroid");
Center = [Center3D(3), round(stats(1).Centroid)];
j = figure();
set(j, 'Position', get(0, 'Screensize'));
bg = rot90(squeeze(Membrane(:, 512, :)), -1);
j1 = imagesc(bg);
colormap('gray');
hold on;
overlay = rot90(squeeze(SpheroidEdge(:, 512, :)), -1);
overlay(Center(1)-2:Center(1)+2, Center(2)-4:Center(2)+4) = 1;
yellowOverlay = cat(3, ones(size(overlay)), ones(size(overlay)), zeros(size(overlay)));
j2 = imagesc(yellowOverlay);
set(j2, 'AlphaData', overlay);
axis image;
hold on;
saveas(j, append(OutputFolder, filesep, 'Membrane_edge.svg'));






%%% particle image with arrow
figure()
imagesc(rot90(squeeze(Particles(:, 512, :)), -1));
colormap('gray')
title('Click on a point in the image');
[x2, y2] = ginput(1);
hold off

z = figure();
set(z, 'Position', get(0, 'Screensize'));
imagesc(rot90(squeeze(Particles(:, 512, :)), -1))
colormap('gray')
hold on 
x1 = Center(2);
y1 = Center(1);
theta = atan2d(y2 - y1, x2 - x1);
quiver(x1, y1, x2 - x1, y2 - y1, 0, 'w', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
r = 25;
theta_vals = linspace(0, deg2rad(theta), 20);
arc_x = x1 + r * cos(theta_vals);
arc_y = y1 + r * sin(theta_vals);
plot(arc_x, arc_y, 'w', 'LineWidth', 1.5);
text(x1 + 1.2 * r, y1 - r , '\theta', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'w');
x_mid1 = (x1 + x2) / 2;
y_mid1 = (y1 + y2) / 2;
dx1 = x2-x1;
dy1 = y2-y1;
magnitude = sqrt(dx1^2 + dy1^2);
dx1 = dx1 / magnitude;
dy1 = dy1 / magnitude;
perp_dx1 = -dy1;  % Rotate by 90 degrees to get perpendicular direction
perp_dy1 = dx1;
offset = 25;
x_text1 = x_mid1 + offset * perp_dx1;
y_text1 = y_mid1 + offset * perp_dy1;
text(x_text1, y_text1, 'r_1', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'w');
axis image;
hold off;
saveas(z, append(OutputFolder, filesep, 'Particle_arrow.svg'));

%%% membrane image with arrow
thetaDeg = deg2rad(theta);
dx = cos(thetaDeg);
dy = sin(thetaDeg);
stepSize = 1;
x = x1;
y = y1;
binaryImage = rot90(squeeze(SpheroidEdge(:, 512, :)), -1);
while x > 1 && x < size(binaryImage, 2) && y > 1 && y < size(binaryImage, 1)
    xi = round(x);
    yi = round(y);

    if binaryImage(yi, xi) > 0
        break;
    end

    x = x + stepSize * dx;
    y = y + stepSize * dy;
end
x3 = xi;
y3 = yi;

q = figure();
set(q, 'Position', get(0, 'Screensize'));
bg = rot90(squeeze(Membrane(:, 512, :)), -1);
q1 = imagesc(bg);
colormap('gray');
hold on;
overlay = rot90(squeeze(SpheroidEdge(:, 512, :)), -1);
overlay(Center(1)-2:Center(1)+2, Center(2)-4:Center(2)+4) = 1;
yellowOverlay = cat(3, ones(size(overlay)), ones(size(overlay)), zeros(size(overlay)));
q2 = imagesc(yellowOverlay);
set(q2, 'AlphaData', overlay);
axis image;
hold on;
x1 = Center(2);
y1 = Center(1);
theta = atan2d(y3 - y1, x3 - x1);
quiver(x1, y1, x3 - x1, y3 - y1, 0, 'w', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
r = 25;
theta_vals = linspace(0, deg2rad(theta), 20);
arc_x = x1 + r * cos(theta_vals);
arc_y = y1 + r * sin(theta_vals);
plot(arc_x, arc_y, 'w', 'LineWidth', 1.5);
text(x1 + 1.2 * r, y1 - r , '\theta', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'w');
x_mid2 = (x1 + x3) / 2;
y_mid2 = (y1 + y3) / 2;
dx2 = x3-x1;
dy2 = y3-y1;
magnitude = sqrt(dx2^2 + dy2^2);
dx2 = dx2 / magnitude;
dy2 = dy2 / magnitude;
perp_dx2 = -dy2;  % Rotate by 90 degrees to get perpendicular direction
perp_dy2 = dx2;
offset = 25;
x_text2 = x_mid2 + offset * perp_dx2;
y_text2 = y_mid2 + offset * perp_dy2;

text(x_text2, y_text2, 'r_2', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'w');
axis image;
hold off;
saveas(q, append(OutputFolder, filesep, 'Spheroid_arrow.svg'));


%%% particle image with arrow and edge
thetaDeg = deg2rad(theta);
dx = cos(thetaDeg);
dy = sin(thetaDeg);
stepSize = 1;
x = x1;
y = y1;
binaryImage = rot90(squeeze(SpheroidEdge(:, 512, :)), -1);
while x > 1 && x < size(binaryImage, 2) && y > 1 && y < size(binaryImage, 1)
    xi = round(x);
    yi = round(y);

    if binaryImage(yi, xi) > 0
        break;
    end

    x = x + stepSize * dx;
    y = y + stepSize * dy;
end
x3 = xi;
y3 = yi;

u = figure();
set(u, 'Position', get(0, 'Screensize'));
bg = rot90(squeeze(Particles(:, 512, :)), -1);
u1 = imagesc(bg);
colormap('gray');
hold on;
overlay = rot90(squeeze(SpheroidEdge(:, 512, :)), -1);
overlay(Center(1)-2:Center(1)+2, Center(2)-4:Center(2)+4) = 1;
yellowOverlay = cat(3, ones(size(overlay)), ones(size(overlay)), zeros(size(overlay)));
u2 = imagesc(yellowOverlay);
set(u2, 'AlphaData', overlay);
axis image;
hold on;
x1 = Center(2);
y1 = Center(1);
theta = atan2d(y3 - y1, x3 - x1);
quiver(x1, y1, x2 - x1, y2 - y1, 0, 'w', 'LineWidth', 1.5, 'MaxHeadSize', 0.5); %Center to Particle
%quiver(x1, y1, x3 - x1, y3 - y1, 0, 'w', 'LineWidth', 1.5, 'MaxHeadSize', 0.5); % Center to Edge
quiver(x3, y3, x2-x3, y2-y3, 0, 'g', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
r = 25;
theta_vals = linspace(0, deg2rad(theta), 20);
arc_x = x1 + r * cos(theta_vals);
arc_y = y1 + r * sin(theta_vals);
plot(arc_x, arc_y, 'w', 'LineWidth', 1.5);
text(x1 + 1.2 * r, y1 - r , '\theta', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'w');
text(x_text1, y_text1, 'r_1', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'w');
x_mid3 = (x2 + x3) / 2;
y_mid3 = (y2 + y3) / 2;
dx3 = x3-x2;
dy3 = y3-y2;
magnitude = sqrt(dx3^2 + dy3^2);
dx3 = dx3 / magnitude;
dy3 = dy3 / magnitude;
perp_dx3 = -dy3;  % Rotate by 90 degrees to get perpendicular direction
perp_dy3 = dx3;
offset = 25;
x_text3 = x_mid3 + offset * perp_dx3;
y_text3 = y_mid3 + offset * perp_dy3;
text(x_text3, y_text3, 'Depth', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'g');
axis image;
hold off;
saveas(u, append(OutputFolder, filesep, 'Particles_arrow_edge.svg'));

%%% save intensity profile
k = figure()
set(k, 'Position', get(0, 'Screensize'));
plot(IntensityInFunctionOfDepth(:,1), IntensityInFunctionOfDepth(:,2))
xlim([-10 100])
xlabel('Penetration depth (µm)')
ylabel('Intensity (a.u)')
saveas(k, append(OutputFolder, filesep, 'Intensity_profile.svg'));