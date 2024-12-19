%% main
%This code aims to be used by the end-user for analysis of HIV particles in
%cells as well as segmentation of various cell structures. This code will
%mostly call the main functions needed while most of the work will occur in
%the background. This is to make this code as understandable and short as
%possible
clear
close all
clc

%% User Input
file.path = 'D:\Documents\Unif\PostDoc\2024 - Data\05 - May\Maria 2\dataTest\NC'; %give empty brackets [], to open file selection
file.ext = '.tif'; %expected extension of the movie(s);
info.runMethod = 'load'; %'load'or 'run', if load is chosen it will try to load previously calculated data(e.g localized particles)
info.fitMethod = 'phasor'; %'Gauss' or 'phasor'
info.zMethod   = 'Intensity';
%if it exist. run will always re-run the analysis and erase previous data.

info.pxSizeXY = 454.5; 
info.pxSizeZ  = 1000;
info.FWHM_px = 3; %in pixel
%Give info about the channels, the word needs to be lowercase with no typos
%care that the
chan.ch00 = 'Nuclei';
chan.ch02 = 'Transfection';
chan.ch03 = 'ignore';
chan.ch04 = 'ignore';

%% Data Loading

OrgMov = Core.OSLocMovie(file,info);
OrgMov.loadData(chan);


%% Particle detection
%rough detection of particle
detectParam.delta = 6; % ROI around the molecule for testing
detectParam.chi2  = 60;% Threshold for detections ([24-80],24 for single molecules, up to 80 for brighter objects
OrgMov.findCandidatePos(detectParam,'cell');
if ~isempty(OrgMov.candidatePos)
    OrgMov.showCandidate(1,5,'cell');
end
%% Localization
%Fitting for accurate localization of the detected particle
OrgMov.SRLocalizeCandidate('cell');

%% consolidation of particle position
%check that particle are detected in more than one z-slice
OrgMov.consolidatePlanes();

%% Super-resolve in 3D
%need to update the position of the plane
OrgMov.superResolve('cell');


%% Segmentation of Lamina
membrane = 'lamina';
sensitivity = 0.4;
OrgMov.segmentLamina(sensitivity,true);
OrgMov.showMembrane(membrane);
%% Segmentation of NUP
sensitivity = 0.4;%put higher when signal is worst (e.g. expansion images)
OrgMov.segmentNUP(sensitivity,true);
membrane = 'NUP';
OrgMov.showMembrane(membrane);
%% Segmentation of Lipid
membrane = 'lipid';
%HIVData.segmentLipid();
%HIVData.showMembrane(membrane);
%% Segmentation of red lipid
OrgMov.segmentRedLipid();
OrgMov.showMembrane('lipid',1);

%% fit lipid
OrgMov.getMembranePos('lipid');
OrgMov.showMembrane('lipid');
%% get membrane position 
membrane = 'lamina';
OrgMov.getMembranePos(membrane);
OrgMov.showMembrane(membrane,1);
% get second membrane position
membrane = 'NUP';
OrgMov.getMembranePos(membrane);
OrgMov.showMembrane(membrane,2);

%%
OrgMov.showAllMembranes;

%% Get Distance between membranes
distance = OrgMov.getMembraneToMembraneDistance('lamina','NUP');

%% Get distance between particles and membrane

a = OrgMov.getHIVToMembraneDistance('NUP');




