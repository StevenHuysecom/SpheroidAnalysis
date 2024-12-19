function [imStack,simPos,simElip] = simulateImages(nImages,emitter,detector,bkg)

% Information about detector
x_size = detector.xSize;
pix_size = detector.pxSize; %[nm/pix]

% information about emitters
em_n = emitter.num;
em_mean_int = emitter.meanInt;
em_int_sigma = emitter.intSigma;
FWHM_nm = emitter.FWHM_nm;
posRange = emitter.posRange;
maxEllip = emitter.maxSigma;
% information about normal bg
mean_bg = bkg.mean;
SNR     = bkg.SNR;
d_range   = 'uint16';

imStack = zeros(x_size,x_size,nImages);
simPos  = zeros(emitter.num,2,nImages);
simElip = zeros(emitter.num,nImages);

if nImages > 10
h = waitbar(0, 'Simulations of images...');
end
for i = 1:nImages
    
% calculations for emitters positions and int
[ em_pos ] = EmitterSim.getRandPos(posRange(2)-posRange(1), em_n );
em_pos = em_pos+posRange(1);
[simPos(1:em_n,1,i),ind] = sort(em_pos(:,1));
simPos(1:em_n,2,i) = em_pos(ind,2);

int_model.name     = 'normal';
int_model.mean_int = em_mean_int;
int_model.sigma    = em_int_sigma;
[ em_int ] = EmitterSim.getIntensity( int_model, em_n );

% small calculations for image generation
y_size = x_size;
xv    = 1:1:x_size;
yv    = 1:1:y_size;
[X,Y] = meshgrid(xv,yv);
im = uint16(zeros(size(X)));

% calculations for psf
FWHM_pix = FWHM_nm / pix_size; %[pix]
sigma_pix = FWHM_pix / (2*((2*log(2))^0.5));
[sigX,sigY] = EmitterSim.getRandSigma(sigma_pix, maxEllip, em_n);
simElip(1:em_n,i) = sigY./sigX;
psf_model.name = 'gaussian';
psf_model.sigma_x = sigma_pix;
psf_model.sigma_y = psf_model.sigma_x;

[ G ] = EmitterSim.getPSF( X, Y, round(mean(xv)), round(mean(yv)), psf_model);
G = G.*(em_mean_int);
max_int = round(max(G(:)));

%Generate the noise of the background in the same way Boris was generating
%the noise in test2DGradTracking so the results are comparable.
switch emitter.noiseType
    case 'Gaussian'
        noiseProp.maxCount = max_int;
        noiseProp.S2N = SNR;
        noiseProp.bkg = mean_bg;
        ROI = zeros(x_size);
        bg_im = ImageNoise.generateNoise(ROI,'Gaussian',noiseProp);
       
    case 'Poisson'
        noiseProp.maxCount = max_int;
        noiseProp.S2N = SNR;
        noiseProp.bkg = mean_bg;
        ROI = zeros(x_size);
        bg_im = ImageNoise.generateNoise(ROI,'Poisson',noiseProp);
        
    otherwise
        bg_im = ones(size(im))*mean_bg;
end

bg_im = uint16(bg_im);

for em_i = 1:em_n
    % get possition and intensity of the emitter
    x0_pix = em_pos(em_i,1);
    y0_pix = em_pos(em_i,2);
    n_counts = em_int(em_i);
    psf_model.sigma_x = sigX(em_i);
    psf_model.sigma_y = sigY(em_i);
    % get psf pdf of emitter
    [ psf ] = EmitterSim.getPSF( X, Y, x0_pix, y0_pix, psf_model );
    % get out only area of interest
    [ roi_lims ] = EmitterSim.getROI(x0_pix, y0_pix, 20, x_size, y_size);
    psf_roi = psf(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2));
    % sample psf
    [ im_roi ] = EmitterSim.samplePSF( psf_roi, n_counts, false );
    % update total image
    im(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2)) = ...
              im(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2)) + im_roi;
end
ccd_frame  = im+bg_im;
imStack(:,:,i) = ccd_frame;
if nImages > 10
waitbar( i/nImages, h, sprintf('Simulations of images... - %d / %d percent achieved',round(i/nImages*100),100));
end
end

if nImages > 10
close(h)
end

end
