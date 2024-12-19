classdef OSLocMovie < Core.OSParticleMovie
    %MPLocMovie hold all the method and info necessary for localization,
    %expension of this object will be first centered around Tracking but
    %it will also be easy to extend to STORM/PALM.
    
    properties (SetAccess = 'protected')
       
    end
    
    methods
        
        function obj = OSLocMovie(file,info)
            
            obj  = obj@Core.OSParticleMovie(file,info);
           
        end
  
        function [locPos] = getLocPos(obj,frames)
             %Extract the position of the candidate of a given frame
            [idx] = Core.OrganoidSegmentation.checkFrame(frames, obj.raw.nFrames);
            locPos = obj.corrLocPos{idx};
            
            if isempty(locPos)
                
                warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                
            end
        end
                
        function superResolve(obj,chan)
            disp('super resolving positions ... ');
            if isempty(obj.particles)
                warning('no HIV files found or no particles found, aborting localization');
            else
                
                %Check if some particle were super resolved already:
                [run,SRList] = obj.existZResParticles(obj.info.runMethod,obj.raw.path,'.mat');

                if run
                    data2Resolve = obj.particles.List;
                    nPlanes = obj.raw.nPlanes;
                    nParticles = sum(obj.particles.nParticles);
                    pxSize = obj.info.pxSizeXY;
                    SRList = table(zeros(nParticles,1),...
                            zeros(nParticles,1), zeros(nParticles,1), zeros(nParticles,1),...
                            zeros(nParticles,1), zeros(nParticles,1), zeros(nParticles,1),...
                            zeros(nParticles,1), zeros(nParticles,1), zeros(nParticles,1),...
                            'VariableNames',{'row','col','z','intZ','rowM',...
                            'colM','zM','intensity','SNR','t'});
                    nFrames = length(data2Resolve);
                    h = waitbar(0,'SuperResolving position...');

                    for i = 1:nFrames

                        frameData = data2Resolve{i};
                        frameData2Store = table(zeros(size(frameData)),...
                            zeros(size(frameData)),zeros(size(frameData)),zeros(size(frameData)),...
                            zeros(size(frameData)),zeros(size(frameData)),zeros(size(frameData)),...
                            zeros(size(frameData)),zeros(size(frameData)),zeros(size(frameData)),...
                            'VariableNames',{'row','col','z','intZ','rowM','colM',...
                            'zM','intensity','SNR','t'});

                        fData = obj.getFrame(i,chan); 
                        for j = 1:length(frameData)

                            partData = frameData{j};


                            switch obj.info.zMethod
                                case 'Intensity'
                                    if nPlanes ==1
                                        row  = partData.row(3)*pxSizeXY;
                                        col  = partData.col(3)*pxSizeXY;
                                        z    = partData.z(3);
                                        rowM = partData.row(3)*pxSizeXY;
                                        colM = partData.col(3)*pxSizeXY;
                                        zM   = partData.z(3);
                                        intZ = partData.intensity;
                                        data = table(row,col,z,intZ,rowM,colM,zM,...
                           'VariableNames',{'row','col','z','intZ','rowM','colM','zM'});

                                    else


                                    [data] = obj.resolveXYZInt(partData(:,{'row','col','z','ellip','plane'}),fData);
                                    end
                                otherwise
                                    error('Unknown method for z extraction requested')

                            end

                            frameData2Store(j,{'row','col','z','intZ','rowM','colM','zM'}) = data;
                            frameData2Store.intensity(j) = partData.intensity(3);
                            frameData2Store.SNR(j) = partData.SNR(3);
                            frameData2Store.t(j) = i;

                        end
                    startIdx = find(SRList.row==0,1);   
                    SRList(startIdx:startIdx+height(frameData2Store)-1,:) = frameData2Store;   
                    waitbar(i/nFrames,h,['SuperResolving positions: frame ' num2str(i) '/' num2str(nFrames) ' done']);
                    end
                    close(h);
                    %clean up the list
                    SRList(isnan(SRList.row),:) = [];
                end

                obj.particles.SRList = SRList;
                particle = obj.particles;
                %Save the data
                fileName = sprintf('%s%sparticle.mat',obj.raw.path,'\');
                save(fileName,'particle');
                disp('========> DONE ! <=========');
                
            end
            
        end       
    end
    
    methods (Static)
        
                    
    end
    
    
    methods (Access = protected)
             
        function [data] = resolveXYZInt(obj,partData,frameData)
            planes2Fit =5;
            nPlanes = size(frameData,3);
            pxSize = obj.info.pxSizeXY;
            ROIRad = ceil(obj.info.FWHM_px/2+1);
            planes  = partData(~isnan(partData.plane),:).plane;

            bf = partData.plane(3);
            planePos = obj.raw.movInfo(1).planePos;
            
            %Get ROI XZ, YZ scaled to same pixel size
            [Mag] = Core.HIVLocMovie.getZPhasorMag(partData,ROIRad,frameData);
            
            if abs(bf-nPlanes)<3
                endIdx = nPlanes;
                startIdx = endIdx-(planes2Fit-1);
            elseif bf<3
                startIdx = 1;
                endIdx = planes2Fit;
            else
                startIdx = bf-floor(planes2Fit/2);
                endIdx   = bf+floor(planes2Fit/2);
                
            end

            domain = planePos(startIdx:endIdx);
            data   = [Mag.x]+[Mag.y];
            data   = data(startIdx:endIdx);
            
            guess.sig = 2*obj.info.FWHM_px*pxSize/1000;
            guess.mu  = planePos(bf);
            guess.minMaxDomain =[min(domain) max(domain)];
            [Res,fit] = SimpleFitting.gauss1D(data,domain,guess);

            z = Res(2);
            %if the z position is out of bound we do not consider the data
            if or(z<min(domain),z>max(domain))
                z   = NaN;                           
                row = NaN;
                col = NaN;
                zM   = NaN;                           
                rowM = NaN;
                colM = NaN;
                intZ = NaN;
            else
              
                row = partData.row(3)*pxSizeXY;
                col = partData.col(3)*pxSizeXY;
                zM = z/mean(diff(obj.raw.movInfo(1).planePos));                      
                rowM = partData.row(3);
                colM = partData.col(3);
                intZ = fit(3);

            end
           
            %store the data
            data = table(row,col,z,intZ,rowM,colM,zM,...
                   'VariableNames',{'row','col','z','intZ','rowM','colM','zM'});
        end
    end
    
    methods (Static)
        
        function [run, SRList] = existZResParticles(runMethod,Path, ext)
            SRList = [];
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.OrganoidSegmentation.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'particle')==true)
                        particle = load([file2Analyze(1).folder filesep 'particle.mat']);
                        particle = particle.particle;
                        if isfield(particle,'SRList')
                            if ~isempty(particle.SRList)
                                run = false;
                                SRList = particle.SRList;
                            else
                                run = true;
                            end
                        else
                            run = true;
                        end
                    else
                
                        run = true;
                        
                    end
                    
                case 'run'
                    
                     run = true;
   
            end
        end
        function [Mag] = getZPhasorMag(partData,ROIRad,volIm)

            %Possible improvement : Translate the coordinate of the best
            %focus into the otherplanes to extract the exact value where
            %the PSF should be taken    
            imSize = size(volIm);
            pos = [round(nanmean(partData.row)),round(nanmean(partData.col))];

            ROIs = Misc.getROIs(pos,ROIRad,imSize(1:2));

            ROI = volIm(ROIs(1):ROIs(2),ROIs(3):ROIs(4),:);

            Mag = struct('x',zeros(1,size(ROI,3)),'y',zeros(1,size(ROI,3)));
            for i =1:size(ROI,3)
                [~,~,~,Mag(i).x,Mag(i).y] = Localization.phasor(ROI(:,:,i));
            end

        end 

    end
end