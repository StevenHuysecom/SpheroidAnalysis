classdef LysosomeSegmentation < handle
    %LYSOSOMESEGMENTATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        raw
        info
        channels
        results
    end
    
    methods
        function obj = LysosomeSegmentation(raw,info)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.raw = raw;
            obj.info = info;
        end
        
        function set.raw(obj,raw)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            assert(isstruct(raw),'raw is expected to be a structure');
            assert(and(isfield(raw,'path'),isfield(raw,'ext')),'raw is expected to be a structure with 2 fields: path and ext');
            assert(ischar(raw.path), 'Input path needs to be a char or string');
            assert(isfolder(raw.path),'Input path should be a folder containing one dataset to analyze');
            obj.raw = raw;
        end

        function loadData(obj,chan)
            path = obj.raw.path;
            ext  = obj.raw.ext;
            %let us check that there is no channel data existing
            if ~obj.existChannel
                disp('no channel data found, starting extraction ...')
                %get all file of appropriate extension in the file
                fileList = Core.OrganoidSegmentation.getFileInPath(path,ext);
                %get the different channel from the data
                channel  = obj.retrieveChannel2(fileList,chan);       
                %save the channel as matlab variable for the future
                filename = [obj.raw.path filesep 'channels.mat'];
                save(filename,'channel','-v7.3');
                disp('==========> DONE <==========')
            else
                disp('channel data found, loading from existing file ...')
                filename = [path filesep 'channels.mat'];
                tmp = load(filename);
                field = fieldnames(tmp);
                channel = tmp.(field{1});
                disp('==========> DONE <==========')
            end
            nFrames = size(channel.membrane,3);
            obj.raw.nFrames = nFrames;
            obj.raw.nPlanes = nFrames;
            %store the data back into the object
            obj.channels = channel;
        end

        function loadDataBioform(obj, chan)

            FieldNames = fieldnames(chan);
            for l = 1:size(FieldNames, 1)
                if exist(append(obj.raw.path, filesep, chan.(FieldNames{l,1}), '.mat'))
                    S = load(append(obj.raw.path, filesep, chan.(FieldNames{l,1}), '.mat'));
                    Name = fieldnames(S);
                    obj.channels.(chan.(FieldNames{l,1})) = S.(Name{1,1});
                end
            end
        end
        
        function channels = getChannel(obj,chan)
            switch nargin 
                case 1
                    channels = obj.channels;
                case 2
                    channels = obj.channels.(chan);
            end
        end

        function showChannel(obj)
            channel = obj.getChannel;
            field = fieldnames(channel);
            
            nField = length(field);
            sliceToShow = round(size(channel.(field{1}),3)/2);
            frameToShow = round(size(channel.(field{1}),4)/2);
            figure
            for i = 1:nField
                subplot(1,nField,i)
                currChan = channel.(field{i});
                try
                    imagesc(currChan(:,:,sliceToShow,frameToShow))
                catch
                    
                end
                axis image
                colormap('hot');
                title(field{i})
            end
        end
    end
end

