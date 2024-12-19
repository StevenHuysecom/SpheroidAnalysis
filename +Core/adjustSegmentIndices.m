function new_z_stack = adjustSegmentIndices(z_stack,idx)
        [rows, cols, num_slices] = size(z_stack);
        new_z_stack = z_stack;

        %% go up from idx
        for k = idx:num_slices
            current_slice = new_z_stack(:,:,k);
            prev_slice = new_z_stack(:,:,k-1);
            adj_slice = zeros(size(current_slice));
    
            segment_indices = unique(current_slice(:));
            segment_indices(segment_indices == 0) = [];
    
            for i = 1:length(segment_indices)
                segment_index = segment_indices(i);
                segment_mask = (current_slice == segment_index);
    
                overlap = prev_slice(segment_mask);
                
                previous_indices = overlap(:);
                previous_indices(previous_indices == 0) = [];
                if isempty(previous_indices)
                    neighbour_idx = max(prev_slice(:))+1;
                else
                    neighbour_idx = mode(previous_indices);
                end
                adj_slice(segment_mask == 1) = neighbour_idx;

            end
            new_z_stack(:,:,k) = adj_slice;
        end

        %% go down from idx
        for k = idx-1:-1:1
            current_slice = new_z_stack(:,:,k);
            prev_slice = new_z_stack(:,:,k+1);
            adj_slice = zeros(size(current_slice));
    
            segment_indices = unique(current_slice(:));
            segment_indices(segment_indices == 0) = [];
    
            for i = 1:length(segment_indices)
                segment_index = segment_indices(i);
                segment_mask = (current_slice == segment_index);
    
                overlap = prev_slice(segment_mask);
                
                previous_indices = overlap(:);
                previous_indices(previous_indices == 0) = [];
                if isempty(previous_indices)
                    neighbour_idx = max(prev_slice(:))+1;
                else
                    neighbour_idx = mode(previous_indices);
                end
                adj_slice(segment_mask == 1) = neighbour_idx;

            end
            new_z_stack(:,:,k) = adj_slice;
        end
end

