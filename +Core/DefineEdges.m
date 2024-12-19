function [ZStackEdges] = DefineEdges(A)
    
    kernel = ones(3, 3, 3);
    kernel(2, 2, 2) = 0;
    
    ZStackEdges = A; % Start with the original matrix
    
    unique_vals = unique(A);
    for val = unique_vals'
        bin_matrix = (A == val);
        
        neighbor_sum = convn(bin_matrix, kernel, 'same');
        surrounded = (neighbor_sum == 26);
        
        % Set pixels not fully surrounded to zero
        ZStackEdges(~surrounded & A == val) = 0;
    end
end

