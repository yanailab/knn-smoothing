% K-nearest neighbor smoothing for UMI-filtered scRNA-Seq data
% (Matlab implementation)

% Author: Maayan Baron <Maayan.Baron@nyumc.org>
% Copyright (c) 2017 New York University

function [mat_smooth] = knn_smooth(raw_mat,k)
%This function smooth computes the a smooth matrix using K-nearest neighbor
%by Wagner et al. 2017. The output is a smoothed matrix, not normilazed or
%transformed. 
%                           INPUT ARGUMENTS
%                           ---------------
%   [mat_smooth] = knn_smooth(raw_mat,k) computes the smoothed expression of
%   raw_mat where rows are genes/features and columns are samples/cells using
%   k+1 neighbours.
%
%                         OUTPUT ARGUMENTS
%                         ----------------
%   The results of the is a the smoothed matrix mat_smooth.


mat_smooth = raw_mat;
num_of_steps = ceil(log2(k+1));
disp_text = ['number of steps: ' num2str(num_of_steps)];
disp(disp_text)
for s = 1:num_of_steps
    k_step = min(2^s-1,k);
    mat_tpm = median(sum(mat_smooth))*bsxfun(@rdivide,mat_smooth,sum(mat_smooth));
    mat_trans = sqrt(mat_tpm)+sqrt(mat_tpm+1);
    dist_matrix = squareform(pdist(mat_trans'));
    disp_texp = ['calculating distance matrix ' num2str(s) '/' num2str(num_of_steps) ' times'];
    disp(disp_texp)
    
    %mat_smooth = [];
    for cell = 1:size(mat_smooth,2)
        [~,c_sort_cell_idx] = sort(dist_matrix(cell,:));
        mat_smooth(:,cell) = sum(raw_mat(:,c_sort_cell_idx(1:k_step+1)),2);
    end
end
disp done!
end
