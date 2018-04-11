% K-nearest neighbor smoothing for high-throughput scRNA-Seq data
% (Matlab implementation.)

% Author: Maayan Baron <Maayan.Baron@nyumc.org>
% Copyright (c) 2017, 2018 New York University

function [mat_smooth] = knn_smooth(raw_mat,k,varargin)
% This function smooth computes the smoothed matrix using K-nearest neighbor
% by Wagner et al. (2018). The output is a smoothed matrix, not normalized or
% transformed. 
% Dependencies: Randomized Singular Value Decomposition (rsvd) function: 
% https://www.mathworks.com/matlabcentral/fileexchange/47835-randomized-singular-value-decomposition
%
%                           INPUT ARGUMENTS
%                           ---------------
%   [mat_smooth] = knn_smooth(raw_mat,k) computes the smoothed expression of
%   raw_mat where rows are genes/features and columns are samples/cells using
%   k neighbours.
%   varargin:
%   'num_of_pc' (default = 10) perform PCA before computing distances
%      
%                         OUTPUT ARGUMENTS
%                         ----------------
%   The smoothed expression matrix.

%default
num_of_pc = 10;

% only want 1 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 1
    error('Requires at most 1 optional inputs');
end

% set defaults for optional inputs
optargs = num_of_pc;

% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs = varargin;

% Place optional args in memorable variable names
[num_of_pc] = optargs{:};

mat_smooth = raw_mat;
num_of_steps = ceil(log2(k+1));
disp_text = ['number of steps: ' num2str(num_of_steps)];
disp(disp_text)
for s = 1:num_of_steps
    k_step = min(2^s-1,k);
    mat_tpm = median(sum(mat_smooth))*bsxfun(@rdivide,mat_smooth,sum(mat_smooth));
    mat_trans = sqrt(mat_tpm)+sqrt(mat_tpm+1);
    [~,genes_with_exp] = sort(sum(mat_smooth'),'descend');
    disp_texp = ['preforming pca ' num2str(s) '/' num2str(num_of_steps) ' times'];
    disp(disp_texp)
    [~,~,V] = rsvd(mat_trans',num_of_pc);
    score = mat_trans'*V;
    disp_texp = ['preforming pca - done! ' num2str(s) '/' num2str(num_of_steps) ' times'];
    disp(disp_texp)
    dist_matrix = squareform(pdist(score));
    disp_texp = ['calculating distance matrix ' num2str(s) '/' num2str(num_of_steps) ' times'];
    disp(disp_texp)
    for cell = 1:size(mat_smooth,2)
        [~,c_sort_cell_idx] = sort(dist_matrix(cell,:));
        mat_smooth(:,cell) = sum(raw_mat(:,c_sort_cell_idx(1:k_step+1)),2);
    end
end
disp done!
end