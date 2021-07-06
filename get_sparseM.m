% get_sparseM.m
%
% Returns the matrix.csr format from the sparseM R package of the input mat.
% Only the slots necessary for the quantile regression algorithm are
% returned.
%
% sparseM package by Roger Koenker and Pin Ng.
% https://cran.r-project.org/package=SparseM
%
% Input:
% mat- matrix to be converted to sparseM format
%
% Output:
% ra- see sparseM documentation
% ja- see sparseM documentation
% ia- see sparseM documentation
%
% Evan Corden
% 4/4/2020
function [ra, ja, ia] = get_sparseM(mat)
% Get indices
[row_indices,ja,ra] = find(mat);

% Sort all by row
combined = sortrows([row_indices ja ra], 1);
row_indices = combined(:,1);
ja = combined(:,2);
ra = combined(:,3);

num_rows = size(mat,1);
nnz = length(ja);
ia = zeros(1,num_rows+1);
ia(1) = 1; % first element is always 1

j = 1; % index in row_indices
tmp = row_indices(j);
delta = 0;
for i=1:num_rows
    while tmp==i
        delta = delta + 1;
        
        j = j + 1;
        if j<=nnz
            tmp = row_indices(j);
        else
            tmp = -1;
        end
    end
    
    ia(i+1) = ia(i) + delta;
    delta=0;
end
% Must be horizontal for the Fortran function
ja = ja';
ra = ra';
end

