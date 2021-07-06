% sfnMessage.m
%
% Translated from sfnMessage.r in quantreg R package 
% version 5.5 by Roger Koenker and coauthors
%
% https://cran.r-project.org/package=quantreg
%
% April 1, 2020
%
% Evan Corden
function msg = sfnMessage(ierr)
switch ierr
    case 1
        msg = "insufficient storage (work space) when calling extract";
    case 2
        msg = "nnzd > nnzdmax";
    case 3
        msg = "insufficient storage in iwork when calling ordmmd";
    case 4
        msg = "insufficient storage in iwork when calling sfinit";
    case 5
        msg = "nnzl > nnzlmax when calling sfinit";
    case 6
        msg = "nsub > nsubmax when calling sfinit";
    case 7
        msg = "insufficient work space in iwork when calling symfct";
    case 8
        msg = "inconsistancy in input when calling symfct";
    case 9
        msg = "tmpsiz > tmpmax when calling bfinit; increase tmpmax";
    case 10
        msg = "nonpositive diagonal encountered blkfct() matrix is not positive definite";
    case 11
        msg = "insufficient work storage in tmpvec when calling blkfct";
    case 12
        msg = "insufficient work storage in iwork  when calling blkfct";
    case 13
        msg = " nnzd > nnzdmax in e,je when calling amub";
    case 14
        msg = "nnzd > nnzdmax in g,jg when calling amub";
    case 15
        msg = "nnzd > nnzdmax in h,jh when calling aplb";
    case 17
        msg = "tiny diagonals replaced with Inf when calling blkfct";
    otherwise
        msg = "impossible error condition";
end