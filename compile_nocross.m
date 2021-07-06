% compile_nocross.m
%
% Compiles fortran code for quantreg_nocross.m
% Creates srqfncGate.mex Fortran executable file.
%
% Evan Corden 4/6/2020
mex -R2018a srqfncGate.F srqfnc.f extract.f boundc.f chlfct.f cholesky.f sparskit2.f bckslv.f daxpy.f ddot.f