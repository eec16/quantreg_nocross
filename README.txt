README.txt

A MATLAB translation of the R codes accompanying:
Bondell, H. D., Reich, B. J., and Wang, H. (2010).
Non-crossing quantile regression curve estimation. Biometrika 97, 825-838.

Also translates code from:
Base R source code
quantreg R package
sparseM R package

Fortran scripts other than srqfncGate.F come from sources such 
as LAPACK embedded in R or R packages cited above.

Translator: 
Evan Corden

7/3/2021

Add both root and fortran_code folders to MATLAB path or move MEX executable
to root folder

Contents:
Main function:
quantreg_nocross.m

Helper functions:
get_sparseM.m
bandwidth_rq.m
r_quantile.m
sfnMessage.m

/fortran_code folder:
compile_nocross.m- MATLAB code to compile FORTRAN code into MATLAB MEX function.
    Uses the -R2018a option
bckslv.f
boundc.f
chlfct.f
daxpy.f
ddot.f
extract.f
sparskit2.f
srqfnc.f
srqfncGate.F- Gateway function for MATLAB MEX function
srqfncGate.mexa64- pre-compiled MATLAB MEX executable in 64-bit Linux format
    
