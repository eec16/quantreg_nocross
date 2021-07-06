% Translation of the bandwidth.rq function from quantreg 
% R package Version 5.85. 
%
% https://cran.r-project.org/package=quantreg
% 
% R code is commented above MATLAB code
% 3/4/2021
% Translated by Evan Corden
function bw = bandwidth_rq(p, n, hs, alpha)
%bandwidth.rq <- function(p, n, hs = TRUE, alpha = 0.05)
% Handle default parameter values
if nargin < 3
    hs = 1;
end

if nargin < 4
    alpha = 0.05;
end

% {
% 
% 	# Bandwidth selection for sparsity estimation two flavors:
% 
% 	#	Hall and Sheather(1988, JRSS(B)) rate = O(n^{-1/3})
% 
% 	#	Bofinger (1975, Aus. J. Stat)  -- rate = O(n^{-1/5})
% 
% 	# Generally speaking, default method, hs=TRUE is preferred.

% 	x0 <- qnorm(p)
x0 = norminv(p);
% 
% 	f0 <- dnorm(x0)
f0 = normpdf(x0);
% 
% 	if(hs)
% 
%             n^(-1/3) * qnorm(1 - alpha/2)^(2/3) *
% 
%                 ((1.5 * f0^2)/(2 * x0^2 + 1))^(1/3)
% 
% 	else n^-0.2 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^ 0.2
% 
% }

if hs
    bw = n^(-1/3) * norminv(1-alpha/2)^(2/3)*((1.5 * f0^2)/(2 * x0^2 + 1))^(1/3);
else
    bw = n^-0.2 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^ 0.2;
end