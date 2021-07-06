% quantreg_nocross.m
%
% A MATLAB translation of the R codes accompanying:
% Bondell, H. D., Reich, B. J., and Wang, H. (2010).
% Non-crossing quantile regression curve estimation. Biometrika 97, 825-838.
%
% Input:
% y - dependent data set
% x - independent data set. This code automatically adds a constant to
% this data set
% taus - vector of quantiles between 0 and 1
%
% Output:
% bhat - quantile regression coefficent matrix of num variables X num
% quantiles
% se_bhat - matrix of standard errors of the regression coefficients
%
% Original R code is commented above each line of MATLAB code
%
% Evan Corden 3/5/2021
% Quantile Regression with non-crossing condition
function [bhat,se_bhat,cov_bhat,taus]=quantreg_nocross(y, x, taus)
%
% Input validation checks
%   if (length(taus)<2)
%	  stop("At least 2 quantile levels should be specified. If only using a single quantile, you should use rq.")
%
%    taus = sort(taus)
%    if (max(taus)>=1)
%	  stop("All quantile levels should be between 0 and 1, not including the boundary.")
%    if (min(taus)<=0)
%	  stop("All quantile levels should be between 0 and 1, not including the boundary.")
if (length(taus)<2)
    disp("At least 2 quantile levels should be specified. If only using a single quantile, you should use rq.")
    return
end

taus = sort(taus);
if (max(taus)>=1)
    disp("All quantile levels should be between 0 and 1, not including the boundary.")
    return
end
if (min(taus)<=0)
    disp("All quantile levels should be between 0 and 1, not including the boundary.")
    return
end

n = size(x,1); % Number of observations
p = size(x,2); % Number of explanatory variables
m = length(taus); % Number of quantiles
pp = p+1;

% Y = rep(y, m)
Y = y;
for i =1:m-1
    Y = [Y; y];
end

% Y = rep(y, m)
xtemp = x;

% x = matrix(0,nrow=n,ncol=p)
x = zeros(n, p);
% shifts = apply(xtemp,2,min)
shifts = min(xtemp,[], 1);
% scalings = rep(0,p)
scalings = zeros(1, p);

%     for (i in 1:p)
%     {
%      	x[,i] = xtemp[,i] - shifts[i]
% 	scalings[i] = max(x[,i])
% 	x[,i] = x[,i]/scalings[i]
%     }
for i=1:p
    x(:,i) = xtemp(:,i) - shifts(i);
    scalings(i) = max(x(:,i));
    x(:,i) = x(:,i)/scalings(i);
end

%x = cbind(rep(1,n),x)
x = [ones(n,1) x];

%D = diag(m)
%D[lower.tri(D)] = 1
D = ones(m);
D = tril(D);

%X = D %x% x
X = kron(D,x);

%X2 = cbind(X, -X)
X2 = [X -X];

% sX = as.matrix.csr(X2)
sX = sparse(X2);

%K = m*pp
K = m*pp;

%R1 = (diag(m) %x% rep(1,1) %x% t(c(1, rep(0, p))))[-(1:1),]
R1 = kron(kron(eye(m),1), [1 zeros(1,p)]);
R1 = R1(2:end,:);

%R2 = rep(1,1)
R2 = 1;

%     for (j in 1:p)
%     {
%     	R2 = cbind(R2, rep(1,1) %x% (diag(1) %x% rep(1,1)))
%     }
for j=1:p
    R2 = [R2 kron(1, kron(1, 1))];
end

%R2 = (diag(m) %x% R2)[-(1:1),]
R2 = kron(eye(m),R2);
R2 = R2(2:end,:);

%sR = as.matrix.csr(rbind(diag(2*K), cbind(R1, -R2)))
sR = sparse([eye(2*K); [R1 -R2]]);

%r2 = rep(0, 2*K + (m-1)*(1))
r2 = zeros(1,2*K + (m-1));

%tau.v = rep(taus, each=n)
tauv = repelem(taus,n);

% rhs = t(X2)%*%(1-tau.v)
rhs = X2' * (1-tauv)';

% coeff1 =  myrq.fit.sfnc(sX, Y, sR, r2, tau=tau.v, rhs=rhs,tmpmax=100000)$coef
tmpmax = 100000;
[coeff1, ierr, maxiter, time] = sfnc(sX, Y, sR, r2, tauv, rhs, tmpmax);

% coeff = coeff1[1:K]-coeff1[-(1:K)]
coeff = coeff1(1:K) - coeff1(K+1:end);
% gamma.m = matrix(coeff, ncol=m)
gamma_m = reshape(coeff, [], m);
% D = diag(m)
% D[upper.tri(D)]=1
D = ones(m);
D = triu(D);

% bhat.temp = gamma.m %*% D
bhat_temp = gamma_m * D;

% cov.bhat = array(0,c(pp,pp,m))
cov_bhat = zeros(pp,pp,m);

% se.bhat = matrix(0,ncol=m,nrow=pp)
se_bhat = zeros(pp,m);

% transform.mat = ...
% rbind(c(1,-shifts/scalings),cbind(rep(0,p),diag(as.vector(1/scalings),...
% nrow = length(as.vector(1/scalings)), ncol = length(as.vector(1/scalings)))))
base = [zeros(p,1) diag(1./scalings)];
top = [1 -shifts./scalings];
transform_mat = [top; base];
%     bhat = transform.mat%*%bhat.temp
bhat = transform_mat*bhat_temp;

%return;
%     for (j in 1:m)
%     {
% 	cov.bhat[,,j] = se.constr(x,y,bhat.temp[,j],taus[j])
%    	cov.bhat[,,j] = transform.mat%*%cov.bhat[,,j]%*%t(transform.mat)
% 	se.bhat[,j] = sqrt(diag(cov.bhat[,,j]))
%     }
for j=1:m
    cov_bhat(:,:,j) = se_constr(x,y,bhat_temp(:,j),taus(j));
    cov_bhat(:,:,j) = transform_mat*cov_bhat(:,:,j)*transform_mat';
    se_bhat(:,j) = sqrt(diag(cov_bhat(:,:,j)));
end
%
%     vars<-c("intercept",paste("x",1:p,sep=""))
%     rownames(bhat)<-rownames(se.bhat)<-vars
%     colnames(bhat)<-colnames(se.bhat)<-taus
%     dimnames(cov.bhat)<-list(vars,vars,taus)
%
%     constr.qr.fit = NULL
%     constr.qr.fit$bhat = bhat
%     constr.qr.fit$se.bhat = se.bhat
%     constr.qr.fit$cov.bhat = cov.bhat
%     constr.qr.fit$taus = taus
%
%     return(constr.qr.fit)
end

% myrq.fit.sfnc <- function (x, y, R, r, tau, rhs, nsubmax, tmpmax,nnzlmax,...
%   cachsz = 64,small = 1e-08, maxiter = 100, warn.mesg = TRUE)
function [coef, ierr, maxiter, time] = sfnc(x, y, R, r, tau, rhs, tmpmax)
% constants taken from R-code arguments
cachsz = 64;
small = 1e-08;
maxiter = 100;

%y <- -y
y = -y;

%r <- -r
r = -r;

%n1 <- length(y)
n1 = length(y);

%m <- x@dimension[2]
m = size(x,2);

if (n1 ~= size(x,1))
    disp("The design matrix A1' and response vector y are not compatible")
    return
end
% n2 <- length(r)
n2 = length(r);

%     if (n2 != R@dimension[1])
%         stop("The constraint matrix A2' and the constraint right-hand-side are not compatible")
if n2 ~= size(R,1)
    disp('The constraint matrix A2 and the constraint right-hand-side are not compatible')
    return
end

%maxn1n2 <- max(n1, n2)
maxn1n2 = max([n1 n2]);

%u <- rep(1, length = n1)
u = ones(1, n1);

%     if(length(tau)==1)     x1 <- rep(1 - tau, length = n1)
%     else x1=1-tau
if (length(tau)==1)
    x1 = 1 - tau;
    for i=1:n1-1
        x1 = [x1 (1-tau)];
    end
else
    x1 = 1 - tau;
end

%x2 <- rep(1, length = n2)
x2 = ones(1, n2);

% wwm <- vector("numeric", 6 * m)
wwm = zeros(1, 6*m);

%wwm[1:m] <- rhs
wwm(1:m) = rhs;

%nnzx <- x@ia[x@dimension[1] + 1] - 1
nnzx = nnz(x);

%nnzR <- R@ia[R@dimension[1] + 1] - 1
nnzR = nnz(R);

% nnzdmax <- max(nnzx, nnzR)
nnzdmax = max([nnzx nnzR]);

%iwmax <- 7 * m + 3
iwmax = 7*m + 3;

%ao1 <- t(x)
ao1 = x';

%ao2 <- t(R)
ao2 = R';

% e <- ao1 %*% x
e = ao1 * x;

%g <- ao2 %*% R
g = ao2 * R;

%h <- e + g
h = e + g;

%nnzemax <- e@ia[e@dimension[1] + 1] - 1
nnzemax = nnz(e);

%nnzgmax <- g@ia[g@dimension[1] + 1] - 1
nnzgmax = nnz(g);

%nnzhmax <- h@ia[h@dimension[1] + 1] - 1
nnzhmax = nnz(h);

nnzlmax = 4 * nnzdmax;

nsubmax = nnzhmax;

%s <- u - x1
s = u - x1;

%chol.o <- chol(e, tmpmax = tmpmax, nsubmax = nsubmax, nnzlmax = nnzlmax)
%b <- backsolve(chol.o, ao1 %*% y)
% See page 26 of sparseM documentation
b = pinv(full(e)) * (ao1*y);

%r1 <- y - x %*% b
r1 = y - x*b;

%z1 <- ifelse(abs(r1) < small, (r1 * (r1 > 0) + small), r1 * (r1 > 0))
if abs(r1) < small
    z1 = r1 .* (r1 > 0) + small;
else
    z1 = r1 .* (r1 > 0);
end

%w <- z1 - r1
w = z1 - r1;

%z2 <- rep(1, n2)
z2 = ones(1,n2);

%wwn1 <- matrix(0, n1, 10)
wwn1 = zeros(n1,10);

%wwn1[, 1] <- z1
wwn1(:,1) = z1;

%wwn1[, 2] <- w
wwn1(:,2) = w;

%wwn2 <- matrix(0, n2, 7)
wwn2 = zeros(n2,7);

%wwn2[, 2] <- z2
wwn2(:,2) = z2;

%     srqfnc.o <- .Fortran("srqfnc", n1 = as.integer(n1), m = as.integer(m),
%         nnzx = as.integer(nnzx), x = as.double(x@ra), jx = as.integer(x@ja),
%         ix = as.integer(x@ia), ao1 = as.double(ao1@ra), jao1 = as.integer(ao1@ja),
%         iao1 = as.integer(ao1@ia), n2 = as.integer(n2), nnzR = as.integer(nnzR),
%         R = as.double(R@ra), jR = as.integer(R@ja), iR = as.integer(R@ia),
%         ao2 = as.double(ao2@ra), jao2 = as.integer(ao2@ja), iao2 = as.integer(ao2@ia),
%         nnzdmax = as.integer(nnzdmax), d = double(nnzdmax), jd = integer(nnzdmax),
%         id = integer(m + 1), dsub = double(nnzhmax + 1), jdsub = integer(nnzhmax +
%             1), nnzemax = as.integer(nnzemax), e = as.double(e@ra),
%         je = as.integer(e@ja), ie = as.integer(e@ia), nnzgmax = as.integer(nnzgmax),
%         g = double(nnzgmax), jg = integer(nnzgmax), ig = integer(m +
%             1), nnzhmax = as.integer(nnzhmax), h = double(nnzhmax),
%         jh = integer(nnzhmax), ih = integer(m + 1), nsubmax = as.integer(nsubmax),
%         lindx = integer(nsubmax), xlindx = integer(m + 1), nnzlmax = as.integer(nnzlmax),
%         lnz = double(nnzlmax), xlnz = integer(m + 1), iw = integer(m *
%             5), iwmax = as.integer(iwmax), iwork = integer(iwmax),
%         xsuper = integer(m + 1), tmpmax = as.integer(tmpmax),
%         tmpvec = double(tmpmax), maxn1n2 = as.integer(maxn1n2),
%         ww1 = double(maxn1n2), wwm = as.double(wwm), wwn1 = as.double(wwn1),
%         wwn2 = as.double(wwn2), cachsz = as.integer(cachsz),
%         level = as.integer(8), x1 = as.double(x1), x2 = as.double(x2),
%         s = as.double(s), u = as.double(u), c1 = as.double(y),
%         c2 = as.double(r), sol = as.double(b), small = as.double(small),
%         ierr = integer(1), maxiter = as.integer(maxiter), time = double(7),
%         PACKAGE = "quantreg")[c("sol", "ierr", "maxiter", "time")]
% Get sparseM formats for all relevant sparse matrices
[x_ra, x_ja, x_ia] = get_sparseM(x);
[ao1_ra, ao1_ja, ao1_ia] = get_sparseM(ao1);
[R_ra, R_ja, R_ia] = get_sparseM(R);
[ao2_ra, ao2_ja, ao2_ia] = get_sparseM(ao2);
[e_ra, e_ja, e_ia] = get_sparseM(e);

% Fortran function call
[coef,ierr,maxiter,time] = srqfncGate(n1, m, nnzx, x_ra, x_ja, x_ia,...
    ao1_ra, ao1_ja, ao1_ia, n2, nnzR,R_ra, R_ja, R_ia, ao2_ra,       ...
    ao2_ja, ao2_ia, nnzdmax, zeros(1, nnzdmax), zeros(1, nnzdmax),   ...
    zeros(1,m+1), zeros(1,nnzhmax+1), zeros(1,nnzhmax+1), nnzemax,   ...
    e_ra, e_ja, e_ia, nnzgmax, zeros(1,nnzgmax), zeros(1,nnzgmax),   ...
    zeros(1,m+1), nnzhmax, zeros(1,nnzhmax), zeros(1,nnzhmax),       ...
    zeros(1,m+1), nsubmax, zeros(1,nsubmax), zeros(1,m+1), nnzlmax,  ...
    zeros(1,nnzlmax), zeros(1,m+1), zeros(1,m*5), iwmax,             ...
    zeros(1,iwmax), zeros(1,m+1), tmpmax, zeros(1,tmpmax), maxn1n2,  ...
    zeros(1,maxn1n2), wwm, wwn1, wwn2, cachsz, 8, x1, x2, s, u, y', r,...
    b', small, 0, maxiter, zeros(1,7));

%     ierr <- srqfnc.o$ierr
%     if (ierr == 13)
%         stop("Increase nnzh.factor")
if ierr == 13
    error('Increase nnzh factor');
end
%     if (!(ierr == 0) && warn.mesg)
%         warning(sfnMessage(ierr))
if ierr ~= 0
    warning(sfnMessage(ierr))
end
%     list(coef = -srqfnc.o$sol, ierr = ierr, it = srqfnc.o$maxiter,
%         time = sum(srqfnc.o$time))
coef = -coef;
end

% se.constr = function(x,y,coef,tau)
function cov = se_constr(x,y,coef,tau)
% {
%         require(quantreg)
%         n = nrow(x)
n = size(x,1);
%         p = ncol(x)
p = size(x,2);
%         h = bandwidth.rq(tau, n)
h = bandwidth_rq(tau,n);
%         if (tau + h > 1)
%             stop("tau + h > 1:  error in summary.rq")
if tau+h>1
    error("tau + h > 1: error in summary.rq")
end
%         if (tau - h < 0)
%             stop("tau - h < 0:  error in summary.rq")
if tau-h<0
    error("tau - h < 0: error in summary.rq")
end
%         uhat = c(y - x %*% coef)
uhat = y - (x*coef);
%         h = (qnorm(tau + h) - qnorm(tau - h)) * min(sqrt(var(uhat)),
%             (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
h = (norminv(tau + h) - norminv(tau - h)) * min([sqrt(var(uhat))        ...
    (r_quantile(uhat,0.75) - r_quantile(uhat, 0.25))./1.34]);
%         f = dnorm(uhat/h)/h
f = normpdf(uhat./h)/h;
%         fxxinv = diag(p)
fxxinv = diag(ones(1,p));
%         fxxinv = backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p, drop = FALSE],
%             fxxinv)
[~,R] = qr(sqrt(f).*x);
R = R(1:p,1:p);
fxxinv = R\fxxinv;
%         fxxinv = fxxinv %*% t(fxxinv)
fxxinv = fxxinv * fxxinv';
%         cov = tau * (1 - tau) * fxxinv %*% crossprod(x) %*% fxxinv
cov = tau .* (1-tau).*fxxinv*(x'*x)*fxxinv;
%         cov
% }
end
