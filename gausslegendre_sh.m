function [Mu_s,W_s] = gausslegendre_sh(Mu,W);

% This function provide the weights and ascissas for the integral between 0
% and 1
% -------------------------------------------------------------------------

% Assign a and b (Integral limits)
a = 0;
b = 1;

% Compute alpha and beta
Am = [1 a;1 b];
bm = [-1;1];
xm = Am\bm;
alpha = xm(1);
beta = xm(2);

% Evaluate M
siz = size(Mu);
M = siz(1);

% Initialize the matrices
Mu_s = zeros(M,1);
W_s = zeros(M,1);

% Compute Mu_s and W_s
Mu_s = (1/beta)*(Mu - alpha*ones(M,1));
W_s = (1/beta)*W;


    
    