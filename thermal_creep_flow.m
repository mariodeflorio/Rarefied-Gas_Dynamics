%%
clear; close all; clc;
format long
%--------------------------------------------------------------------------
%{
  1D thermal creep channel flow in the BGK approximation
  via X-TFC

  Authors:
  Mario De Florio, PhD 
  Enrico Schiassi, PhD
%}
%%
%--------------------------------------------------------------------------
%% Input

N = 24;    % Discretization order for mu (0,1]
M = 400;    % Discretization order for tau (and x)
m = 90;    % number of basis function
a = 1; % Semi-thicnkess of channel
alpha = 1; % Accomodation coefficient (0,1]

type_AF = 2; % activation function: 1 = Chebyshev NN; 2 = Legendre NN.
%--------------------------------------------------------------------------

%Parameters for simplicity
beta = (1 - alpha)/alpha;
theta = (1/alpha);

%%
%--------------------------------------------------------------------------
tau_0 = 0; tau_f = a; % tau span

tau = linspace(0,a,M)';

%--------------------------------------------------------------------------
%% PROBLEM FORMULATION
%--------------------------------------------------------------------------
%% Select weights

[Mu,W] = gausslegendre(N); % Both column vector

% integral between 0 and 1:
% Note: The weigths and abscissas are shifted
[mu,w] = gausslegendre_sh(Mu,W);

% % Re-order weights and bias in 0-1
for i = 1:N
    Mu_ss(i,1) = mu(N - i + 1);
    W_ss(i,1) = w(N - i + 1);
end

% change Back
mu = Mu_ss;
w = W_ss;

% Define the vector associated to the Psi
% Note; The Psi function is function of zeta, mapped in the interval [0,1]
psi = zeros(N,1);
g = zeros(N,1);
for i = 1:N
    
    % compute psi and f
    psi(i) = (1/sqrt(pi))*(1/mu(i))*exp(-((-log(mu(i)))^2));
    g(i) = (alpha/2)*((-log(mu(i)))^2 - 0.5);
end

%--------------------------------------------------------------------------

%% ToC Formulation
%% Definition new variable x (for basis expansion)
% To check the collocation method
x0 = -1; xf = 1; % x span (for orthogonal polynomials)
c = (xf - x0)/(tau_f - tau_0);
x = x0 + c*(tau - tau_0);
%--------------------------------------------------------------------------

%% Basis functions evaluation (Orthogonal polynomials)

switch type_AF
    case 1
        [H, HD] = CP(x, m + 1);
        
    case 2
        [H, HD] = LeP(x, m + 1);
end

% Restore matrices (to start with the 2nd order polynom)
I = 2:m+1;
H = H(:,I);  HD = HD(:,I);
H0=H(1,:); Hf=H(end,:);

%--------------------------------------------------------------------------
%% Building A*xi=B for the LS solution

%% Construction of 'A' matrix

% Construction of summation on k of the inhomogeneous term
sum_k = 0;
for k = 1:N
    sum_k = sum_k + w(k)*psi(k)*2*theta*g(k);
end

xi_pos =  zeros(N*m,1);
xi_neg =  zeros(N*m,1);

A = zeros(2*M*N,2*m*N);
B = zeros(2*M*N,1);

for i = 1:N
    
    for k = 1:N
        if k == i
            A_abcd_k(:,(2*m *(k-1)) + 1 :( 2*m * k ) ) = [ -c*log(mu(i))*HD + H + beta*Hf - ( 1 + beta)*H0 - w(k)*psi(k)*(H + 2*beta*Hf - ( 1 + 2*beta)*H0) ,   - theta*Hf + (1 + beta)*H0  - w(k)*psi(k)*(H - 2*theta*Hf + ( 1 + 2*beta)*H0) ;
                 beta*Hf - beta*H0  - w(k)*psi(k)*(H + 2*beta*Hf - (1 + 2*beta)*H0 ) ,  c*log(mu(i))*HD + H - theta*Hf + beta*H0 - w(k)*psi(k)*(H - 2*theta*Hf + (1 + 2*beta)*H0 ) ];
        else
            A_abcd_k(:,(2*m *(k-1)) + 1 :( 2*m * k ) ) = [ - w(k)*psi(k)*(H + 2*beta*Hf - ( 1 + 2*beta)*H0) , - w(k)*psi(k)*(H - 2*theta*Hf + (1+2*beta)*H0) ;
                 - w(k)*psi(k)*(H + 2*beta*Hf - ( 1 + 2*beta)*H0)  ,  - w(k)*psi(k)*(H - 2*theta*Hf + (1 + 2*beta)*H0 ) ];
        end
    end
    
    A ((2*M *(i-1)) + 1 :( 2*M * i ),: ) = A_abcd_k;
    
    b(1:2*M,1) = - theta*g(i) + sum_k;
    B ((2*M *(i-1)) + 1 :( 2*M * i ),: ) = b;
    
end

%--------------------------------------------------------------------------
%% Xi computation via Least-Squares
tic
xi = (A'*A)\(A'*B);
toc

% Computation of xi_pos and xi_neg vectors
for r = 0:(N-1)
    xi_pos((r*m)+1 : (r+1)*m) = xi(   2*r*m + 1 : (2*r+1)*m);
    xi_neg((r*m)+1 : (r+1)*m) = xi( (2*r+1)*m+1 : (2+2*r)*m);
end

% Computation of the solution

Y_p = zeros(M,N);
Y_n = zeros(M,N);
YD_p = zeros(M,N);
YD_n = zeros(M,N);

for i = 0:(N-1)
    
    Y_p(:,i+1) = ( H + beta*Hf - (1 + beta)*H0 ) * xi_pos((i*m)+1 : (i+1)*m)  +  ( (1 + beta)*H0 - theta*Hf) * xi_neg((i*m)+1 : (i+1)*m)  + theta*g(i+1) ;
    Y_n(:,i+1) =  beta*(Hf - H0)* xi_pos((i*m)+1 : (i+1)*m)  + (H + beta*H0 - theta*Hf ) * xi_neg((i*m)+1 : (i+1)*m) + theta*g(i+1) ;
    
    YD_p(:,i+1) = HD * xi_pos((i*m)+1 : (i+1)*m);
    YD_n(:,i+1) = HD * xi_neg((i*m)+1 : (i+1)*m);
    
end

%%
% Macroscopic velocity profile
Y_0 = 0;
for i = 1:N
    Y_0 = Y_0 + w(i)*psi(i)*(Y_p(:,i) + Y_n(:,i)) ;
end

q = Y_0;

figure(1)
hold on
grid on
xlabel('\tau')
ylabel('q')
plot(tau,q,'LineWidth',2);
title('Macroscopic velocity profile in half-channel')


%% Flow rate

MM = 1e4;

tau_int = linspace(0,a,MM)';

x0 = -1; xf = 1; % x span (for orthogonal polynomials)
c = (xf - x0)/(tau_f - tau_0);
x = x0 + c*(tau_int - tau_0);
%--------------------------------------------------------------------------

%% Activation functions evaluation 

switch type_AF
    case 1
        [h, hd] = CP(x, m + 1);
        
    case 2
        [h, hd] = LeP(x, m + 1);
end

% Restore matrices (to start with the 2nd order polynom)
I = 2:m+1;
h = h(:,I);  hd = hd(:,I);
h0=h(1,:); hf=h(end,:);

Y_p = zeros(MM,N);
Y_n = zeros(MM,N);
YD_p = zeros(MM,N);
YD_n = zeros(MM,N);

for i = 0:(N-1)
    
    Y_p(:,i+1) = ( h + beta*hf - (1 + beta)*h0 ) * xi_pos((i*m)+1 : (i+1)*m)  +  ( (1 + beta)*h0 - theta*hf) * xi_neg((i*m)+1 : (i+1)*m)  + theta*g(i+1) ;
    Y_n(:,i+1) =  beta*(hf - h0)* xi_pos((i*m)+1 : (i+1)*m)  + (h + beta*h0 - theta*hf ) * xi_neg((i*m)+1 : (i+1)*m) + theta*g(i+1) ;
    
    YD_p(:,i+1) = hd * xi_pos((i*m)+1 : (i+1)*m);
    YD_n(:,i+1) = hd * xi_neg((i*m)+1 : (i+1)*m);
    
end

%%
% Macroscopic velocity profile
Y_0 = 0;
for i = 1:N
    Y_0 = Y_0 + w(i)*psi(i)*(Y_p(:,i) + Y_n(:,i)) ;
end

% Flow rate

Q = -(1/(a^2)) * trapz(tau_int,Y_0);

fprintf('The flow rate Q is %.7f \n' , Q) 









