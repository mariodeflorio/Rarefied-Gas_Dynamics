%%
clear; close all; clc;
format long
%--------------------------------------------------------------------------
%{
  1D Poiseuille flow in the BGK approximation
  via X-TFC

  Authors:
  Mario De Florio, PhD 
  Enrico Schiassi, PhD
%}
%--------------------------------------------------------------------------
%% Input
start = tic;
N = 24;    % Discretization order for mu (0,1]
M = 400;    % Discretization order for tau (and x)
m = 90;    % number of basis function

a = 1; % Semi-thicnkess of channel
alpha = 1; % Accomodation coefficient (0,1]

type_AF = 1; % activation function: 1 = Chebyshev NN; 2 = Legendre NN.
%--------------------------------------------------------------------------

% Parameters for simplicity
phi = (1 - alpha)/alpha;
gamma = (1/alpha);

%%
%--------------------------------------------------------------------------
tau_0 = 0; tau_f = a; % tau span
eta_d = [0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1];
eta_d = a*eta_d;


tau = linspace(0,a,M)';


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
for i = 1:N
    
    % compute psi and f
    psi(i) = (1/sqrt(pi))*(1/mu(i))*exp(-((-log(mu(i)))^2));
    kappa(i) = alpha*(-log(mu(i)))^2 - a*log(mu(i))*(2 - alpha);
    
end

%--------------------------------------------------------------------------

%% X-TFC Formulation
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
    sum_k = sum_k + w(k)*psi(k)*2*gamma*kappa(k);
end

beta_pos =  zeros(N*m,1);
beta_neg =  zeros(N*m,1);

A = zeros(2*M*N,2*m*N);
B = zeros(2*M*N,1);

for i = 1:N
    
    for k = 1:N
        if k == i
            A_abcd_k(:,(2*m *(k-1)) + 1 :( 2*m * k ) ) = [ -c*log(mu(i))*HD + H + phi*Hf - ( 1 + phi)*H0 - w(k)*psi(k)*(H + 2*phi*Hf - ( 1 + 2*phi)*H0) ,   - gamma*Hf + (1 + phi)*H0  - w(k)*psi(k)*(H - 2*gamma*Hf + ( 1 + 2*phi)*H0) ;
                 phi*Hf - phi*H0  - w(k)*psi(k)*(H + 2*phi*Hf - (1 + 2*phi)*H0 ) ,  c*log(mu(i))*HD + H - gamma*Hf + phi*H0 - w(k)*psi(k)*(H - 2*gamma*Hf + (1 + 2*phi)*H0 ) ];
        else
            A_abcd_k(:,(2*m *(k-1)) + 1 :( 2*m * k ) ) = [ - w(k)*psi(k)*(H + 2*phi*Hf - ( 1 + 2*phi)*H0) , - w(k)*psi(k)*(H - 2*gamma*Hf + (1+2*phi)*H0) ;
                 - w(k)*psi(k)*(H + 2*phi*Hf - ( 1 + 2*phi)*H0)  ,  - w(k)*psi(k)*(H - 2*gamma*Hf + (1 + 2*phi)*H0 ) ];
        end
    end
    
    A ((2*M *(i-1)) + 1 :( 2*M * i ),: ) = A_abcd_k;
    
    b(1:2*M,1) = - gamma*kappa(i) + sum_k;
    B ((2*M *(i-1)) + 1 :( 2*M * i ),: ) = b;
    
end

%--------------------------------------------------------------------------
%% Xi computation via Least-Squares

startleast = tic;

beta = (A'*A)\(A'*B);

LEASTSQUARE = toc(startleast)

% Computation of xi_pos and xi_neg vectors
for r = 0:(N-1)
    beta_pos((r*m)+1 : (r+1)*m) = beta(   2*r*m + 1 : (2*r+1)*m);
    beta_neg((r*m)+1 : (r+1)*m) = beta( (2*r+1)*m+1 : (2+2*r)*m);
end

% Computation of the solution

Y_p = zeros(M,N);
Y_n = zeros(M,N);
YD_p = zeros(M,N);
YD_n = zeros(M,N);

for i = 0:(N-1)
    
    Y_p(:,i+1) = ( H + phi*Hf - (1 + phi)*H0 ) * beta_pos((i*m)+1 : (i+1)*m)  +  ( (1 + phi)*H0 - gamma*Hf) * beta_neg((i*m)+1 : (i+1)*m)  + gamma*kappa(i+1) ;
    Y_n(:,i+1) =  phi*(Hf - H0)* beta_pos((i*m)+1 : (i+1)*m)  + (H + phi*H0 - gamma*Hf ) * beta_neg((i*m)+1 : (i+1)*m) + gamma*kappa(i+1) ;
    
    YD_p(:,i+1) = HD * beta_pos((i*m)+1 : (i+1)*m);
    YD_n(:,i+1) = HD * beta_neg((i*m)+1 : (i+1)*m);
    
end

%%

Y_0 = 0;
for i = 1:N
    Y_0 = Y_0 + w(i)*psi(i)*(Y_p(:,i) + Y_n(:,i)) ;
end

% Macroscopic velocity profile
q = zeros(M,1);
for i = 1:length(tau)
     q(i) = 0.5*(1 - a^2 + (tau(i))^2) - Y_0(i);
end


figure(1)
hold on
grid on
xlabel('\tau')
ylabel('q')
plot(tau,q,'LineWidth',2);
title('Macroscopic velocity profile in half-channel')

% Loss
Lossp = zeros(M,N);
Lossn = zeros(M,N);

for i=1:N
    
    sum_loss_k = 0;
    
    for k = 1:N
        sum_loss_k = sum_loss_k + w(k)*psi(k)*(Y_p(:,k) + Y_n(:,k));
    end
    
    Lossp(:,i) =  -c*log(mu(i))*YD_p(:,i) + Y_p(:,i) - sum_loss_k;
    Lossn(:,i) =  c*log(mu(i))*YD_n(:,i) + Y_n(:,i) - sum_loss_k;
    
 end


figure(2)
hold on
grid on
xlabel('\tau')
ylabel('Loss')
set(gca, 'YScale','log')
plot(tau,mean(abs(Lossp),2),'*','Color','r');
title(sprintf('Loss function - positive flux')) ;

%
figure(3)
hold on
grid on
xlabel('\tau')
ylabel('Loss')
set(gca,'YScale','log')
plot(tau,mean(abs(Lossn),2),'*');
title(sprintf('Loss function - negative flux')) ;


%% Flow rate

funInt = griddedInterpolant(tau,q,'spline');
fun = @(tau) funInt(tau);

int = integral(fun,0,a);

Q = 2*(-1/(2*a^2))*int ;
fprintf('The flow rate Q is %.7f \n' , Q) 

