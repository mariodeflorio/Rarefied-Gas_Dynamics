function [X,W] = gausslegendre(order)

% Generate the abscissa and weights
vector = (1:order-1)./sqrt((2*(1:order-1)).^2-1);
[w1,xi1] = eig(diag(vector,-1)+diag(vector,1));
xi1 = diag(xi1);
w1 = w1(1,:)'.^2;


X1 = [xi1];
W1 = 2*[w1(:).']'; %weights normalize such that Sum_n wn = 2
%W1 = [w1(:).']'; %weights normalize such that Sum_n wn = 1

% Adapting to the integration that starts with the positive mu
k = 1;
% initialize
W = zeros(order,1);
X = zeros(order,1);

for i = length(W1):-1:1

   W(k) = W1(i);
   X(k) = X1(i);
    
   k = k + 1;
 
end

